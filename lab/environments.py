# -*- coding: utf-8 -*-
#
# lab is a Python API for running and evaluating algorithms.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import datetime
import logging
import math
import multiprocessing
import os
import pkgutil
import random
import sys

from lab import tools
from lab.steps import Sequence


class Environment(object):
    def __init__(self):
        self.exp = None
        self.main_script_file = 'run'

    def write_main_script(self):
        raise NotImplementedError

    def build_linked_resources(self, run):
        """
        Only if we are building an argo experiment, we need to add all linked
        resources to the resources list.
        """
        pass

    def start_exp(self):
        """
        Run main experiment step.
        """
        raise NotImplementedError

    def run_steps(self):
        raise NotImplementedError


class LocalEnvironment(Environment):
    """
    Environment for running experiments locally on a single machine.
    """
    def __init__(self, processes=None):
        """
        If given, *processes* must be in the range [1, ..., #CPUs].
        If omitted, it will be set to #CPUs.
        """
        Environment.__init__(self)
        cores = multiprocessing.cpu_count()
        if processes is None:
            processes = cores
        if not 1 <= processes <= cores:
            raise ValueError("processes must be in the range [1, ..., #CPUs].")
        self.processes = processes

    def write_main_script(self):
        dirs = [repr(os.path.relpath(run.path, self.exp.path))
                for run in self.exp.runs]
        replacements = {'DIRS': ',\n'.join(dirs),
                        'PROCESSES': str(self.processes)}

        script = pkgutil.get_data('lab', 'data/local-job-template.py')
        for orig, new in replacements.items():
            script = script.replace('"""' + orig + '"""', new)

        self.exp.add_new_file('', self.main_script_file, script, permissions=0o755)

    def start_exp(self):
        tools.run_command(['./' + self.main_script_file], cwd=self.exp.path)

    def run_steps(self, steps):
        Sequence.run_steps(steps)


class OracleGridEngineEnvironment(Environment):
    """Abstract base class for grid environments."""

    DEFAULT_QUEUE = None             # must be overridden in derived classes
    TEMPLATE_FILE = 'grid-job-header-template'  # can be overridden in derived classes
    MAX_TASKS = float('inf')         # can be overridden in derived classes
    DEFAULT_PRIORITY = 0             # can be overridden in derived classes
    HOST_RESTRICTIONS = {}           # can be overridden in derived classes
    DEFAULT_HOST_RESTRICTION = ""    # can be overridden in derived classes

    def __init__(self, queue=None, priority=None, host_restriction=None,
                 email=None, randomize_task_order=True, extra_options=None):
        """
        .. note::

            Previously, you had to provide the ``--all`` option on the
            commandline to queue steps sequentially on the grid engine.
            Now, the selected steps will be queued sequentially if at
            least one of them itself submits runs to the queue.

            For correct sequential execution, this class writes job
            files to <cache_dir>/grid-steps/<timestamp>-<exp-name> and
            makes them depend on one another. The driver.{log,err}
            files in this directory can be inspected if something goes
            wrong. Since the job files call the main experiment script
            during execution, it mustn't be changed during the
            experiment.

        *queue* must be a valid queue name on the grid.

        *priority* must be in the range [-1023, 0] where 0 is the highest
        priority. If you're a superuser the value can be in the range
        [-1023, 1024].

        If *email* is provided and the steps run on the grid, a message
        will be sent when the experiment finishes.

        If *randomize_task_order* is True (this is the default since
        version 1.5), tasks for runs are started in a random order.
        This is useful to avoid systematic noise due to e.g. one of
        the configs being run on a machine with heavy load. Note
        that due to the randomization, run directories may be
        pristine while the experiment is running even though the
        grid engine marks the runs as finished.

        Use *extra_options* to pass additional options. Example that
        allocates 16 cores per run on maia::

            MaiaEnvironment(extra_options='#$ -pe smp 16')

        """
        Environment.__init__(self)
        if queue is None:
            queue = self.DEFAULT_QUEUE
        if priority is None:
            priority = self.DEFAULT_PRIORITY
        if host_restriction is None:
            host_restriction = self.DEFAULT_HOST_RESTRICTION

        self.queue = queue
        self.host_spec = self._get_host_spec(host_restriction)
        assert priority in xrange(-1023, 1024 + 1)
        self.priority = priority
        self.runs_per_task = 1
        self.email = email
        self.randomize_task_order = randomize_task_order
        self.extra_options = extra_options or '## (not used)'

        # When submitting an experiment job, wait for this job name.
        self.__wait_for_job_name = None
        self._job_name = None

    @classmethod
    def _escape_job_name(cls, name):
        if name[0].isdigit():
            name = 'j' + name
        return name

    def _get_common_job_params(self):
        return {
            'logfile': 'driver.log',
            'errfile': 'driver.err',
            'priority': self.priority,
            'queue': self.queue,
            'host_spec': self.host_spec,
            'notification': '#$ -m n',
            'extra_options': self.extra_options,
        }

    def write_main_script(self):
        num_tasks = int(math.ceil(len(self.exp.runs) / float(self.runs_per_task)))
        if num_tasks > self.MAX_TASKS:
            logging.critical('You are trying to submit a job with %d tasks, '
                             'but only %d are allowed.' %
                             (num_tasks, self.MAX_TASKS))
        job_params = self._get_common_job_params()
        job_params.update(name=self._escape_job_name(self.exp.name),
                          num_tasks=num_tasks)
        header = pkgutil.get_data('lab', 'data/' + self.TEMPLATE_FILE) % job_params

        body_params = dict(num_tasks=num_tasks, run_ids='')
        if self.randomize_task_order:
            run_ids = [str(i + 1) for i in xrange(num_tasks)]
            random.shuffle(run_ids)
            body_params['run_ids'] = ' '.join(run_ids)
        body = pkgutil.get_data('lab', 'data/grid-job-body-template') % body_params

        filename = self.exp._get_abs_path(self.main_script_file)
        with open(filename, 'w') as file:
            logging.debug('Writing file "%s"' % filename)
            file.write('%s\n\n%s' % (header, body))

    def start_exp(self):
        submitted_file = os.path.join(self.exp.path, 'submitted')
        if os.path.exists(submitted_file):
            tools.confirm('The file "%s" already exists so it seems the '
                          'experiment has already been submitted. Are you '
                          'sure you want to submit it again?' % submitted_file)
        submit = ['qsub']
        if self.__wait_for_job_name:
            submit.extend(['-hold_jid', self.__wait_for_job_name])
        if self._job_name:
            # The name set in the job file will be ignored.
            submit.extend(['-N', self._job_name])
        submit.append(self.main_script_file)
        tools.run_command(submit, cwd=self.exp.path)
        # Write "submitted" file.
        with open(submitted_file, 'w') as f:
            f.write('This file is created when the experiment is submitted to '
                    'the queue.')

    def _get_script_args(self):
        """
        Retrieve additional commandline parameters given when the experiment
        is called by the user and pass them again when the step is called by
        the grid.
        """
        # Remove step names from the back of the commandline to avoid deleting
        # custom args by accident.
        commandline = list(reversed(sys.argv[1:]))
        if '--all' in commandline:
            commandline.remove('--all')
        for step_name in self.exp.args.steps:
            commandline.remove(step_name)
        return list(reversed(commandline))

    def _get_job_name(self, step):
        return self._escape_job_name(
            '%s-%02d-%s' % (self.exp.name, self.exp.steps.index(step) + 1,
                            step.name))

    def _get_job_header(self, step):
        job_params = self._get_common_job_params()
        job_params.update(name=self._get_job_name(step), num_tasks=1)
        if step.is_last_step and self.email:
            job_params['notification'] = '#$ -M %s\n#$ -m e' % self.email
        return pkgutil.get_data('lab', 'data/' + self.TEMPLATE_FILE) % job_params

    def _get_job(self, step):
        # Abort if one step fails.
        template = """\
%(job_header)s
if [ -s "%(stderr)s" ]; then
    echo "There was output on stderr. Please check %(stderr)s. Aborting."
    exit 1
fi

cd %(cwd)s
%(python)s %(script)s %(args)s %(step_name)s
"""
        return template % {
            'cwd': os.getcwd(),
            'python': sys.executable or 'python',
            'script': sys.argv[0],
            'args': ' '.join(self._get_script_args()),
            'step_name': step.name,
            'stderr': 'driver.err',
            'job_header': self._get_job_header(step)}

    def _get_host_spec(self, host_restriction):
        if not host_restriction:
            return '## (not used)'
        else:
            hosts = self.HOST_RESTRICTIONS[host_restriction]
            return '#$ -l hostname="%s"' % '|'.join(hosts)

    def run_steps(self, steps):
        timestamp = datetime.datetime.now().isoformat()
        job_dir = os.path.join(self.exp.cache_dir,
                               'grid-steps',
                               timestamp + '-' + self.exp.name)
        tools.overwrite_dir(job_dir)

        # Build the job files before submitting the other jobs.
        logging.info('Building job scripts')
        for step in steps:
            if step._funcname == 'build':
                script_step = step.copy()
                script_step.kwargs['only_main_script'] = True
                script_step()

        prev_job_name = None
        for number, step in enumerate(steps, start=1):
            job_name = self._get_job_name(step)
            # We cannot submit a job from within the grid, so we submit it
            # directly.
            if step._funcname == 'run':
                self.__wait_for_job_name = prev_job_name
                self._job_name = job_name
                step()
            else:
                step.is_last_step = (number == len(steps))
                with open(os.path.join(job_dir, job_name), 'w') as f:
                    f.write(self._get_job(step))
                submit = ['qsub']
                if prev_job_name:
                    submit.extend(['-hold_jid', prev_job_name])
                submit.append(job_name)
                tools.run_command(submit, cwd=job_dir)
            prev_job_name = job_name


class GkiGridEnvironment(OracleGridEngineEnvironment):
    """Environment for Freiburg's AI group."""

    DEFAULT_QUEUE = 'opteron_core.q'
    MAX_TASKS = 75000

    def run_steps(self, steps):
        if 'xeon' in self.queue:
            logging.critical('Experiments must be run stepwise on xeon, '
                             'because mercurial is missing there.')
        OracleGridEngineEnvironment.run_steps(self, steps)


def _host_range(prefix, from_num, to_num):
    return ['%s%02d*' % (prefix, num) for num in xrange(from_num, to_num + 1)]


class MaiaEnvironment(OracleGridEngineEnvironment):
    """Environment for Basel's AI group."""

    DEFAULT_QUEUE = '"all.q@ase*"'
    DEFAULT_HOST_RESTRICTION = ''

    # Note: the hosts in the following host restrictions are part of the
    # queue 'all.q' and not part of the default queue '"all.q@ase*"'.
    # Use them like this:
    # MaiaEnvironment(queue='all.q', host_restrictions='maia-six')
    HOST_RESTRICTIONS = {
        'maia-quad': _host_range('uni', 1, 32) + _host_range('ugi', 1, 8),
        'maia-six': _host_range('uni', 33, 72),
    }
