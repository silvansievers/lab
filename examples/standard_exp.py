#! /usr/bin/env python

import os
import platform
import shutil
from subprocess import call
import sys

from lab.environments import LocalEnvironment, MaiaEnvironment
from lab.steps import Step

from downward.experiment import DownwardExperiment
from downward.reports.absolute import AbsoluteReport


NODE = platform.node()
REMOTE = 'cluster' in NODE or NODE.startswith('gkigrid') or NODE in ['habakuk', 'turtur']
SCP_LOGIN = 'maia'
ATTRIBUTES = ['coverage', 'cost', 'quality', 'total_time']

REMOTE_EXPS = '/infai/sieverss/experiments/executed'
LOCAL_EXPS = '/home/silvan/experiments/executed'

REMOTE_REPO = '/infai/sieverss/repos/fd-master-clone-rev-08f33fca24e6'
LOCAL_REPO = '/home/silvan/work/downward/master'

REMOTE_PYTHON = 'python2.7'
LOCAL_PYTHON = 'python2.7'

if REMOTE:
    EXPS = REMOTE_EXPS
    REPO = REMOTE_REPO
    PYTHON = REMOTE_PYTHON
    CACHE_DIR = os.path.expanduser('~/experiments/cache_dir')
else:
    EXPS = LOCAL_EXPS
    REPO = LOCAL_REPO
    PYTHON = LOCAL_PYTHON
    CACHE_DIR = os.path.expanduser('~/experiments/cache_dir')


class StandardDownwardExperiment(DownwardExperiment):
    def __init__(self, path=None, repo=None, environment=None,
                 combinations=None, limits=None, attributes=None,
                 derived_properties=None, priority=0, queue=None,
                 processes=2, email=None, cache_dir=CACHE_DIR, **kwargs):
        if path is None:
            path = os.path.splitext(os.path.basename(sys.argv[0]))[0]

        expname = os.path.basename(path)

        remote_exppath = os.path.join(REMOTE_EXPS, path)
        local_exppath = os.path.join(LOCAL_EXPS, path)

        if REMOTE:
            exppath = remote_exppath
            repo = repo or REMOTE_REPO
            environment = environment or MaiaEnvironment(priority=priority,
                                                         queue=queue,
                                                         email=email)
        else:
            exppath = local_exppath
            repo = repo or LOCAL_REPO
            environment = environment or LocalEnvironment(processes=processes)

        DownwardExperiment.__init__(self, path=exppath, environment=environment,
                                    repo=repo, combinations=combinations,
                                    limits=limits, cache_dir=cache_dir, **kwargs)

        self.set_path_to_python(PYTHON)

        if attributes is None:
            attributes = ATTRIBUTES

        # Add report steps
        abs_report_file = os.path.join(self.eval_dir, '%s-abs.html' % expname)
        self.add_report(AbsoluteReport(attributes=attributes, colored=True, derived_properties=derived_properties),
                        name='report-abs', outfile=abs_report_file)

        self.add_step(Step('remove-eval-dir', shutil.rmtree, self.eval_dir, ignore_errors=True))

        if not REMOTE:
            # Copy the results to local directory
            self.add_step(Step('scp-eval-dir', call, [
                'scp', '-r',
                '%s:%s-eval' % (SCP_LOGIN, remote_exppath),
                '%s-eval' % local_exppath]))


def get_exp(suite, configs, combinations=None, limits=None, attributes=None):
    exp = StandardDownwardExperiment(combinations=combinations, limits=limits,
                                     attributes=attributes)

    exp.add_suite(suite)
    for nick, config in configs:
        exp.add_config(nick, config)
    return exp
