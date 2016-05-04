#! /usr/bin/env python

"""Solve some tasks with A* and the LM-Cut heuristic."""

import os.path
import platform

from lab.environments import LocalEnvironment, MaiaEnvironment

from downward.experiment import FastDownwardExperiment
from downward.reports.absolute import AbsoluteReport


SUITE = ['gripper:prob01.pddl']
ATTRIBUTES = ['coverage', 'expansions']

if 'cluster' in platform.node():
    REPO = os.path.expanduser('~/projects/downward')
    BENCHMARKS_DIR = os.path.expanduser('~/projects/benchmarks')
    ENV = MaiaEnvironment(priority=-10)
else:
    REPO = os.path.expanduser('~/projects/Downward/downward')
    BENCHMARKS_DIR = os.path.expanduser('~/projects/Downward/benchmarks')
    ENV = LocalEnvironment(processes=2)
CACHE_DIR = os.path.expanduser('~/lab')

exp = FastDownwardExperiment(environment=ENV, cache_dir=CACHE_DIR)
exp.add_suite(BENCHMARKS_DIR, SUITE)
exp.add_algorithm(
    'lmcut', REPO, 'tip', ['--search', 'astar(lmcut())'])

# Make a report (AbsoluteReport is the standard report).
exp.add_report(
    AbsoluteReport(attributes=ATTRIBUTES), outfile='report.html')

# Parse the commandline and show or run experiment steps.
exp()
