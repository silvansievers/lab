#! /usr/bin/env python

import os

import standard_exp


EXPPATH = os.path.join(standard_exp.EXPS, 'preprocess-all')
SUITE = [
    'appn-adl', 'appn-strips', 'briefcaseworld',
    'fsc-blocks', 'fsc-grid-a1', 'fsc-grid-a2', 'fsc-grid-r', 'fsc-hall', 'fsc-visualmarker',
    'ged2-ds1', 'ged2-ds2nd', 'ged3-ds1', 'ged3-ds2nd', 'gedp-ds2ndp',
    't0-adder', 't0-coins', 't0-comm', 't0-grid-dispose', 't0-grid-push',
    't0-grid-trash', 't0-sortnet', 't0-sortnet-alt', 't0-uts'
]

exp = standard_exp.StandardDownwardExperiment(path=EXPPATH)
exp.add_suite(SUITE)
exp.add_config('unused', ['unused'])
del exp.steps[3:]

# Parse the commandline and show or run experiment steps.
exp()
