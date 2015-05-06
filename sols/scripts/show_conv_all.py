#!/usr/bin/python

import numpy as np
from math import ceil
import os
import re


DENSE_DIR = '../bench/dense/'
HEUR_DIR = '../bench/heur/'

SHOW_FULL = True

DENSE_DIR += 'full/' if SHOW_FULL else 'lim/'
HEUR_DIR += 'full/' if SHOW_FULL else 'lim/'


def pinch(st, l, r):
    return st.partition(l)[2].rpartition(r)[0]


def loadDense(fname):
    ''' '''

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open file: %s' % fname
        return None

    max_model = float(fp.readline().split(':')[1])

    fp.close()

    return int(ceil(max_model))


def loadHeur(fn):
    ''' '''
    try:
        fp = open(fn, 'r')
    except IOError:
        print 'Failed to open file: ' + fn + '...'
        return None

    models = {}

    for line in fp:
        if '[' in line:    # new vertex-model
            cell, model = line.strip().split(':')
            cell = int(cell)
            qbit_list = pinch(model, '[', ']')
            model = [int(qbit) for qbit in qbit_list.split(',')]
            models[cell] = model

    fp.close()

    max_model = max(map(len, models.values()))

    return max_model


def loadDenseSols(fnames):
    ''' '''

    max_models = []
    for fname in fnames:
        max_model = loadDense(fname)
        if not max_model is None:
            max_models.append(max_model)

    if not max_models:
        return None
    return np.mean(max_models)


def loadHeurSols(fnames):
    ''' '''

    max_models = []
    for fname in fnames:
        max_model = loadHeur(fname)
        if not max_model is None:
            max_models.append(max_model)

    if not max_models:
        return None

    return np.mean(max_models)


def main():
    ''' '''

    # set up file structure
    dense_dirs = next(os.walk(DENSE_DIR))[1]

    heur_dirs = next(os.walk(HEUR_DIR))[1]

    dense_sols = {key: [] for key in dense_dirs}
    heur_sols = {key: [] for key in heur_dirs}

    print '*'*40
    print 'Reading Dense Solutions...\n'
    # get dense solutions
    rgx = re.compile('^sol[0-9]+$')
    for key in dense_sols:
        print '%-15s' % key,
        root = DENSE_DIR + key + '/conv_dir/'
        fnames = filter(rgx.match, os.listdir(root))
        fnames = map(lambda s: root + '/' + s, fnames)
        sol = loadDenseSols(fnames)
        dense_sols[key] = sol
        if not sol is None:
            print '%.2f' % sol
        else:
            print 'None'

    print '\n\n' + '*'*40
    print 'Reading Heuristic Solutions...\n'
    # get heur solutions
    rgx = re.compile('^solution[0-9]+$')
    for key in heur_sols:
        print '%-15s' % key,
        fnames = filter(rgx.match, os.listdir(HEUR_DIR + key))
        fnames = map(lambda s: HEUR_DIR + key + '/' + s, fnames)
        sol = loadHeurSols(fnames)
        heur_sols[key] = sol
        if not sol is None:
            print '%.2f' % sol
        else:
            print 'None'


if __name__ == '__main__':
    main()
