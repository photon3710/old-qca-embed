#!/usr/bin/python

from __future__ import division
import pylab as plt
import os
import re

DENSE_DIR = '../bench/dense/'
HEUR_DIR = '../bench/heur/'

SHOW_FULL = True

DENSE_DIR += 'full/' if SHOW_FULL else 'lim/'
HEUR_DIR += 'full/' if SHOW_FULL else 'lim/'

IMG_DIR = '../../img/'

FS = 14
DY = .1


def pinch(st, l, r):
    return st.partition(l)[2].rpartition(r)[0]


def strToQbit(st):
    if 'None' in st:
        return None
    else:
        try:
            return tuple(map(int, pinch(st, '(', ')').split(',')))
        except:
            return None


def loadDense(fn):

    try:
        fp = open(fn, 'r')
    except IOError:
        print('Failed to open file: %s' % fn)
        return None

    qb_flag = False
    pt_flag = False
    hd_flag = False
    dq_flag = False

    qbits = {}
    paths = {}
    dis_qbits = []

    for line in fp:
        if len(line) < 3:
            continue
        if '<header>' in line:
            hd_flag = True
            continue
        elif '<dis_qbits>' in line:
            dq_flag = True
            continue
        elif '<qubits>' in line:
            qb_flag = True
            continue
        elif '<paths>' in line:
            pt_flag = True
            continue
        elif '</header>' in line:
            hd_flag = False
            continue
        elif '</dis_qbits>' in line:
            dq_flag = False
            continue
        elif '</qubits>' in line:
            qb_flag = False
            continue
        elif '</paths>' in line:
            pt_flag = False
            continue

        # parse header
        if hd_flag:
            if 'M' in line:
                M = int(line.split('=')[1])
            elif 'N' in line:
                N = int(line.split('=')[1])
            elif 'L' in line:
                L = int(line.split('=')[1])

        # disabled qubits
        if dq_flag:
            q = strToQbit(line.strip())
            dis_qbits.append(q)

        # parse qubits
        if qb_flag:
            cell, qbit = line.strip().split(':')
            cell = int(cell)
            qbit = strToQbit(qbit)
            qbits[cell] = qbit

        # parse paths
        elif pt_flag:
            key, path_st = line.strip().split(':')
            key = tuple(map(int, key.split(';')))
            path = map(lambda x: strToQbit(x), path_st.split(';'))
            paths[key] = path

    return qbits, paths, dis_qbits, M, N, L


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

    return models


def loadDenseSols(fnames):
    ''' '''

    out = [0, []]

    for fname in fnames:

        sol = loadDense(fname)
        if sol is None:
            continue
        qbits, paths, dis_abits, M, N, L = sol
        chain_lengths = map(lambda x: len(x)-2, paths.values())
        num_qbits = len(qbits)+sum(chain_lengths)
        out[1].append(num_qbits)

    # compute mean qubit count
    if len(out[1]) > 0:
        mean_qbits = sum(out[1])/len(out[1])
        out[0] = mean_qbits
    else:
        return None

    return out


def loadHeurSols(fnames):
    ''' '''

    out = [0, []]

    for fname in fnames:

        models = loadHeur(fname)
        if models is None:
            continue
        num_qbits = sum(map(len, models.values()))
        out[1].append(num_qbits)

    # compute mean qubit count
    if len(out[1]) > 0:
        mean_qbits = sum(out[1])/len(out[1])
        out[0] = mean_qbits
    else:
        return None

    return out


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
        fnames = filter(rgx.match, os.listdir(DENSE_DIR + key))
        fnames = map(lambda s: DENSE_DIR + key + '/' + s, fnames)
        sol = loadDenseSols(fnames)
        dense_sols[key] = sol
        if not sol is None:
            print '%.2f \t %d' % (sol[0], len(sol[1]))
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
            print '%.2f \t %d' % (sol[0], len(sol[1]))
        else:
            print 'None'

    ## show results

    keys = filter(lambda k: not dense_sols[k] is None, dense_sols.keys())
    key_order = sorted(keys, key=lambda k: dense_sols[k][0])

    axes = plt.figure().add_subplot(111)

    i = 0

    for key in key_order:
        i += 1
        ds = dense_sols[key]
        hs = heur_sols[key]

        # mean runtimes
        if not ds is None:
            axes.plot(ds[0], i+DY, 'bx', markersize=10, markeredgewidth=2)
        if not hs is None:
            axes.plot(hs[0], i-DY, 'go', markersize=10, markeredgewidth=2)

        # trials
        if not ds is None:
            axes.plot(ds[1], [i+DY]*len(ds[1]), 'b.')
        if not hs is None:
            axes.plot(hs[1], [i-DY]*len(ds[1]), 'g.')

    y_labels = [''] + key_order + ['']

    plt.legend(['Dense Placement', 'Heuristic'], fontsize=FS,
               loc='lower right', numpoints=1)
    axes.set_yticklabels(y_labels)
    axes.tick_params(axis='x', which='major', labelsize=FS)
    plt.xlabel('Qubits Used', fontsize=FS)
    #plt.title('Qubit Usage', fontsize=FS)
    sfname = 'bench-' + ('full' if SHOW_FULL else 'lim') + '.eps'
    plt.savefig(IMG_DIR+sfname, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
