#---------------------------------------------------------
# Name: gen_test.py
# Purpose: Test code for circuit generator
# Author:    Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

from __future__ import division

import os
import re
import numpy as np
import pylab as plt
from random import random

from dwave_sapi import find_embedding, get_chimera_adjacency
from dense_placement.embed import denseEmbed, setChimeraSize, \
    getCouplerFlags, setQbitAdj

from generator import generateCircuit2, GtoCoef
from auxil import coefToConn
from time import time

import sys

FLAG_SOL = False
ONCE_FLAG = False    # continue after a single found solution

SHOW = False

NUM_RUNS = 10
NUM_TRIALS = 10
TIMEOUT = 300    # seconds

MN = 14     # number of rows and columns of tiles
M, N, L = MN, MN, 4
A_size = M*N*L*2
ndis = 0

SIZE_RANGE = {'full': [0, 400],
              'lim': [0, 500]}

WRITE_DIR = '../sols/%dgen/%d/' % (A_size, ndis)


def getFname(direc):
    '''Get unused filename for writing'''

    # check for directory existence
    if not os.path.exists(direc):
        os.makedirs(direc)

    # check for free filename
    regex = re.compile('^out[0-9]+$')
    old_outs = filter(regex.match, os.listdir(direc))
    old_ext = map(lambda x: int(x[3::]), old_outs)
    old_max = max(old_ext) if len(old_ext) > 0 else -1
    fname = direc + 'out%d' % (old_max+1)
    return fname


def writeToFile(OUT, append=False, typ='dense'):
    '''write generated circuit results to file'''

    write_op = 'a' if append else 'w'

    write_dir = WRITE_DIR + typ + '/'

    try:
        fname = getFname(write_dir)
        fp = open(fname, write_op)
    except:
        print 'Failed to open file: %s' % fname
        try:
            fname = getFname('.')
            print 'Saving to default: %s' % fname
            fp = open(fname, write_op)
        except:
            print 'Failed...'
            return None

    print 'Writing to: {0}'.format(fname)
    for key in OUT.keys():
        fp.write('<%s\>\n\n' % str(key))
        for d in OUT[key]:
            fp.write('%d\t%.2f\t%d\t%.5f\t%.3f\n' %
                     (d[0], d[1], d[2], d[3], d[4]))
        fp.write('\n\n')

    fp.close()


def procDenseEmbed(embed):
    '''Process formatted dense placement solution'''

    cell_map, paths, t = embed

    chain_lengths = map(lambda x: len(x)-2, paths.values())

    qubits = len(cell_map) + sum(chain_lengths)

    return qubits, chain_lengths, t


def procHeurEmbed(embed):
    '''Process an embedding from the heuristic algorithm'''
    global A_size

    chain_lengths = []
    flags = [False]*A_size
    fail = False

    embed, t = embed
    for chain in embed:
        for qbit in chain:
            if flags[qbit]:
                fail = True
                break
            flags[qbit] = True
        if fail:
            break
        chain_lengths.append(len(chain))

    if fail:
        return -1, chain_lengths, t
    else:
        return sum(chain_lengths), chain_lengths, t


def formatCoef(h, J):
    '''Format coefficient matrices for input to heuristic algorithm'''

    S = {}
    N = len(h)

    for i in xrange(N):
        for j in xrange(N):
            if abs(J[i, j]) > 0.01:
                S[(i, j)] = 1
            else:
                S[(i, j)] = 0

    return S, N


def runDense(h, J, max_count):
    '''Run the dense placement embedding algorithm'''

    source = coefToConn(h, J)
    count, num_success, num_trials = 0, 0, 0

    good_embeds = []

    while count < max_count:
        num_trials += 1
        t = time()
        try:
            cell_map, paths = denseEmbed(source, write=False)
            good_embeds.append([cell_map, paths, time()-t])
            num_success += 1
            print 'Solution %d found...' % num_success
            if ONCE_FLAG:
                break
        except Exception as e:
            if type(e).__name__ == 'KeyboardInterrupt':
                raise KeyboardInterrupt
            print '*'
            pass
        count = num_success if FLAG_SOL else num_trials

    return good_embeds


def runHeuristic(S, S_size, max_count):
    '''run the heuristic embedding algorithm'''

    global M, N, L, A_size

    A = get_chimera_adjacency(M, N, L)

    trial_num = 0
    success_num = 0
    count = 0

    good_embeds = []

    while count < max_count:

        t = time()
        embeddings = find_embedding(S, S_size, A, A_size,
                                    verbose=0, tries=1, timeout=TIMEOUT)
        trial_num += 1

        if len(embeddings) == S_size:    # successful embedding
            success_num += 1
            good_embeds.append([embeddings, time()-t])
            print 'solution ' + str(success_num) + ' found...'
            if ONCE_FLAG:
                break
        else:
            print '*'

        count = success_num if FLAG_SOL else trial_num

    print '\n\nEmbedded %d of %d attempts \n\n' % (success_num, trial_num)

    return good_embeds


def main(typ):
    ''' '''

    global M, N, L, ndis

    dis_coup = []
    dis_qbits = list(set(map(int, [random()*M*N*L*2 for _ in xrange(ndis)])))

    setChimeraSize(M, N, L)
    CF = getCouplerFlags(dis_coup, dis_qbits)
    setQbitAdj(CF)

    OUT = {}

    try:
        for full_adj in [True, False]:

            if full_adj:
                OUT['full'] = []
                out = OUT['full']
                size_range = SIZE_RANGE['full']
            else:
                OUT['lim'] = []
                out = OUT['lim']
                size_range = SIZE_RANGE['lim']

            print ('\n'+'*'*50)*2 + '\n** ',
            print 'FULL ADJACENCY' if full_adj else 'LIM ADJACENCY'

            for run in xrange(NUM_RUNS):

                print '\n' + '*'*40 + '\nRun %d\n\n' % (run+1)
                while True:
                    G = generateCircuit2(full_adj=full_adj)
                    if len(G) >= size_range[0] and len(G) <= size_range[1]:
                        break
                h, J = GtoCoef(G)
                if typ == 'dense':
                    good_embeds = runDense(h, J, max_count=NUM_TRIALS)
                else:
                    S, S_Size = formatCoef(h, J)
                    good_embeds = runHeuristic(S, S_Size, max_count=NUM_TRIALS)

                QUBITS = []
                TIMES = []

                if good_embeds:
                    for embed in good_embeds:
                        if typ == 'dense':
                            qubits, lengths, t = procDenseEmbed(embed)
                        else:
                            qubits, lengths, t = procHeurEmbed(embed)
                        QUBITS.append(qubits)
                        TIMES.append(t)
                else:
                    QUBITS = []

                if QUBITS:
                    mean_qubits = np.mean(QUBITS)
                    min_qubits = max(0, np.min(QUBITS))
                    mean_time = np.mean(TIMES)
                    p = len(good_embeds)*1./NUM_TRIALS
                    out.append([len(h), mean_qubits, min_qubits, mean_time, p])
                else:
                    out.append([len(h), -1, -1, -1, 0.])

            if SHOW:
                # show mean qubits
                c = ['g', 'r'] if full_adj else ['b', 'm']
                X, Y = [], []
                for d in out:
                    if d[2] == -1:
                        continue
                        plt.plot(d[0], 0, c[1]+'x',
                                 markersize=5, markeredgewidth=2)
                    else:
                        X.append(d[0])
                        Y.append(d[2])
                plt.plot(X, Y, c[0]+'x', markersize=5, markeredgewidth=2)

        if SHOW:
            plt.legend(['Full Adjacency', 'Limited Adjacency'],
                       numpoints=1, loc='upper left')
            plt.xlabel('Number of Cells')
            plt.ylabel('Average Qubit Usage')
            #plt.title('Average Qubit Usage vs. Number of Cells')
            plt.show()

        writeToFile(OUT, typ=typ)

    except KeyboardInterrupt:
        print 'Keyboard interrupt...'
        writeToFile(OUT, typ=typ)

if __name__ == '__main__':

    try:
        typ = sys.argv[1].lower()
        assert typ in ['dense', 'heur'], "Invalid embedding algorithm"
    except:
        print('Using default algorithm: dense')
        typ = 'dense'

    main(typ)
