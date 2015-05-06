#!/usr/bin/python

from __future__ import division
from dwave_sapi import find_embedding, get_chimera_adjacency

import sys
import os
import datetime
import time

from parse_qca import parseQCAFile
from auxil import generateAdjDict, convertToNearestNeighbour


MAX_COUNT = 100
FLAGSOL = True

RUN_FULL = False

root_dir = '../sols/temp/'

root_dir += 'full/' if RUN_FULL else 'lim/'


info_file = root_dir+'info'


def formatHeuristic(adjacency):
    ''' formats the values from the adjacency dict into a structure needed to
    run find_embeddings in the SAPI module...

    find_embedding(S, S_size, A, A_size, verbose=1)'''

    S = {}
    N = len(adjacency.keys())

    index_map = adjacency.keys()

    keys = []

    for i in xrange(N):
        adj = adjacency[index_map[i]]
        for j in xrange(1, len(adj)):
            cell_ind, Ek = adj[j]
            keys.append((i, index_map.index(cell_ind)))

    for i in xrange(N):
        for j in xrange(N):
            if (i, j) in keys or (j, i) in keys:
                S[(i, j)] = 1
            else:
                S[(i, j)] = 0

    return S


def runHeuristic(adjacency, max_count, flagSol=False):
    '''run the heuristic method given the adjacency dict'''

    M = 8
    N = M
    L = 4

    S = formatHeuristic(adjacency)
    S_size = len(adjacency.keys())

    A = get_chimera_adjacency(M, N, L)
    A_size = M * N * L * 2

    trial_num = 0
    success_num = 0
    count = 0

    key_map = adjacency.keys()
    good_embeds = []
    times = []

    while count < max_count:

        t1 = time.clock()
        embeddings = find_embedding(S, S_size, A, A_size, verbose=0, tries=1)
        t2 = time.clock()
        trial_num += 1

        if len(embeddings) == S_size:    # successful embedding
            success_num += 1
            # restucture embeddings
            mapped_embed = {}
            for i in xrange(len(key_map)):
                mapped_embed[key_map[i]] = embeddings[i]
            good_embeds.append(mapped_embed)
            print 'solution ' + str(success_num) + ' found...'

            times.append(t2-t1)

        count = success_num if flagSol else trial_num

    print '\n'*2+'Embedded ' + str(success_num) + ' of ' + str(trial_num) + ' attempts' + '\n'*2

    return good_embeds, times, [success_num, trial_num]


def writeToFile(filename, embeddings, times, counts):
    ''' '''

    print '*'*40
    print 'Writing to file...\n'

    dir_name = filename.rpartition('/')[2]

    solution_dir = root_dir+dir_name+'/session'
    count = 0

    if not dir_name in os.listdir(root_dir):
        os.mkdir(root_dir+dir_name)

    print 'Building solution directory...'
    while count < 5:
        try:
            print 'Trying: ' + solution_dir + str(count) + '...'
            os.makedirs(solution_dir+str(count))
        except:
            count += 1
            continue
        solution_dir += str(count)+'/'
        break

    print '\n' + 'Solution Directory: ' + solution_dir

    fname_log = solution_dir + 'log'

    fp_log = open(fname_log,'w')

    # write log header

    fp_log.write('Source file: ' + filename + '\n')

    time_str = datetime.datetime.today().isoformat()
    time_str = time_str.replace('T', ' ').partition('.')[0]

    fp_log.write('Timestamp: ' + time_str + '\n')
    fp_log.write('Trials: ' + str(max(counts)) + '\n')
    fp_log.write('Successes: ' + str(min(counts)) + '\n')

    for i in xrange(len(embeddings)):
        embedding = embeddings[i]
        fname_sol = solution_dir + 'solution' + str(i)
        fp_log.write('\nSolution ' + str(i) +'\n')
        fp_log.write('file: ' + fname_sol + '\n')
        fp_log.write('time (seconds): ' + str(times[i]) + '\n')
        fp_sol = open(fname_sol, 'w')
        for key in embedding:
            fp_sol.write(str(key) + ' : ' + str(embedding[key])+'\n')
        fp_sol.close()

    fp_log.close()

    print '\nWriting complete...'



def run_embed(filename, max_count, flagSol=False):

    print '*'*40
    # parse QCA file
    cells, spacing = parseQCAFile(filename)

    # generate adjacency list
    adjacency, drivers = generateAdjDict(cells, spacing)
    
    if not RUN_FULL:
        adjacency = convertToNearestNeighbour(adjacency, drivers)

    if False:
        print '*'*40+'\n'*3
        print 'Drivers: ' + str(drivers)
        print '*'*40+'\n'*3
        print 'Adjacency: '
        for key in adjacency.keys():
            print str(key) + ' : ' + str(adjacency[key])

    print '*'*40
    print '\n'

    # format entries for dense placement algorithm

    # format entries for heuristic algorithm
    embeddings, times, counts = runHeuristic(adjacency, max_count, flagSol)

    writeToFile(filename, embeddings, times, counts)




if __name__ == '__main__':

    if len(sys.argv) == 1:
        print 'Error: enter QCADesigner file name'
        sys.exit()

    try:
        fname = str(sys.argv[1])
    except:
        print 'Invalid input... should not have happened'
        sys.exit()

    run_embed(fname, MAX_COUNT, FLAGSOL)
