from __future__ import division

import re
import pylab as plt
import numpy as np


FS = 14
LOG = False

DENSE_ROOT = '../bench/dense/full/ser-add/sol'
HEUR_ROOT = '../bench/heur/full/ser-add/solution'


def pinch(st, l, r):
    return st.partition(l)[2].rpartition(r)[0]


def strToQbit(s):
    ''' convert a formatted string to a qubit tuple '''
    r, c, h, l = map(int, s.strip()[1:-1].split(','))
    return (r, c, h, l)


def readSub(data):
    ''' '''

    read_flag = False
    rx = re.compile('<.*>')
    dat = []
    while data:
        line = data.pop(0)
        if len(line) < 3 or line.strip()[0] == '#':
            continue
        if not read_flag:
            m = rx.search(line)
            if not m is None:
                name = m.group(0)[1:-1]
                read_flag = True
                continue
        if read_flag:
            if '</%s>' % name in line:
                break
            dat.append(line.strip('\n'))

    return dat


def loadDense(fname):
    '''Load a Dense Placement solutions file'''

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open filename: %s' % fname
        return None

    data = fp.readlines()
    fp.close()

    header_dat = readSub(data)
    dis_qbit_dat = readSub(data)
    qbit_dat = readSub(data)
    path_dat = readSub(data)

    # parse header data
    M = int(header_dat[0].split()[-1])
    N = int(header_dat[1].split()[-1])
    L = int(header_dat[2].split()[-1])

    # parse dis_qbit data
    dis_qbits = dis_qbit_dat
    dis_qbits = []

    # parse qbit data
    qbits = {}
    for line in qbit_dat:
        cell, qbit = line.split(':')
        cell = int(cell)
        qbit = strToQbit(qbit)
        qbits[cell] = qbit

    # parse paths data
    paths = {}
    for line in path_dat:
        key, path = line.split(':')
        key = tuple(map(int, key.strip().split(';')))
        path = map(strToQbit, path.split(';'))
        paths[key] = path

    sol_dict = {}
    sol_dict['M'] = M
    sol_dict['N'] = N
    sol_dict['L'] = L
    sol_dict['dis_qbits'] = dis_qbits
    sol_dict['qbits'] = qbits
    sol_dict['paths'] = paths

    return sol_dict


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


def main():

    dense_hist = {}
    heur_hist = {}

    for i in [10]:

        dense = loadDense(DENSE_ROOT + str(i))
        heur = loadHeur(HEUR_ROOT + str(i))

        dense_lens = map(lambda x: len(x)-1, dense['paths'].values())
        heur_lens = map(len, heur.values())

        max_len = max(dense_lens + heur_lens)

        for l in xrange(1, max_len+1):
            if not l in dense_hist:
                dense_hist[l] = 0
            if not l in heur_hist:
                heur_hist[l] = 0
            dense_hist[l] += dense_lens.count(l)
            heur_hist[l] += heur_lens.count(l)

    max_len = max(dense_hist.keys() + heur_hist.keys())
    dense_hist_l = []
    heur_hist_l = []
    for k in xrange(max_len+1):
        if k in dense_hist:
            dense_hist_l.append(dense_hist[k])
        else:
            dense_hist_l.append(0)
        if k in dense_hist:
            heur_hist_l.append(heur_hist[k])
        else:
            heur_hist_l.append(0)

    dense_mean = map(lambda x: x/sum(dense_hist_l), dense_hist_l)
    heur_mean = map(lambda x: x/sum(heur_hist_l), heur_hist_l)

    axes = plt.figure().add_subplot(111)

    # Chain Lengths
    bar_width = .2
    dx = .1
    if LOG:
        dense_mean = np.log10(1+100*np.array(dense_mean, dtype=float))
        heur_mean = np.log10(1+100*np.array(heur_mean, dtype=float))

    axes.bar([i-dx for i in xrange(max_len+1)], dense_mean,
             width=bar_width, color='blue')
    axes.bar([i+bar_width for i in xrange(max_len+1)], heur_mean,
             width=bar_width, color='white')

    plt.xlabel('Chain/Model Size', fontsize=FS)
    if LOG:
        plt.ylabel('Occurrence probability: log$(1+100P)$', fontsize=FS)
    else:
        plt.ylabel('Occurrence probability', fontsize=FS)
    #plt.title('Distribution of Qubit Group Sizes',fontsize=FS)
    axes.tick_params(axis='both', which='major', labelsize=FS)
    plt.legend(['Dense Placement', 'Heuristic'], fontsize=FS)
    plt.xlim([0, 15])
    plt.show()

    print dense_mean
    print heur_mean


if __name__ == '__main__':
    main()
