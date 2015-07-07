#!/usr/bin/python

from __future__ import division
import pylab as plt
import os
import sys
from math import ceil

FILE_ROOT = '../512gen/0/'
IMG_ROOT = '../../img/'

SHOW_MIN = False

FS = 18

SAVE = False

BW = 20
MAX_RANGE = 350

REL = False


def load_outfile(fname):
    ''' '''

    print 'reading file: %s' % fname

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open file: %s' % fname
        return None

    data = {'lim': [], 'full': []}

    for line in fp:
        if len(line) < 3:
            continue
        elif '<lim\>' in line:
            dat = data['lim']
        elif '<full\>' in line:
            dat = data['full']
        else:
            d = line.strip().split()
            n = int(d[0])
            av = float(d[1])
            mi = int(d[2])
            try:
                t = float(d[3])
                p = float(d[4])
                dat.append([n, av, mi, t, p])
            except KeyError:
                dat.append([n, av, mi])

    #print '\t(%d,%d)' % (len(data['lim']), len(data['full']))

    return data


def main(file_dir, img_dir, bw):
    ''' '''

    fnames = os.listdir(file_dir)
    fnames = map(lambda s: file_dir+s, fnames)

    all_data = {'lim': [], 'full': []}

    for fname in fnames:
        data = load_outfile(fname)
        if not data is None:
            all_data['lim'] += data['lim']
            all_data['full'] += data['full']

    # plot qubit usage

    plt.figure('qubit-usage')
    plt.clf()
    ind = 2 if SHOW_MIN else 1
    for key in ['full', 'lim']:
        c = ['g', 'r'] if key == 'full' else ['w', 'm']
        s = 'x' if key == 'full' else 'o'
        X, Y = [], []
        for d in all_data[key]:
            #print '%s: \t %s' % (key, str(d))
            if d[ind] == -1:
                continue
                plt.plot(d[0], 0, c[1]+s, markersize=5, markeredgewidth=2)
            else:
                X.append(d[0])
                Y.append(d[ind])
        if REL:
            Y = [Y[i]*1./X[i] for i in xrange(len(Y))]
        if key == 'full':
            plt.plot(X, Y, c[0]+'x', markersize=5, markeredgewidth=2)
        else:
            plt.plot(X, Y, c[0]+'o', markersize=5, markeredgewidth=2,
                     markeredgecolor='blue')
        X = range(1, max(X)+1)
    plt.plot(X, X, 'k', ls='--', linewidth=2, alpha=.5)
    plt.xlim(X[0], X[-1])

    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc='best', fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    ylab = ('Minimum' if SHOW_MIN else 'Average') + 'Qubit Usage'
    plt.ylabel(ylab, fontsize=FS)
    plt.tick_params(axis='both', labelsize=FS)
    #plt.title('Average Qubit Usage vs. Number of Cells')
    if SAVE:
        plt.savefig(img_dir+'-gen-qbits.eps', bbox_inches='tight')
    plt.show(block=False)

    for key in all_data:
        print('{0}: {1} trials'.format(key, len(all_data[key])))

    # plot success probability

    plt.figure('success-prob')
    plt.clf()

    for key in ['full', 'lim']:

        c = ['g', 'r'] if key == 'full' else ['w', 'm']
        max_bin = max(map(lambda x: x[0], all_data[key]))
        num_bins = int(ceil(max_bin/bw))+1

        success = [0 for _ in xrange(num_bins)]
        fail = [0 for _ in xrange(num_bins)]

        for d in all_data[key]:
            if d[1] > 0:
                success[int(d[0]/bw)] += 1
            else:
                fail[int(d[0]/bw)] += 1

        rate = []
        for i in xrange(num_bins):
            if success[i]+fail[i] > 0:
                rate.append(success[i]/(success[i]+fail[i]))
            else:
                rate.append(0)
        X = map(lambda x: bw*(.5+x), range(num_bins))
        if key is 'full':
            plt.plot(X, rate, c[0]+'x', markersize=8, markeredgewidth=3)
        else:
            plt.plot(X, rate, c[0]+'o', markersize=8, markeredgewidth=3,
                     markeredgecolor='blue')
    loc = 'best'
    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc=loc, fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Success probability', fontsize=FS)
    plt.xlim([0, MAX_RANGE])
    plt.tick_params(axis='both', labelsize=FS)
    if SAVE:
        plt.savefig(img_dir+'-gen-prob.eps', bbox_inches='tight')
    plt.show(block=False)

    # plot run-times

    plt.figure('run-times')
    plt.clf()

    for key in ['full', 'lim']:

        c = 'gx' if key == 'full' else 'wo'
        X = []
        Y = []
        for d in all_data[key]:
            try:
                n, av, mi, t, p = d
            except:
                continue
            if t < 0:
                continue
            X.append(n)
            Y.append(t)
        if key is 'full':
            plt.plot(X, Y, c, markersize=5, markeredgewidth=2)
        else:
            plt.plot(X, Y, c, markersize=5, markeredgewidth=2,
                     markeredgecolor='blue')

    loc = 'best'
    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc=loc, fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Run-time (s)', fontsize=FS)
    plt.tick_params(axis='both', labelsize=FS)
    if SAVE:
        plt.savefig(img_dir+'-gen-runtimes.eps', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':

    try:
        typ = sys.argv[1].strip().lower()
    except:
        typ = 'dense'

    try:
        bw = int(sys.argv[2])
    except:
        bw = BW

    if typ == 'dense':
        file_dir = FILE_ROOT + 'dense/'
        img_dir = IMG_ROOT + 'dense'
    else:
        file_dir = FILE_ROOT + 'heur/'
        img_dir = IMG_ROOT + 'heur'

    main(file_dir, img_dir, bw)
