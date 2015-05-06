#!/usr/bin/python


import pylab as plt
import sys

FNAME_QBITS = '../summ/bench-qbits'
FNAME_PROB = '../summ/bench-prob'
FNAME_MODELS = '../summ/bench-models'
IMG_DIR = '../../img/'

QBITS = [20, 31, 73, 98, 120, 126, 192, 273]
FS = 14

SAVE = False


def loadData(fname):
    ''' '''

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open file: %s' % fname
        sys.exit()

    N = len(QBITS)
    data = {}
    data['dense'] = {'lim': [0]*N, 'full': [0]*N}
    data['heur'] = {'lim': [0]*N, 'full': [0]*N}

    i = 0
    for line in fp:
        if '#' in line or len(line) < 3:
            continue
        d = filter(None, line.strip().split())
        data['dense']['full'][i] = float(d[0])
        data['dense']['lim'][i] = float(d[2])
        data['heur']['full'][i] = float(d[1])
        data['heur']['lim'][i] = float(d[3])
        i += 1

    return data


def main():
    ''' '''

    ## QBITS

    data = loadData(FNAME_QBITS)

    N = len(data['dense']['lim'])

    for i in xrange(N):

        if data['dense']['lim'][i] > 0:
            plt.plot(QBITS[i], data['dense']['lim'][i], 'rx',
                     markersize=10, markeredgewidth=2)
        if data['dense']['full'][i] > 0:
            plt.plot(QBITS[i], data['dense']['full'][i], 'g^',
                     markersize=10, markeredgewidth=2)
        if data['heur']['lim'][i] > 0:
            plt.plot(QBITS[i], data['heur']['lim'][i], 'bv',
                     markersize=10, markeredgewidth=2)
        if data['heur']['full'][i] > 0:
            plt.plot(QBITS[i], data['heur']['full'][i], 'mo',
                     markersize=10, markeredgewidth=2)
    
    X = [0,200]
    plt.plot(X, X, 'k', ls='--', linewidth=2, alpha=.5)
    plt.xlim(X[0], X[-1])
    
    plt.legend(['DP: Limited', 'DP: Full', 'Heur: Limited', 'Heur: Full'],
               loc='lower right', fontsize=FS, numpoints=1)
    plt.tick_params(axis='both', which='major', labelsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Average Qubit Usage', fontsize=FS)
    if SAVE:
        plt.savefig(IMG_DIR+'summ-qbits.eps', bbox_inches='tight')
    plt.show()

    ## PROBS

    data = loadData(FNAME_PROB)

    N = len(data['dense']['lim'])

    for i in xrange(N):

        if data['dense']['lim'][i] > 0:
            plt.plot(QBITS[i], 100./data['dense']['lim'][i], 'rx',
                     markersize=10, markeredgewidth=2)
        if data['dense']['full'][i] > 0:
            plt.plot(QBITS[i], 100./data['dense']['full'][i], 'g^',
                     markersize=10, markeredgewidth=2)
        if data['heur']['lim'][i] > 0:
            plt.plot(QBITS[i], 100./data['heur']['lim'][i], 'bv',
                     markersize=10, markeredgewidth=2)
        if data['heur']['full'][i] > 0:
            plt.plot(QBITS[i], 100./data['heur']['full'][i], 'mo',
                     markersize=10, markeredgewidth=2)

    plt.legend(['DP: Limited', 'DP: Full', 'Heur: Limited', 'Heur: Full'],
               loc='upper right', fontsize=FS, numpoints=1)
    plt.tick_params(axis='both', which='major', labelsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Single Trial Success Probability', fontsize=FS)
    plt.ylim([-.02, 1.05])
    if SAVE:
        plt.savefig(IMG_DIR+'summ-probs.eps', bbox_inches='tight')
    plt.show()

    ## MODELS

    data = loadData(FNAME_MODELS)

    print data

    N = len(data['dense']['lim'])

    for i in xrange(N):

        if data['dense']['lim'][i] > 0:
            plt.plot(QBITS[i], data['dense']['lim'][i], 'rx',
                     markersize=10, markeredgewidth=2)
        if data['dense']['full'][i] > 0:
            plt.plot(QBITS[i], data['dense']['full'][i], 'g^',
                     markersize=10, markeredgewidth=2)
        if data['heur']['lim'][i] > 0:
            plt.plot(QBITS[i], data['heur']['lim'][i], 'bv',
                     markersize=10, markeredgewidth=2)
        if data['heur']['full'][i] > 0:
            plt.plot(QBITS[i], data['heur']['full'][i], 'mo',
                     markersize=10, markeredgewidth=2)

    plt.legend(['DP: Limited', 'DP: Full', 'Heur: Limited', 'Heur: Full'],
               loc='upper right', fontsize=FS, numpoints=1)
    plt.tick_params(axis='both', which='major', labelsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Average Maximum Model Size', fontsize=FS)
    #plt.ylim([-.02, 1.05])
    if SAVE:
        plt.savefig(IMG_DIR+'summ-models.eps', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
