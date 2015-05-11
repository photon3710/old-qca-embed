#!/usr/bin/python

#---------------------------------------------------------
# Name: wire_test.py
# Purpose: Energy gap computation for different wire lengths
# Author:    Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

from __future__ import division
import pylab as plt
import numpy as np
from time import time

from solve import solve
from build import createWire
from auxil import generateAdjDict, adjToCoef, EK0, GAMMA


N = 9
M = N
C = -1.

FS = 16

IMG_ROOT = '../img/'
SAVE = False


def main():
    '''Main loop'''

    Egap = np.zeros([N, N], dtype=float)

    Pol = [-1, 1] if C > 0 else [-1, -1]

    for n in xrange(1, N+1):
        for m in xrange(1, M+1):
            t = time()
            print 'Running (%d, %d) of (%d,%d)...' % (n, m, N, M),
            # generate wire and map to coefs
            cells, spacing = createWire([n, m], [-1, C, -1], Pol)
            adj, driver = generateAdjDict(cells, spacing)
            h, J = adjToCoef(adj)
            # get lowest spectrum
            gs, es, spec = solve(h, J, gamma=GAMMA*EK0, minimal=True)
            egap = spec[1]-spec[0]
            Egap[n-1, m-1] = egap
            print ' %.4f s' % (time()-t)

    plt.imshow(-100*np.flipud((Egap-Egap[0, 0])/Egap), cmap=plt.cm.jet,
               interpolation='none', extent=[1, N, 1, M])
    # colorbar
    cb = plt.colorbar()
    for l in cb.ax.get_yticklabels():
        l.set_fontsize(FS)
    plt.xlabel('N', fontsize=FS)
    plt.ylabel('M', fontsize=FS)
    plt.title('')
    xtp = [(i+.5)*(N-1)/N + 1/N for i in xrange(1, N+1)]
    plt.xticks(xtp, range(1, N+1))
    ytp = [(i+.5)*(M-1)/M + 1/M for i in xrange(1, M+1)]
    plt.yticks(ytp, range(1, M+1))
    plt.tick_params(axis='both', labelsize=FS)
    if SAVE:
        plt.savefig(IMG_ROOT + 'wire_eg.eps', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
