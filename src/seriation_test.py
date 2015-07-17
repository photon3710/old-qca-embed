#---------------------------------------------------------
# Name: seriation.py
# Purpose: Test code for matrix seriation
# Author:    Jacob Retallick
# Created: 13.07.2015
# Last Modified: 13.07.2015
#---------------------------------------------------------

from parse_qca import parse_qca_file
from auxil import convert_to_lim_adjacency, convert_to_full_adjacency,\
    hash_problem
import sys
import numpy as np


N_TRIALS = 10

fname = None
#fname = '../dat/bench_circ/serial_adder_bench'


def main(fname):
    ''' '''

    # load QCA file
    cells, spacing, zones, J, feedback = parse_qca_file(fname, one_zone=True)
    J = convert_to_lim_adjacency(cells, spacing, J)

    N = len(cells)

    h = np.array([0]*N, dtype=float)
    g = np.array([0]*N, dtype=float)

    ninds = range(N)

    for i in xrange(N_TRIALS):

        # shuffle J
        np.random.shuffle(ninds)

        # generate scale
        sc = 200*np.random.random()

        hh = sc*h[ninds]
        JJ = sc*J[ninds, :][:, ninds]
        gg = sc*g[ninds]

        val, scale, inds = hash_problem(hh, JJ, gg)

        print val, scale

if __name__ == '__main__':
    if fname is None:
        try:
            fname = sys.argv[1]
        except:
            print('No filename given...')
            sys.exit()
    main(fname)
