#!/usr/bin/python

#---------------------------------------------------------
# Name: qca_sim.py
# Purpose: General solver handling for QCA circuits
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 10.06.2015
#---------------------------------------------------------

from auxil import convert_to_full_adjacency, convert_to_lim_adjacency, \
    construct_zone_graph, gen_pols

from parse_qca import parse_qca_file
from rp_solve import rp_solve

from solution import Solution
from zone import Zone, write_zones_to_xml

import sys
import numpy as np
from pprint import pprint


def rp_solver(h, J, gam=0.):
    '''Solve the low-energy spectra of a zone using rp_solve'''

    h = h.tolist()[0]
    gam = np.max(np.abs(J))*gam

    E, states, Eps, pstates, state_pols, quick_fix = rp_solve(h, J, gam=gam)

    out = {'Es': E,
           'states': states,
           'Eps': Eps,
           'pstates': pstates,
           'state_pols': state_pols}

    return out


SOLVERS = {'rp': rp_solver,
           'dwave': None,
           'default': rp_solver}

assert 'default' in SOLVERS, 'Default solver has not be set...'


def qca_sim(fn, **kwargs):
    '''Simulate all possible outcomes for each zone of a given qca circuit'''

    # parse QCADesigner file
    cells, spacing, zones, J, feedback = parse_qca_file(fn, show=False)

    # check for specification of adjacency type
    if 'adj' in kwargs:
        if kwargs['adj'].lower() == 'full':
            J = convert_to_full_adjacency(cells, spacing, J)
        elif kwargs['adj'].lower() == 'lim':
            J = convert_to_lim_adjacency(cells, spacing, J)

    # set the solver type
    if 'solver' in kwargs:
        try:
            solver = SOLVERS[kwargs['solver'].lower()]
        except:
            print('No solver specified. Using default...')
            solver = SOLVERS['default']
    else:
        solver = SOLVERS['default']

    # set up zone formulation
    Gz = construct_zone_graph(cells, zones, J, feedback, show=False)
    Zones = {key: Zone(key, Gz, J, cells, feedback) for key in Gz.nodes()}
    # solve every zone for every possible set of inputs
    solution = Solution(Gz)
    for i in xrange(len(zones)):
        for j in xrange(len(zones[i])):
            key = (i, j)
            cases, outs, z_order = Zones[key].solve_all(solver)
            solution.add_zone(Zones[key], outs, z_order)

    # write solution to file
    solution.write_to_file('./qca_sim_test')

    # run solution inputs (single sets)
    input_inds = solution.get_inputs()

    input_pols = gen_pols(len(input_inds))  # set of all possible input pols

    print 'start'
    out = {}    # dict of outputs lists for each input polarization list
    for pols in input_pols:
        input_pol = {(0,0): [(pols, )]}
        out[pols] = solution.run_input(input_pol)

    pprint(out)


if __name__ == '__main__':

    try:
        if False:
            fn = sys.argv[1]
        else:
            fn = 'test_circuits/split'
    except:
        print('No filename entered...')
        sys.exit()

    qca_sim(fn)
