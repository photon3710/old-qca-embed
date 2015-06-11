#!/usr/bin/python

#---------------------------------------------------------
# Name: qca_sim.py
# Purpose: General solver handling for QCA circuits
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 10.06.2015
#---------------------------------------------------------

from new_auxil import convert_to_full_adjacency, convert_to_lim_adjacency, \
    construct_zone_graph

from new_parse_qca import parse_qca_file

#from solution import Solution
from zone import Zone

import sys
from pprint import pprint


SOLVERS = {'rp': None,
           'dwave': None,
           'default': None}

assert 'default' in SOLVERS, 'Default solver has not be set...'


def qca_sim(fn, **kwargs):
    '''Simulate all possible outcomes for each zone of a given qca circuit'''

    # parse QCADesigner file
    cells, spacing, zones, J = parse_qca_file(fn)

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
    Gz = construct_zone_graph(cells, zones, J, show=True)
    Zones = {key: Zone(key, Gz, J, cells) for key in Gz.nodes()}
    
    for z in sorted(Zones):
        print(str(Zones[z]))

    # solve every zone for every possible set of inputs

    # construct solution object


if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print('No filename entered...')
        sys.exit()
        
    qca_sim(fn)
