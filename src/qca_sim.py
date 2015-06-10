#!/usr/bin/python

#---------------------------------------------------------
# Name: qca_sim.py
# Purpose: General solver handling for QCA circuits
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 10.06.2015
#---------------------------------------------------------

from new_parse_qca import parse_qca_file
from new_auxil import convert_to_full_adjacency, convert_to_lim_adjacency

import sys


SOLVERS = {'rp': None,
           'dwave': None,
           'default': None}

assert 'default' in SOLVERS, 'Default solver has not be set...'


class Solution:

    def __init__(self):
        '''Initialise a Solution object'''
        pass

    def run_input(self, pols):
        ''''''
        pass

    def run_input_seq(self, pols):
        ''' '''
        pass
    
    def get_input_list(self):
        ''' '''
        pass
    
    def get_ouput_list(self):
        ''' '''
        pass
    
class Zone:
    
    def __init__(self, J):
        '''Initialise a Zone object'''        
        
        in_inds = []
        out_inds = []
        drivers = []
        fixed = []
        spectra = {}
    
    def solve(self):
        ''' '''
        pass
    


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
    
    # solve every zone for every possible set of inputs
    
    # construct solution object
    

if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print('No filename entered...')
        sys.exit()
