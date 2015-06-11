#!/usr/bin/python

#---------------------------------------------------------
# Name: zone.py
# Purpose: Container file for Zone class and associated functions
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 10.06.2015
#---------------------------------------------------------

import numpy as np


class Zone:
    '''A Zone object contains all the relevant information needed about a
    single clock zone of a circuit. This includes the current h and J
    coefficients, references and coupling matrices to all Zones which input
    or ouput to/from the zone, and the set of all spectra for each possible
    input configuration'''

    def __init__(self, zone, Gz, J, cells):
        '''Initialise a Zone object.

        inputs: zone    : key for the zone (for Gz)
                Gz      : directed zone graph from construct_zone_graph
                J       : kink-energies for all included cell interactions in
                         the circuit (including input and fixed cells)
                cells   : list of cell dicts for all cells in the circuit
        '''

        self.key = zone     # for record keeping

        self.inds = Gz.node[zone]['inds']    # indices of all simulated cells
        self.fixed = Gz.node[zone]['fixed']  # indices of all fixed cells
        self.drivers = Gz.node[zone]['drivers']  # indices of all input cells
        self.outputs = Gz.node[zone]['outputs']  # indices of all output cells

        self.N = len(self.inds)     # number of active cells in the zone

        # set up h and J containers
        self.h_fixed = np.zeros([1, self.N], dtype=float)
        self.h_driver = np.zeros([1, self.N], dtype=float)
        self.h = np.zeros([1, self.N], dtype=float)
        self.J = J[self.inds, :][:, self.inds]

        # get coupling to driver cells
        self.C_driver = J[self.drivers, :][:, self.inds]

        # get coupling to fixed cells
        self.C_fixed = J[self.fixed, :][:, self.inds]

        # compute h contribution from fixed cells
        pol_fixed = np.matrix([cells[ind]['pol'] for ind in self.fixed])
        self.h_fixed = .5*pol_fixed*self.C_fixed

        # set up interactions with input zones
        self.C_ins = {}     # dict of input zone coupling matrices
        self.outs = []      # list of ouput zones (only for searching)

        for in_zone in Gz.predecessors(zone):
            in_inds = Gz.node[in_zone]['inds']
            C = J[in_inds, :][:, self.inds]
            self.C_ins[in_zone] = C

        # store output zone keys
        for out_zone in Gz.successors(zone):
            self.outs.append(out_zone)

    def __str__(self):
        '''Default string covnersion for debugging'''

        outstring = '\n'*2+'*'*30 + '\n'
        
        # key
        outstring += 'key: {0}\n'.format(self.key)
        
        # indices
        outstring += 'inds: {0}\n'.format(self.inds)
        outstring += 'fixed: {0}\n'.format(self.fixed)
        outstring += 'drivers: {0}\n'.format(self.drivers)
        outstring += 'outputs: {0}\n'.format(self.outputs)
        
        # internal coupling matrices: non-zero flags only
        outstring += '\nC_driver:\n{0}\n'.format(1*(self.C_driver != 0))
        outstring += '\nC_fixed:\n{0}\n'.format(1*(self.C_fixed != 0))
        outstring += '\nJ:\n{0}\n'.format(1*(self.J != 0))
        
        # h coefficients
        outstring += '\nh_fixed:\n{0}\n'.format(np.round(self.h_fixed, 3))
        outstring += '\nh_driver:\n{0}\n'.format(np.round(self.h_driver, 3))
        outstring += '\nh:\n{0}\n'.format(np.round(self.h, 3))
        
        # input zones
        for z in self.C_ins:
            outstring += '\n{0}:\n{1}\n'.format(z, np.round(self.C_ins[z], 3))
            
        # output zones
        outstring += '\noutput-zones: {0}'.format(self.outs)
        
        return outstring