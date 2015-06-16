#!/usr/bin/python

#---------------------------------------------------------
# Name: zone.py
# Purpose: Container file for Zone class and associated functions
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 11.06.2015
#---------------------------------------------------------

import numpy as np
import itertools
import xml.etree.ElementTree as ET  # xml input and output
from auxil import gen_pols
from networkx import DiGraph

from pprint import pprint

class Zone:
    '''A Zone object contains all the relevant information needed about a
    single clock zone of a circuit. This includes the current h and J
    coefficients, references and coupling matrices to all Zones which input
    or ouput to/from the zone, and the set of all spectra for each possible
    input configuration'''

    def __init__(self, *args):
        '''Initialise a Zone object.

        inputs: zone    : key for the zone (for Gz)
                Gz      : directed zone graph from construct_zone_graph
                J       : kink-energies for all included cell interactions in
                         the circuit (including input and fixed cells)
                cells   : list of cell dicts for all cells in the circuit
        '''
        if args[0] == 'from_file':
            self.read_from_file(args[1])

        elif type(args[1]) is DiGraph:
            zone = args[0]
            Gz = args[1]
            J = args[2]
            cells = args[3]

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
        else:
            print 'Invalid arguments for Zone constructor'
            raise TypeError


    def __str__(self):
        '''Default string conversion for debugging'''

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

    def get_drivers(self):
        '''Return the indices of the driver cells'''
        return self.drivers

    def get_outputs(self):
        '''Return the indices of the output cells'''
        return self.outputs

    def get_indices(self):
        '''Return the indices of simulated cells'''
        return self.inds

    def solve_all(self, solver=None, **kwargs):
        '''Solve all possible input configurations for the zone using the
        specified solver

        inputs: solver  : A solver function which takes h and J values and
                         returns some form of an output to be used externally
                kwargs  : arguments to pass to solver

        outputs: cases  : A list of all considered configurations of driver and
                         zone input polarizations
                 outs   : A dict of the corresponding output dicts from solver
                 z-order: order of zones as they appear in outs keys
        '''

        # relative indices in each input-zone of cells which couple in.
        c_inds = {}
        for z in self.C_ins:
            z_inds = list(set(np.nonzero(self.C_ins[z])[0].tolist()))
            c_inds[z] = z_inds

        # input zone order
        z_order = sorted(c_inds)

        # set of all driver inputs
        driver_inds = gen_pols(len(self.drivers))
        if len(driver_inds) == 0:
            driver_inds = [()]

        # set of all inputs for each input-zone
        inzone_inds = [gen_pols(len(c_inds[z])) for z in c_inds]

        # combine into list product for easy looping
        cases = list(itertools.product(driver_inds, *inzone_inds))

        # run each configuration, solve
        outs = {}

        for case in cases:
            print case
            # compute driver contribution to h
            self.h_driver = .5*np.matrix(case[0])*self.C_driver
            # for each input zone, map pol key to state and get h contr.
            self.h = self.h_fixed + self.h_driver
            for i in xrange(1, len(case)):
                key = z_order[i-1]
                pol = np.zeros([1, len(self.C_ins[key])], dtype=float)
                pol[0, c_inds[key]] = case[i]
                self.h += .5*np.asmatrix(pol)*self.C_ins[key]

            outs[case] = solver(self.h, -.5*self.J, **kwargs)

        return cases, outs, z_order

    def write_to_file(self, parent):
        '''Add zone information to given xml parent node'''
        this_zone = ET.SubElement(parent, 'zone', attrib = \
            {'clk' : str(self.key[0]), 'i' : str(self.key[1])})

        ET.SubElement(this_zone, 'J', attrib = \
            {'J' : str(np.array(self.J).tolist())})

        ET.SubElement(this_zone, 'h', attrib = \
            {'h' : str(np.array(self.h).tolist())})

        cell_indices = {
            'inds' : str(self.inds),
            'fixed' : str(self.fixed),
            'drivers' : str(self.drivers),
            'outputs' : str(self.outputs)
        }
        ET.SubElement(this_zone, 'cell_indices', attrib = cell_indices)

    def read_from_file(self, node):
        '''Construct a Zone object from its xml node'''
        self.zone = int(node.get('clk')), int(node.get('i'))

        cell_indices = node.find('cell_indices')
        self.inds = string_to_list(cell_indices.get('inds'))
        self.fixed = string_to_list(cell_indices.get('fixed'))
        self.drivers = string_to_list(cell_indices.get('drivers'))
        self.outputs = string_to_list(cell_indices.get('outputs'))

        self.J = string_to_2d_np_arr(node.find('J').get('J'))
        self.h = strings_to_2d_np_arr(node.find('h').get('h'))

        self.N = len(self.inds)


### HELPER FUNCTIONS

def string_to_list(string):
    '''inverse of str(list[])'''
    l = []
    for i in string.strip('[]').split(','):
        if i != '':
            l.append(int(i.strip(' ')))
    return l

def string_to_2d_np_arr(string):
    '''inverse function of str(np.array)'''
    arr = []
    for i in string.split('['):
        l = []
        if i.strip(' ') != '':
            for j in i.strip(',[ ]').split(','):
                l.append(float(j.strip('[ ]')))
            arr.append(l)
    return np.asarray(arr)

### WRITING/READING TO XML FILES

def write_zones_to_xml(Zones, file_name):
    '''writes all nodes in Zones to an xml specified by file_name'''
    root = ET.Element('all_zones')
    for key in Zones:
        zone = Zones[key]
        zone.write_to_file(root)
    tree = ET.ElementTree(root)
    tree.write(file_name)

def read_zones_from_xml(file_name):
    '''returns a list of all zone elements stored in an xml'''
    tree = ET.parse(file_name)
    root = tree.getroot()
    zones = []
    for zone in root.findall('zone'):
        zones.append(Zone('from_file', zone))

    return zones
