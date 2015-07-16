#!/usr/bin/python

#---------------------------------------------------------
# Name: solution.py
# Purpose: Container file for Solution class and associated functions
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 12.06.2015
#---------------------------------------------------------

import xml.etree.ElementTree as ET
import string_conversion as sc
from zone import write_zones_to_xml, read_zones_from_xml
import networkx as nx
import numpy as np
from pprint import pprint

ZONE_EXT = '_zones'
XML_EXT = '.xml'

class Solution:
    '''A Solution object is constructed from the results of running a solver on
    a QCA circuit. It is assumed that the solver approach is to find the
    spectra of each clock zone for all possible inputs. The Solution object
    contains methods to analyse these results in order to determine propoerties
    of the whole circuit.'''

    def __init__(self, *args):
        '''Initialise a Solution object'''
        if args[0] == 'from_file':
            self.load_from_file(args[1])

        elif type(args[0]) is nx.DiGraph:
            self.Gz = args[0]        # directed zone graph
            self.zones = []     # list of zone object
            self.zone_dict = {} # dictionary of zone objects
            self.z_orders = []  # order of input zones for sols indexing
            self.sols = []      # list of sol dictionaries for each zone
        else:
            print 'Invalid arguments for Solution constructor'
            raise TypeError

    def add_zone(self, zone, fmt_out, z_order):
        '''Add a solved zone to the Solution object. Variable fmt_out
        must be a dict with keys given by 'cases' countaining dicts of the form
            fmt_out[case]['Es']     -> list of eigenenergies
            fmt_out[case]['states'] -> list of eigenstates
            fmt_out[case]['Eps']    -> list of product state energies
            fmt_out[case]['pstates']-> list of product states
        In general, only the product state parameters need to be populated to
        estimate 'computing in the ground state' for low tunnelling.

        inputs: zone    : Zone object
                fmt_out : formatted output from a solver (see note above)
                z_order : order of input-zones in the keys of fmt_out
        '''

        # check format of fmt_out
        try:
            for case in fmt_out:
                assert 'Es' in fmt_out[case], 'No Es key found'
                assert 'states' in fmt_out[case], 'No states key found'
                assert 'Eps' in fmt_out[case], 'No Eps key found'
                assert 'pstates' in fmt_out[case], 'No pstates key found'
        except Exception as e:
            print(e.message)
            print 'Invalid fmt_out format... see documentation'
            return False

        self.zones.append(zone)
        self.zone_dict[zone.key] = zone
        self.z_orders.append(z_order)
        self.sols.append(fmt_out)

    def write_to_file(self, fn):
        '''Create xml tree for the Solution object and write to file'''
        root = ET.Element('solution')

        Gz = ET.SubElement(root, 'Gz')
        edges = ET.SubElement(Gz, 'edges')
        for edge in self.Gz.edges():
            e = ET.SubElement(edges, 'edge')
            e.text = str(edge)
        nodes = ET.SubElement(Gz, 'nodes')
        for node in self.Gz.nodes():
            n = ET.SubElement(nodes, 'node')
            n.text = str(node)

        write_zones_to_xml(self.zones, fn  + ZONE_EXT)

        z_o = ET.SubElement(root, 'z_orders')
        z_o.text = str(self.z_orders)
        sols = ET.SubElement(root, 'sols')
        for zone in self.sols:
            zone_el = ET.SubElement(sols, 'zone')
            for sol in zone:
                case = ET.SubElement(zone_el, 'case')
                ET.SubElement(case, 'case_key').text = str(sol)
                for key in zone[sol]:
                    ET.SubElement(case,str(key)).text = str(zone[sol][key])

        tree = ET.ElementTree(root)
        tree.write(fn + XML_EXT)

    def load_from_file(self, fn):
        '''load solution from xml'''
        tree = ET.parse(fn + XML_EXT)
        root = tree.getroot()
        Gz = nx.DiGraph()
        for edge in root.iter('edge'):
            e = sc.str_to_2d_tuple(edge.text)
            Gz.add_edge(e[0], e[1])
        for node in root.iter('node'):
            Gz.add_node(sc.str_to_tuple(node.text))
        self.Gz = Gz
        self.zone_dict = {}
        self.zones = read_zones_from_xml(fn + ZONE_EXT)
        for zone in self.zones:
            self.zone_dict[zone.key] = zone

        z_o = root.find('z_orders')
        self.z_orders = sc.str_to_2d_list_of_tup(z_o.text)

        sols = []

        for zone in root.iter('zone'):
            zone_dict = {}
            for case in zone.iter('case'):
                case_dict = {}
                key = sc.str_to_2d_tuple(case.find('case_key').text)

                #states
                case_dict['states'] = sc.str_to_list_of_nparr(case.find('states').text, integer=False)
                #pstates
                case_dict['pstates'] = sc.str_to_list_of_tup(case.find('pstates').text)
                #state_pols
                case_dict['state_pols'] = sc.str_to_list_of_nparr(case.find('state_pols').text, integer=False)
                #Eps
                case_dict['Eps'] = sc.str_to_tuple(case.find('Eps').text, integer=False)
                #Es
                case_dict['Es'] = sc.str_to_nparr(case.find('Es').text, integer=False)

                zone_dict[key] = case_dict
            sols.append(zone_dict)
        self.sols = sols

    def get_inputs(self):
        '''Get list of all input cell indices'''
        input_indices = []
        for zone in self.zones:
            for i in zone.drivers:
                input_indices.append(i)
        return input_indices

    def run_input_single(self, zone_inputs, sequence=False, first_run=True):
        '''Use the solution information to determine the ouput polarizations
        for a single set of input polarizations'''

        # ASSUMES ZONES ARE IN ORDER
        inputs = {}
        final_pol = {}
        n_cells = 0

        # only used for sequences
        defaulted_inputs = None
        visited_zones = []

        for i in range (len(self.zones)):
            zone = self.zones[i]
            zone_sols = self.sols[i]
            visited_zones.append(zone.key)

            # check to see if all inputs have been set (i.e. feedback loops)
            # if it hasn't revert to default polarization
            if sequence and first_run:
                DEFAULT_POL = 1
                for input_zone in zone.C_ins:
                    if input_zone not in visited_zones:
                        new_input = []
                        for j in range(len(zone.C_ins[input_zone])):
                            if any(zone.C_ins[input_zone][j]):
                                new_input.append((DEFAULT_POL,))

                        new_pols = []
                        defaulted_inputs = list(inputs[zone.key])
                        while inputs[zone.key]:
                            old_pol = inputs[zone.key].pop()
                            new_out_pol = list(old_pol)
                            for out_pol in new_input:
                                new_pols.append(tuple(new_out_pol + new_input))
                        inputs[zone.key] = new_pols
                        zone_inputs[zone.key] = new_pols


            # gets all possible output polarizations (to next zone or final)
            # with regards to the possible inputs
            out_pols, final_pol[zone.key] = get_output_polarizations\
                (zone, self.zone_dict, zone_sols, zone_inputs[zone.key], n_cells)

            # removes defaulted inputs from zone_inputs
            if defaulted_inputs:
                inputs[zone.key] = list(defaulted_inputs)
                del zone_inputs[zone.key]
                defaulted_inputs = None

            # remove the duplicate polarizations
            for key in out_pols:
                out_pols[key] = list(set(out_pols[key]))

            # flag to identify the first predecessor to the next clock cycle
            first_to_next = False
            # go through all the zones and format alll the future inputs
            for key in zone.outs:
                if key not in inputs:
                    first_to_next = True
                    inputs[key] = [((),out_pols[key][0])]

                if first_to_next:
                    for j in range(1,len(out_pols[key])):
                        inputs[key].append(((),(out_pols[key][j])))
                else:
                    new_pols = []
                    while inputs[key]:
                        old_pol = inputs[key].pop()
                        new_out_pol = list(old_pol)
                        for out_pol in out_pols[key]:
                            new_pols.append(tuple(new_out_pol + [out_pol]))
                    inputs[key] = new_pols

            # append to zone_inputs
            for key in zone.outs:
                next_zone = self.zone_dict[key]
                complete_input = True
                for input_zone in next_zone.C_ins:
                    if input_zone not in visited_zones:
                        complete_input = False

                if complete_input:
                    if key in zone_inputs:
                        if set(zone_inputs[key]) != set(inputs[key]):
##                            print 'writing over zone_inputs[%s] with %s'\
##                                %(str(key), str(inputs[key]))
                            zone_inputs[key] = inputs[key]
                    else:
                        zone_inputs[key] = list(inputs[key])

            n_cells += zone.N+len(zone.fixed)+len(zone.drivers)

        if sequence:
            return final_pol, zone_inputs

        return final_pol

    def run_input_sequence(self, pol_seq):
        '''Use the solution information to deterime the output sequence for a
        sequence of input sets (relevant for circuits with feedback'''
        defaulted_inputs = {}
        final_pol = None
        zone_inputs = pol_seq
        first_run = True

        # run input single until the inputs to zones remain const across runs
        while True:
            final_pol,ret_inputs = self.run_input_single(dict(zone_inputs),\
                sequence=True, first_run=first_run)

            # check returned inputs to previously returned inputs
            if ret_inputs == zone_inputs:
                zone_inputs = ret_inputs
                break

            zone_inputs = ret_inputs
            first_run = False
        return final_pol


### HELPER FUNCTIONS

def get_output_polarizations(zone, zone_dict, zone_sols, zone_inputs, n_cells):
    '''given all the inputs to a zone, returns all the possible outputs'''

    # find the zone-specific indices of cells that interact with the next zone
    indices = {}
    output_indices = []
    for out_index in zone.outputs:
        output_indices.append(out_index - n_cells)

    for next_zone_key in zone.outs:
        next_zone = zone_dict[next_zone_key]
        for j in range(len(next_zone.C_ins[zone.key])):
            if any(next_zone.C_ins[zone.key][j]):
                if next_zone_key in indices:
                    indices[next_zone_key].append(j)
                else:
                    indices[next_zone_key] = [j]

    out_pols = {}
    final_pol = set()
    # retrieve the polarizations of the zone specific indices from all cases
    for zone_input in zone_inputs:
        case = zone_sols[zone_input]
        for an_out_pol in case['pstates']:
            final_pol.add(tuple(an_out_pol[j] for j in output_indices))

            for key in zone.outs:
                out_pol = tuple(an_out_pol[j] for j in indices[key])
                if key in out_pols:
                    out_pols[key].append(out_pol)
                else:
                    out_pols[key] = [out_pol]

    return out_pols, final_pol