#!/usr/bin/python

#---------------------------------------------------------
# Name: solution.py
# Purpose: Container file for Solution class and associated functions
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 12.06.2015
#---------------------------------------------------------

import xml.etree.ElementTree as ET
from zone import write_zones_to_xml, read_zones_from_xml
import networkx as nx
import numpy as np

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
            e = str_to_2d_tuple(edge.text)
            Gz.add_edge(e[0], e[1])
        for node in root.iter('node'):
            Gz.add_node(str_to_tuple_int(node.text))
        self.Gz = Gz

        self.zones = read_zones_from_xml(fn + ZONE_EXT)

        z_o = root.find('z_orders')
        self.z_orders = str_to_2d_list_of_tup(z_o.text)

        sols = []

        for zone in root.iter('zone'):
            zone_dict = {}
            for case in zone.iter('case'):
                case_dict = {}
                key = str_to_2d_tuple(case.find('case_key').text)

                #states
                case_dict['states'] = str_to_list_of_nparr(case.find('states').text)
                #pstates
                case_dict['pstates'] = str_to_list_of_tup(case.find('pstates').text)
                #state_pols
                case_dict['state_pols'] = str_to_list_of_nparr(case.find('state_pols').text)
                #Eps
                case_dict['Eps'] = str_to_tuple_float(case.find('Eps').text)
                #Es
                case_dict['Es'] = str_to_nparr(case.find('Es').text)

                zone_dict[key] = case_dict
            sols.append(zone_dict)
        self.sols = sols

    def get_inputs(self):
        '''Get list of all input cell indices'''
        input_indices = []
        for zone in self.zones:
            input_indices.append(zone.drivers)
        return input_indices

    def run_input_single(self, pols):
        '''Use the solution information to deterime the ouput polarizations
        for a single set of input polarizations'''
        pass

    def run_input_sequence(self, pol_seq):
        '''Use the solution information to deterime the output sequence for a
        sequence of input sets (relevant for circuits with feedback'''
        pass

### HELPER FUNCTIONS

def str_to_list(string):
    '''inverse of str(list[])'''
    l = []
    for i in string.strip('[]').split(','):
        if i != '':
            l.append(int(i.strip(' ')))
    return l

def str_to_nparr(string):
    '''inverse of str(ndarray)'''
    l = []
    for i in string.strip('[]').split(' '):
        if i != '':
            l.append(float(i))
    return np.asarray(l)

def str_to_list_of_tup(string):
    '''inverse of str([(,),(,)]'''
    arr = []
    for i in string.strip('[]').split('),'):
        l = []
        if i.strip('( )') == '':
            l.append(tuple())
        else:
            for j in i.strip('( )').split(','):
                l.append(int(j))
        arr.append(tuple(l))
    return arr

def str_to_2d_list_of_tup(string):
    '''inverse of str([[(,),(,)],[(,)...]])'''
    arr = []
    for i in string.split('],'):
        l = []
        if i.strip('[ ]') == '':
            l.append(tuple())
        else:
            for j in i.strip(',[ ]').split('),'):
                l.append(str_to_tuple_int(j.strip('[ ]')))
        arr.append(l)
    return arr

def str_to_2d_np_arr(string):
    '''inverse function of str(np.array)'''
    arr = []
    for i in string.split('['):
        l = []
        if i.strip(' ') != '':
            for j in i.strip(',[ ]').split(','):
                l.append(float(j.strip('[ ]')))
            arr.append(l)
    return np.asarray(arr)

def str_to_tuple_int(string):
    '''inverse operation of str(tuple)'''
    arr = string.strip('()').split(',')
    l = []
    for i in arr:
        l.append(int(i))
    return tuple(l)

def str_to_tuple_float(string):
    '''inverse operation of str(tuple)'''
    arr = string.strip('()').split(',')
    l = []
    for i in arr:
        l.append(float(i))
    return tuple(l)

def str_to_2d_tuple(string):
    '''inverse operation of str(tuple of tuples)'''
    arr = string.strip('()').split('),')
    l = []
    for i in arr:
        if i != '':
            i = i.strip('( )').split(',')
            il = []
            for j in i:
                if j != '':
                    il.append(int(j.strip('( )')))
            l.append(tuple(il))
    return tuple(l)

def str_to_list_of_nparr(string):
    '''inverse of str([nd.array])'''
    arr = string.strip('[]').split('array([')
    l = []
    for a in arr:
        a = a.strip(']), ')
        npl = []
        if a != '':
            a = a.split(',')
            for n in a:
                if n != '':
                    npl.append(float(n))
        if npl:
            l.append(np.asarray(npl))
    return l
