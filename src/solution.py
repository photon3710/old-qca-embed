#!/usr/bin/python

#---------------------------------------------------------
# Name: solution.py
# Purpose: Container file for Solution class and associated functions
# Author:	Jacob Retallick
# Created: 10.06.2014
# Last Modified: 12.06.2015
#---------------------------------------------------------


class Solution:
    '''A Solution object is constructed from the results of running a solver on
    a QCA circuit. It is assumed that the solver approach is to find the
    spectra of each clock zone for all possible inputs. The Solution object
    contains methods to analyse these results in order to determine propoerties
    of the whole circuit.'''

    def __init__(self):
        '''Initialise a Solution object'''

        self.zones = []
        self.z_orders = []
        self.sols = []

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
