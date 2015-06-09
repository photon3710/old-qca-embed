#!/usr/bin/python

#---------------------------------------------------------
# Name: auxil.py
# Purpose: Auxiliary commonly used functions
# Author:	Jacob Retallick
# Created: 09.06.2014
# Last Modified: 09.06.2015
#---------------------------------------------------------

import numpy as np
import sys


## PHYSICAL PARAMETERS
eps0 = 8.85412e-12  # permittivity of free space
epsr = 12.          # relative permittivity
q0 = 1.602e-19      # elementary charge

### GENERAL FUNCTIONS


def pinch(string, pre, post):
    '''selects the string between the first instance of substring (pre) and
    the last instance of substring (post)'''
    return string.partition(pre)[2].rpartition(post)[0]


### QCA CELL PROCESSING


def getEk(c1, c2, DR=2):
    '''Compute the kink energy using the qdot positions (in nm) of two cells.
    If the cell displacement is greater than DR return False'''

    # check cell-cell range
    dx = c1['x']-c2['x']
    dy = c1['y']-c2['y']
    if dx*dx+dy*dy > DR*DR:
        return False

    qdots_1 = c1['qdots']
    qdots_2 = c2['qdots']

    # compute displacements

    x1 = [qd['x'] for qd in qdots_1]
    y1 = [qd['y'] for qd in qdots_1]

    x2 = [qd['x'] for qd in qdots_2]
    y2 = [qd['y'] for qd in qdots_2]

    X1 = np.array([x1, y1]).T.reshape([4, 1, 2])
    X2 = np.array([x2, y2]).T.reshape([1, 4, 2])

    R = np.sqrt(np.sum(pow(X1-X2, 2), axis=2))

    if np.min(R) == 0:
        print 'qdot overlap detected'
        sys.exit()

    # QCADesigner orders qdots either CW or CCW so same and diff configurations
    # are always alternating indices.

    Q = np.array([1, -1, 1, -1])    # template for charge arrangement

    Q = np.outer(Q, Q)

    Ek = -1e9*q0*np.sum((Q/R))/(8*np.pi*eps0*epsr)

    return Ek
    

def convert_to_full_adjacency(cells, zones, J):
    '''Convert the J matrix from parse_qca to include only full adjacency
    interactions'''
    
    rot_flags = [cell['rot'] for cell in cells] # flags for rotated cells
        
    


def convert_to_lim_adjacency(cells, zones, J):
    '''Convert the J matrix from parse_qca to include only limited adjacency
    interactions'''
    pass