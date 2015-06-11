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

import networkx as nx
import matplotlib.pyplot as plt

## PHYSICAL PARAMETERS
eps0 = 8.85412e-12  # permittivity of free space
epsr = 12.          # relative permittivity
q0 = 1.602e-19      # elementary charge

## QCADESIGNER PARSING PARAMETERS

CELL_FUNCTIONS = {'QCAD_CELL_NORMAL': 0,
                  'QCAD_CELL_INPUT': 1,
                  'QCAD_CELL_OUTPUT': 2,
                  'QCAD_CELL_FIXED': 3}

CELL_MODES = {'QCAD_CELL_MODE_NORMAL': 0,
              'QCAD_CELL_MODE_CROSSOVER': 1,
              'QCAD_CELL_MODE_VERTICAL': 2,
              'QCAD_CELL_MODE_CLUSTER': 3}

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


def comp_E_nn(spacing):
    '''compute the kink energy for two nearest interacting non-rotated cells'''

    A = 0.588672
    E_nn = A/(spacing*epsr)

    return E_nn


def prepare_convert_adj(cells, spacing, J):
    '''Prepares useful variables for converting from the parse_qca J matrix to
    a reduced adjacency form.

    outputs:    Js  : J scaled by the nearest neighbour interaction of two
                     non-rotated cells.
                T   : Array of cell-cell types for each element of J
                        1  -> non-rotated - non-rotated
                        0  -> non-rotated - rotated
                        -1 -> rotated - rotated
                DX  : X displacements in grid-spacings
                DY  : Y displacements in grid-spacings
    '''

    # scale J by the kink energy of two non-rotated adjacent cells
    E_nn = comp_E_nn(spacing)
    Js = np.round(J/E_nn, 4)

    # determine interaction type of each element of J:
    #   1  -> non-rotated - non-rotated
    #   0  -> non-rotated - rotated
    #   -1 -> rotated - rotated

    rot = [cell['rot'] for cell in cells]   # array of rotated flags
    rot = 1*np.array(rot).reshape([-1, 1])

    T = 1-(rot+rot.T)
    #T = T.astype(int)

    # get displacements between each cell

    X = np.array([cell['x'] for cell in cells]).reshape([-1, 1])
    Y = np.array([cell['y'] for cell in cells]).reshape([-1, 1])

    DX = (X.T - X)/spacing
    DY = (Y - Y.T)/spacing

    return Js, T, DX, DY


def convert_to_full_adjacency(cells, spacing, J):
    '''Convert the J matrix from parse_qca to include only full adjacency
    interactions'''

    Js, T, DX, DY = prepare_convert_adj(cells, spacing, J)


def convert_to_lim_adjacency(cells, spacing, J):
    '''Convert the J matrix from parse_qca to include only limited adjacency
    interactions'''

    Js, T, DX, DY = prepare_convert_adj(cells, spacing, J)


def construct_zone_graph(cells, zones, J, show=False):
    '''Construct a DiGraph for all the zones with keys given by (n, m) where
    n is the shell index and m is the zone index within the shell'''

    # create nodes
    G = nx.DiGraph()
    for i_shell in xrange(len(zones)):
        for i_zones in xrange(len(zones[i_shell])):
            key = (i_shell, i_zones)
            kwargs = {'inds': [], 'fixed': [], 'drivers': [], 'outputs': []}

            for ind in zones[i_shell][i_zones]:
                if cells[ind]['cf'] == CELL_FUNCTIONS['QCAD_CELL_INPUT']:
                    kwargs['drivers'].append(ind)
                elif cells[ind]['cf'] == CELL_FUNCTIONS['QCAD_CELL_FIXED']:
                    kwargs['fixed'].append(ind)
                else:
                    kwargs['inds'].append(ind)
                    if cells[ind]['cf'] == CELL_FUNCTIONS['QCAD_CELL_OUTPUT']:
                        kwargs['outputs'].append(ind)

            G.add_node(key, **kwargs)

    # edges
    for shell in xrange(1, len(zones)):
        for i in xrange(len(zones[shell-1])):
            k1 = (shell-1, i)
            for j in xrange(len(zones[shell])):
                k2 = (shell, j)
                C = J[G.node[k1]['inds'], :][:, G.node[k2]['inds']]
                if np.any(C):
                    G.add_edge(k1, k2, C=C)

    plt.figure('Zone-Graph')
    nx.draw_graphviz(G)
    plt.show()

    return G
