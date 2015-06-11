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


def comp_E_nn(spacing, OLD_QCAD = False):
    '''compute the kink energy for two nearest interacting non-rotated cells'''

    A = 0.588672
    if OLD_QCAD:
        A = 0.3827320944

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

    R_MAX = 2
    Js, T, DX, DY = prepare_convert_adj(cells, spacing, J)

    xovers = [i for i in range(len(cells)) if is_xover(cells, DX, DY, i)]
    for i in range(len(cells)):
        # check to see if the cell is involved in a cross over
        for j in range(len(cells)):
            # if it is not a cross over and futher than 2 away, strip J
            if not (i in xovers and j in xovers):
                if (DX[i][j]**2 + DY[i][j]**2 >= R_MAX**2):
                        J[i][j] = 0
                        J[j][i] = 0

    return J

def convert_to_lim_adjacency(cells, spacing, J):
    '''Convert the J matrix from parse_qca to include only limited adjacency
    interactions'''

    c_index = range(len(cells))
    R_MAX = 2
    Js, T, DX, DY = prepare_convert_adj(cells, spacing, J)

    xovers = [i for i in c_index if is_xover(cells, DX, DY, i)]
    invs = [i for i in c_index if is_inv(Js, DX, DY, i)]

    for i in c_index:
        # number of strong interactions of current cell
        si = len([j for j in c_index if Js[i][j] == 1 or Js[i][j] == -1.472])
        for j in c_index:

            dx = DX[i][j]
            dy = DY[i][j]

            if not (i in xovers and j in xovers):
                if (i in invs or j in invs) and si < 2:
                    if dx**2 + dy**2 >= R_MAX**2:
                        J[i][j] = 0
                        J[j][i] = 0

                elif dx**2 + dy**2 >= R_MAX:
                        J[i][j] = 0
                        J[j][i] = 0

    return J


def is_xover(cells, DX, DY, i):
    '''check to see if a cell is involved in a cross over'''

    #find cells directly adjacent horizontally
    hor = [j for j in range(len(DY[i])) if DY[i][j] == 0]
    x_adj = [j for j in hor if abs(DX[i][j]) == 1]

    #find cells directly adjacent vertically
    ver = [j for j in range(len(DX[i])) if DX[i][j] == 0]
    y_adj = [j for j in ver if abs(DY[i][j]) == 1]

    #if the pairs of cells are different, than there is a cross over
    if len(x_adj) == 2:
        if cells[x_adj[0]]['rot'] != cells[x_adj[1]]['rot']:
            return True

    if len(y_adj) == 2:
        if cells[y_adj[0]]['rot'] != cells[y_adj[1]]['rot']:
            return True

    # error message if somehow there is more than 2 adjacent cells in either dir
    if len(x_adj) > 2 or len(y_adj) > 2:
        print 'Error: there are %d cells horizontally adjacent' +\
            ' and %d cells vertically adjacent' % (len(x_adj), len(y_adj))

    return False

def is_inv(Js, DX, DY, i):
    '''check to see if a cell is an inverter cell
    and inverter cell is the cell that has two diagonal interactions, and one
    directly adjacent interaction (labelled IC below)
        c - c
        |    \
    c - c     IC - c
        |    /
        c - c
    '''

    index = range(len(Js[i]))
    # find number of strong and medium bonds
    m = [j for j in index if Js[i][j] == -0.2174 or Js[i][j] == 0.172]
    s = [j for j in index if Js[i][j] == 1 or Js[i][j] == -1.472]

    if len(m) >= 2 and len(s) == 1:
        # in case of weird cases - checks to see that strong interactions
        # are on the opposite side of the medium interactions
        opp = 0
        for j in m:
            if DX[i][j] == (-1) * DX[i][s[0]]:
                opp += 1
            if DY[i][j] == (-1) * DY[i][s[0]]:
                opp += 1
        return opp == 2

    return False

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
                if np.any(J[G.node[k1]['inds'], :][:, G.node[k2]['inds']]):
                    G.add_edge(k1, k2)

    plt.figure('Zone-Graph')
    nx.draw_graphviz(G)
    plt.show()

    return G
