#!usr/bin/python

#---------------------------------------------------------
# Name: auxil.py
# Purpose: Auxiliary commonly used functions
# Author:	Jacob Retallick
# Created: 02.08.2014
# Last Modified: 07.05.2015
#---------------------------------------------------------

import sys
import numpy as np
import scipy.sparse as sp

ADJ_RADIUS = 1.5    # distance to separate full from limited adjacency
TYPEMAP = {'QCAD_CELL_NORMAL': 0,
           'QCAD_CELL_OUTPUT': 1,
           'QCAD_CELL_FIXED': 2,
           'QCAD_CELL_INPUT': 3}

# physical parameters
EK0 = 300e-3        # kink energy in eV
EK_POW = 5          # fall-off power of the kink energy

GAMMA = 0.1         # default tunneling energy, relative to EK

NEAREST_FACT = .5   # portion of EK0 for 'nearest neighbour' coupling


### GENERAL FUNCTIONS


def pinch(string, pre, post):
    '''selects the string between the first instance of substring (pre) and
    the last instance of substring (post)'''
    return string.partition(pre)[2].rpartition(post)[0]


### QCA CELL PROCESSING


def getEk(c1, c2, spacing):
    ''' returns the theoretical value of Ek for two cells'''

    dx = abs(c1['x']-c2['x'])
    dy = abs(c1['y']-c2['y'])

    #print 'DX: ' + str(dx) + '\t DY: ' + str(dy)

    if dx+dy < .5*spacing:
        n1 = str(c1['number'])
        n2 = str(c2['number'])
        print 'Detected cell overlap... cells: ' + n1 + ' and ' + n2
        sys.exit()

    r = np.sqrt(dx*dx+dy*dy)/spacing

    if abs(dx/spacing) < 0.01:    # vertical alignment (divide by zero)
        theta = np.pi/2
    else:
        theta = np.arctan(dy/dx)

    if r > ADJ_RADIUS:
        return False
    else:
        return EK0*np.cos(4*theta)/pow(r, EK_POW)


def generateAdjDict(cells, spacing, verbose=False):
    '''generates a dict which has list elements containing tuples of type
    (cell_index, Ek). For each cell in [cells], finds the indexes of all
    cells satisfying 'adjacency' and appends these indices with the
    corresponding Ek values. First element of adjacency[i] is h_i '''

    N = len(cells)
    adjacency = {}

    drivers = []

    # separate drivers from other cells
    driver_index = []
    for i in xrange(N):
        cell = cells[i]
        if cell['type'] == TYPEMAP['QCAD_CELL_FIXED']:
            drivers.append(cell)
            driver_index.append(i)
        else:
            adjacency[i] = [0.]    # initial h value

    if verbose:
        print 'Driver cell indentified: ' + str(driver_index)

    # assert adjacency key order
    if verbose:
        print 'Cell order: ' + str(adjacency.keys())
        order = sorted(adjacency.keys())
        print 'Sorted order: ' + str(order)
    else:
        order = adjacency.keys()

    M = len(order)
    for i in xrange(M):
        cell = cells[order[i]]

        # find driver contribution
        for driver in drivers:
            Ek = getEk(cell, driver, spacing)
            if Ek is False:
                continue
            else:
                try:
                    pol = driver['pol']
                except KeyError:
                    print 'Driver does not have a polarization... setting as 0'
                    pol = 0.
                adjacency[order[i]][0] += Ek*pol

        # find adjacencies
        adj = []
        for j in xrange(i+1, M):
            cell2 = cells[order[j]]
            Ek = getEk(cell, cell2, spacing)
            if Ek is False:
                continue
            else:
                adj.append((cell2['number'], Ek))
        adjacency[order[i]] += adj

    return adjacency, driver_index


def detectComponents(adj_dict):
    '''Determine the number of majority gates and inverters. Driver cells not
    counted.'''

    cells = adj_dict.keys()

    majs = []
    invs = []

    num_nearest = {cell: 0 for cell in cells}
    num_next = {cell: 0 for cell in cells}

    # count number of nearest and next nearest nieghbours
    for cell in cells:
        adj = adj_dict[cell][1:]
        for c2, Ek in adj:
            if Ek > 0.8 * EK0:
                num_nearest[cell] += 1
                num_nearest[c2] += 1
            elif Ek < 0 and Ek > -0.4:
                num_next[cell] += 1
                num_next[c2] += 1

    for cell in cells:
        if num_nearest[cell] == 4:
            majs.append(cell)
        elif num_nearest[cell] == 1 and num_next[cell] == 2:
            invs.append(cell)

    return majs, invs


def convertToNearestNeighbour(adjacency, driver_index=None):
    '''Converts the output from generateAdjDict to that for restricted
    nearest neighbour coupling'''

    ## count 'nearest neighbour' cells: directly adjacent

    cell_map = adjacency.keys()
    N = len(cell_map)

    # recall adjacency[cell] is a list of tuple pairs (cell_index, Ek) with
    # adj[cell][0] = h_{cell_index}

    num_nearest = [0]*N
    num_next = [0]*N
    ignore = []
    nearest_thresh = NEAREST_FACT*EK0

    # generate full connectivity list
    connect = {i: [] for i in xrange(N)}
    for i in xrange(N):
        adj = list(adjacency[cell_map[i]])
        h = adj.pop(0)
        connect[i].insert(0, h)
        for cell_index, Ek in adj:
            ind = cell_map.index(cell_index)
            connect[i].append([ind, Ek])
            connect[ind].append([i, Ek])

    # count the number of nearest neighbour cells for each cell unless
    for i in xrange(N):
        con = list(connect[i])
        h = con.pop(0)
        for ind, Ek in con:
            if Ek > nearest_thresh:
                num_nearest[i] += 1

    # count the number of next nearest neighbouts
    for i in xrange(N):
        con = list(connect[i])
        h = con.pop(0)
        for ind, Ek in con:
            if abs(Ek) < nearest_thresh:
                num_next[i] += 1
            if num_nearest[ind] > 1:
                ignore.append(i)

    # find inverter cells
    invs = []
    for i in xrange(N):
        if num_nearest[i] == 1 and num_next[i] == 2 and not i in ignore:
            invs.append(cell_map[i])

    ##     if only one nearest neighbour, include next nearest coupling
    ##    else use only nearest neighbours

    new_adjacency = {}

    for i in xrange(N):
        key = cell_map[i]
        if key in invs:
            new_adjacency[key] = adjacency[key]
        else:
            new_adjacency[key] = [adjacency[key].pop(0)]
            for cell_index, Ek in adjacency[key]:
                if Ek > nearest_thresh or cell_index in invs:
                    new_adjacency[key].append([cell_index, Ek])

    return new_adjacency


def adjToCoef(adjacency, verbose=False):
    ''' convert adjacency matrix to h and J coefficients'''

    cell_map = adjacency.keys()
    N = len(cell_map)

    h = []
    J = np.zeros([N, N], dtype=float)

    for i in xrange(N):
        adj = adjacency[cell_map[i]]
        h.append(adj.pop(0))
        for cell_index, Ek in adj:
            J[i, cell_map.index(cell_index)] = -Ek

    if verbose:
        print '%d non driver cells detected' % len(h)

    return h, J


### HAMILTONIAN GENERATION

PAULI = {}
PAULI['x'] = sp.dia_matrix([[0, 1], [1, 0]])
#PAULI['y'] = sp.dia_matrix([[0, 1j], [-1j, 0]])
PAULI['z'] = sp.dia_matrix([[-1, 0], [0, 1]])


def stateToPol(state):
    '''converts a 2^N size state vector to an N size polarization vector'''

    state = np.asmatrix(state)

    ## correct state alignment

    a, b = state.shape
    if a < b:
        state = state.transpose()

    amp = np.abs(np.multiply(state.conjugate(), state))

    N = int(np.log2(np.size(state)))
    SZ = [np.asmatrix(pauli(i+1, N, 'z')) for i in xrange(N)]

    POL = [-float(SZ[i]*amp) for i in xrange(N)]
    return POL


def pauli(index, N, typ):
    '''computes the tensor product sigma_typ(i) '''

    if index < 1 or index > N:
        print 'Invalid tensor product index...must be in range (1...N)'
        sys.exit()

    if not typ in PAULI.keys():
        print "Invalid pauli matrix type... must be in ['x', 'y', 'z']"
        sys.exit()

    p_mat = PAULI[typ]

    if index == 1:
        product = sp.kron(p_mat, sp.eye(pow(2, N-index)))
    else:
        temp = sp.kron(sp.eye(pow(2, index-1)), p_mat)
        product = sp.kron(temp, sp.eye(pow(2, N-index)))

    if typ == 'z':
        return product.diagonal()
    else:
        return sp.tril(product)


def generateHam(h, J, gamma=None):
    ''' computes the Hamiltonian as a sparse matrix.

    inputs:    h - iterable of size N containing on site energies
            J - matrix (type J[i,j]) containing coupling strengths. Needs
                to contain at least the upper triangular values
    '''

    N = len(h)
    N2 = pow(2, N)

    # initialise pauli and data matrices

    #J=sp.triu(J)

    offdiag_flag = False
    if type(gamma) in [int, float]:
        offdiag_flag = True
        gamma = [gamma]*N
        SX = [pauli(i+1, N, 'x') for i in xrange(N)]  # tril sp mat format
        OFFDIAG = sp.dia_matrix((N2, N2), dtype=float)

    SZ = [pauli(i+1, N, 'z') for i in xrange(N)]  # diag np.array

    DIAG = np.zeros([1, N2], dtype=float)

    for i in xrange(N):

        DIAG += SZ[i]*h[i]

        for j in xrange(i+1, N):
            DIAG += J[i, j]*np.multiply(SZ[i], SZ[j])

        if offdiag_flag:
            OFFDIAG = OFFDIAG - gamma[i]*SX[i]

    H = sp.diags(DIAG[0], 0)
    if offdiag_flag:
        H = H+OFFDIAG

    upper = sp.tril(H, -1).getH()

    return H+upper


#############################################################333
## FORMATTING FUNCTIONS


def coefToConn(h, J):
    '''convert the h and J coefficients into a full adjacency list
    for embedding, 0->N indexing '''

    N = len(h)

    D = {i: [] for i in xrange(N)}

    for i in xrange(N):
        d = list(J[i].nonzero()[0])
        for j in d:
            D[i].append(j)
            D[j].append(i)

    return D
