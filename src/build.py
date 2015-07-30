#!usr/bin/python

#---------------------------------------------------------
# Name: build.py
# Purpose: Scripts for generating basic QCA circuits
# Author:    Jacob Retallick
# Created: 02.08.2014
# Last Modified: 07.05.2015
#---------------------------------------------------------

import numpy as np
ADJ_RADIUS = 1.5


def JtoR(J, ds):
    ''' converts a coupler coefficient to a [dx,dy] pair.'''

    if J == 0:
        return [1.1*ADJ_RADIUS, 0]

    theta = 0 if J < 0 else np.pi/4

    r5 = -np.cos(4*theta)/J
    r = pow(r5, .2)

    dx = r*np.cos(theta)
    dy = r*np.sin(theta)

    return [dx*ds, dy*ds]


def addInv(X, Y, spacing, BRANCH=2):
    ''' add an inverter right of the cell with position given by X[-1]
     and Y[-1]. BRANCH is the length of the cell chains in the inverter
     side branch (should be at least 2) '''

    X0, Y0 = X[-1]+spacing, Y[-1]
    X.append(X0)
    Y.append(Y0)

    for i in xrange(BRANCH):
        # lower branch
        X.append(X0+i*spacing)
        Y.append(Y0-spacing)
        # upper branch
        X.append(X0+i*spacing)
        Y.append(Y0+spacing)

    # inverted cell
    X.append(X0+BRANCH*spacing)
    Y.append(Y0)


def createWire(N, C, P, C_int=-1):
    ''' Allocate cell positions and driver polarisations for a QCA wire.

    inputs:    N -> list(int)    : lengths of cell groups with i'th group
                              interacting with (i+1)'th group to give
                              C[i+1]
            C -> list(float): normalised coupling strengths between groups
            P -> list(float): polarizations of driver cells. If P[1] is None
                              the 'right' driver is not included
            C -> float        : internal coupling strength in wire segments

    example:    N=[1,2,4],C=[.177, -.5, -.8, -.2], P=[1,-1]

            -    The wire has a driver on the left with polarisation +1
            -    Right of the driver is a single cell offset to give Ek~.177
            -    Offset of this cell is a group of 2 cells (coupled by -1)
                with offset such that Ek~-.5
            -    Again for a 4 cell group with offset -> Ek~-.8
            -    As P[1] is not None, there is a driver at the end of the
                wire with polarisation -1 offset to give Ek~.2

            --    if P[1] were None, C should not include the last value and
                the driver would not be present
    '''

    ## check/assert valid inputs

    try:
        P = [float(P), None]    # P was a single value (only left driver)
    except:
        # P is an iterable of some sort
        pass

    if len(P) > 2:  # truncate driver polarisation list
        print 'Truncated Polarisation list'
        P = P[:2]

    ND = 1 if (P[1] is None) else 2    # number of drivers

    # N cannot be empty (otherwise no wire to create)

    if len(N) == 0:
        print 'N is empty... no wire to create'
        return [], -1

    # length of C should be len(N)+ND-1

    if len(C) != (len(N)+ND-1):
        print 'Invalid number of coupling strengths'
        return [], -1

    ## inputs good: create cells

    cells = []
    spacing = 20

    ## cell positions

    J = []
    for i in xrange(len(N)):
        J += [C[i]]+[C_int]*(N[i]-1)

    if ND == 2:
        J.append(C[-1])

    R = [JtoR(j, spacing) for j in J]

    X = [0]
    Y = [0]

    for dx, dy in R:
        X.append(dx+X[-1])
        Y.append(dy+Y[-1])

    for i in xrange(len(X)):

        cell = {}
        cell['x'] = X[i]
        cell['y'] = Y[i]
        cell['cx'] = cell['cy'] = .9*spacing
        cell['number'] = i
        cell['type'] = 0

        cells.append(cell)

    cells[0]['pol'] = P[0]
    cells[0]['type'] = 2

    if ND == 2:
        cells[-1]['pol'] = P[1]
        cells[-1]['type'] = 2

    return cells, spacing


def createInv(N, P):
    ''' Allocate cell list for an inverter chain.

    inputs: N -> list(int)     : list of cell chain lengths between each
                              inverter group (N[0] from driver to first)
            P -> float        : polarisation of driver cell
    '''

    ## check/assert valid inputs

    try:
        N = [float(N)]
    except:
        pass

    try:
        P = float(P)
    except:
        try:
            P = float(P[0])
        except:
            print 'Invalid P format'
            return [], -1

    cells = []
    spacing = 20

    X = [0]
    Y = [0]

    for i in xrange(len(N)):
        # wire portion
        for j in xrange(N[i]-(1 if i else 0)):
            X.append(spacing+X[-1])
            Y.append(Y[-1])
        # inverter
        addInv(X, Y, spacing)

    for i in xrange(len(X)):

        cell = {}
        cell['x'] = X[i]
        cell['y'] = Y[i]
        cell['cx'] = cell['cy'] = .9*spacing
        cell['number'] = i
        cell['type'] = 0

        cells.append(cell)

    cells[0]['pol'] = P
    cells[0]['type'] = 2

    return cells, spacing


def createMAJ(N, P):
    ''' Allocate cell list for a Majority gate.

    inputs: N -> list(int)     : list of chain lengths for inputs to the gate
            P -> list(float): polarisation of inputs

    note: indexing goes [top,left,bottom]
    '''

    if len(N) != 3:
        print 'Invalid chain length list'
        return [], -1

    if not all(map(lambda x: x >= 1, N)):
        print 'Each input must have length >= 1'
        return [], -1

    if len(P) != 3:
        print 'Invalid driver polarisation list'
        return [], -1

    X = []
    Y = []

    cells = []
    ds = spacing = 20

    # top branch
    X0, Y0 = N[1]*ds, N[0]*ds
    for i in xrange(N[0]):
        X.append(X0)
        Y.append(Y0-i*ds)

    # left branch
    X0, Y0 = 0, 0
    for i in xrange(N[1]):
        X.append(X0+i*ds)
        Y.append(Y0)

    # bottom branch
    X0, Y0 = N[1]*ds, -N[2]*ds
    for i in xrange(N[2]):
        X.append(X0)
        Y.append(Y0+i*ds)

    for i in xrange(2):
        X.append(X0+i*ds)
        Y.append(0)

    for i in xrange(len(X)):

        cell = {}
        cell['x'] = X[i]
        cell['y'] = Y[i]
        cell['cx'] = cell['cy'] = .9*spacing
        cell['number'] = i
        cell['type'] = 0

        cells.append(cell)

    DRIVERS = [0, N[0], N[0]+N[1]]

    for i in xrange(3):
        d = DRIVERS[i]
        cells[d]['pol'] = P[i]
        cells[d]['type'] = 2

    INPUTS = [N[0]-2, N[0]+N[1]-3, N[0]+N[1]+N[2]-4]  # after removeing driver

    print 'Majority gate cell key:'
    #print 'Driver Cells: ' + ', '.join(map(str,DRIVERS))
    print 'Input Cells: ' + ', '.join(map(str, INPUTS))
    print 'Output Cell: ' + str(cells[-1]['number']-3)

    return cells, spacing
