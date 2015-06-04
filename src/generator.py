#!usr/bin/python

#---------------------------------------------------------
# Name: generator.py
# Purpose: QCA circuit generator
# Author:    Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

from __future__ import division
import networkx as nx
import pylab as plt
import numpy as np
from random import randint, random, shuffle

WIRE_MAX = 15        # maximum length of random wire
INIT_SEED = [1, 3]    # range of number of components initially created
NEW_SEED = [1, 2]    # range of number of components from a single output

MAX_COUNT = 200
MAX_COMP = 14


def GtoCoef(G):
    '''Convertes a graph representation of a circuit to h and J coefficients'''

    N = G.number_of_nodes()

    h = []
    J = np.zeros([N, N], dtype=float)

    nodes = G.nodes()
    nodes_key = {nodes[i]: i for i in xrange(len(nodes))}

    for node in nodes:
        h.append(0)
        if 'h' in G.node[node]:
            h[-1] = G.node[node]['h']
        for key in G.edge[node]:
            n, m = sorted([nodes_key[node], nodes_key[key]])
            if n == m:
                continue
            J[n, m] = G.edge[node][key]['J']

    return h, J


def createInv(index, full_adj=True):
    '''
    Creates an inverter subgraph with nodes with names starting at index.
    Returns the subgraph G, the list of input nodes In, the list of output
    nodes Out, and the next index.
    '''

    G = nx.Graph()
    G.add_nodes_from([index+i for i in xrange(7)])

    # -1 interactions
    edges = [[0, 1], [1, 2], [1, 3], [2, 4], [3, 5]]
    edges = map(lambda x: map(lambda y: y+index, x), edges)
    G.add_edges_from(edges, J=-1)

    # 0.177 interactions
    if full_adj:
        edges = [[0, 2], [0, 3], [1, 4], [1, 5], [4, 6], [5, 6]]
    else:
        edges = [[4, 6], [5, 6]]
    edges = map(lambda x: map(lambda y: y+index, x), edges)
    G.add_edges_from(edges, J=0.177)

    In = [index]
    Out = [index+6]

    return G, In, Out, index+7


def createMAJ(index, full_adj=True):
    '''
    Creates a majority gate subgraph with nodes with names starting at index.
    Returns the subgraph G, the list of input nodes In, the list of output
    nodes Out, and the next index.
    '''

    G = nx.Graph()
    G.add_nodes_from([index+i for i in xrange(5)])

    # -1 interactions
    edges = [[0, 3], [1, 3], [2, 3], [3, 4]]
    edges = map(lambda x: map(lambda y: y+index, x), edges)
    G.add_edges_from(edges, J=-1)

    # 0.177 interactions
    if full_adj:
        edges = [[0, 1], [0, 4], [1, 2], [2, 4]]
        edges = map(lambda x: map(lambda y: y+index, x), edges)
        G.add_edges_from(edges, J=0.177)

    In = [index+i for i in xrange(3)]
    Out = [index+4]

    return G, In, Out, index+5


def routeWire(G, start, end, index, N=None):
    '''Adds nodes to graph G between labelled start and end nodes. If N is
    given, adds a wire of length N (N>=0), otherwise N is a random integer
    between 0 and WIRE_MAX. Returns next index'''

    # check start and end in G
    if not (G.has_node(start) and G.has_node(end)):
        print 'Invalid target nodes...'
        return -1

    if N is None:
        N = randint(0, WIRE_MAX)

    if N == 0:
        G.add_edge(start, end, J=-1)
        return index

    edges = [[start, index]]+[[index+i, index+i+1]
                              for i in xrange(0, N-1)] + [[index+N-1, end]]
    G.add_edges_from(edges, J=-1)

    return index+N


def routeWires(G, start, ends, index):
    '''Routes a single output to a number of inputs'''

    if not ends:
        return index

    if len(ends) == 1:
        return routeWire(G, start, ends[0], index)

    shuffle(ends)
    # generate split point, split target end points
    SPLIT = randint(1, len(ends)-1)
    ends_l = ends[0:SPLIT]
    ends_r = ends[SPLIT::]

    # add split point to Graph
    split_node = index
    G.add_node(split_node)
    index+1

    # route connection

    index = routeWire(G, start, split_node, index)

    # build left/right routes

    index = routeWires(G, split_node, ends_l, index)
    index = routeWires(G, split_node, ends_r, index)

    return index


def newComponent(index, full_adj=True):
    ''' generates a new component '''

    r = random()

    if r < .5:
        G, In, Out, index = createInv(index, full_adj=full_adj)
        typ = 'inv'
    else:
        G, In, Out, index = createMAJ(index, full_adj=full_adj)
        typ = 'maj'

    return G, In, Out, index, typ


def pickM(mlow, mhigh):
    pm = [1./x for x in xrange(mlow, mhigh+1)]
    s = sum(pm)
    pm = map(lambda x: x/s, pm)

    r = random()
    p = 0.
    for i in xrange(len(pm)-1):
        p += pm[i]
        if r < p:
            return mlow+i
    else:
        return mhigh


def generateCircuit(max_count=MAX_COUNT, full_adj=True):
    '''Generates a random QCA circuit'''

    inputs = []
    outputs = []

    G = nx.Graph()
    index = 0
    TYP_COUNT = {'inv': 0, 'maj': 0}

    # initial seed

    N = pickM(*INIT_SEED)
    for n in xrange(N):
        S, In, Out, index, typ = newComponent(index, full_adj=full_adj)
        TYP_COUNT[typ] += 1
        G = nx.union(G, S)
        inputs += In
        outputs += Out

    # generator loop

    go = True

    while go and index < MAX_COUNT:

        # outputs should not be empty

        if inputs and random() < .3:    # fill an input
            # select input
            i = randint(0, len(inputs)-1)
            G.node[inputs[i]]['h'] = (2*random()-1)
            inputs.pop(i)
        else:
            # select output
            i = randint(0, len(outputs)-1)
            out = outputs[i]

            # number of new outputs
            M = pickM(*NEW_SEED)

            for m in xrange(M):
                # route output to input
                chk = (len(inputs)/(len(inputs)+len(outputs)))
                if inputs and random() < chk:
                    i2 = randint(0, len(inputs)-1)
                    index = routeWire(G, out, inputs.pop(i2), index)
                else:    # create new component
                    S, In, Out, index, typ = newComponent(index,
                                                          full_adj=full_adj)
                    TYP_COUNT[typ] += 1
                    G = nx.union(G, S)
                    # select input
                    i2 = randint(0, len(In)-1)
                    index = routeWire(G, out, In.pop(i2), index)
                    inputs += In
                    outputs += Out

        if random() <= (1./(2*len(outputs)+1)):
            go = False

    # fill available inputs

    print 'Unrouted inputs: %d' % len(inputs)

    while inputs:
        inp = inputs.pop()
        G.node[inp]['h'] = random()

    print 'Unrouted outputs: %d' % len(outputs)

    print 'Inverters Used: %d' % TYP_COUNT['inv']
    print 'Majority Gates Used: %d' % TYP_COUNT['maj']

    # keep only largest connected subgraph

    if not nx.is_connected(G):
        print 'Generate Graph not connected... selecting largest subgraph'
        Gs = sorted(nx.connected_component_subgraphs(G),
                    key=lambda x: x.number_of_nodes())
        G = Gs[-1]

    print 'Number of Cells: %d' % G.number_of_nodes()
    return G


def generateCircuit2(full_adj=True):
    '''Generate a QCA circuit by randomly generating inverters and majority
    gates and then connecting outputs to inputs'''

    G = nx.Graph()
    index = 0
    inputs = []
    outputs = []

    # generate components
    N_inv = randint(0, MAX_COMP)
    N_maj = randint(0 if N_inv else 1, MAX_COMP-N_inv)

    for n in xrange(N_inv):
        S, In, Out, index = createInv(index, full_adj=full_adj)
        G = nx.union(G, S)
        inputs += In
        outputs += Out

    for n in xrange(N_maj):
        S, In, Out, index = createMAJ(index, full_adj=full_adj)
        G = nx.union(G, S)
        inputs += In
        outputs += Out

    print 'Number of inverters placed: %d' % N_inv
    print 'Number of majority gates placed: %d' % N_maj
    print 'Inputs : %s' % str(inputs)
    print 'Outputs : %s' % str(outputs)
    ## route outputs to inputs

    # shuffle outputs

    shuffle(outputs)

    # there will be 2*N_maj more inputs than outputs
    # set up the probabilities such that we should expect all the inputs to
    # connect to an output.

    ratio = len(inputs)/len(outputs)

    for output in outputs:
        if not inputs:
            break
        p = ratio/(len(inputs)+1)
        Ins = []
        for i in xrange(len(inputs)):
            if random() < p:
                Ins.append(inputs[i])

        index = routeWires(G, output, Ins, index)
        for In in Ins:
            inputs.remove(In)

    # final graph

    print 'Number of free inputs: %d' % len(inputs)

    for inp in inputs:
        G.node[inp]['h'] = random()

    # keep only largest connected subgraph

    if not nx.is_connected(G):
        print 'Generate Graph not connected... selecting largest subgraph'
        Gs = sorted(nx.connected_component_subgraphs(G),
                    key=lambda x: x.number_of_nodes())
        G = Gs[-1]

    print 'Number of Cells: %d' % G.number_of_nodes()
    return G


if __name__ == '__main__':
    G = generateCircuit()
    print G.number_of_edges()
    print G.number_of_nodes()
    nx.draw_shell(G)
    plt.show()
    h, J = GtoCoef(G)
    print len(h)
    print J.shape
