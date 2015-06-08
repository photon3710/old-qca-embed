#!/usr/bin/python

#---------------------------------------------------------
# Name: parse_qca.py
# Purpose: New version of QCADesigner parsing
# Author:	Jacob Retallick
# Created: 04.06.2014
# Last Modified: 04.06.2015
#---------------------------------------------------------

import sys
import re

import networkx as nx
import pylab as plt
import numpy as np

#from pprint import pprint
from auxil import new_getEk

## mapping for all possible cell functions and modes

CELL_FUNCTIONS = {'QCAD_CELL_NORMAL': 0,
                  'QCAD_CELL_INPUT': 1,
                  'QCAD_CELL_OUTPUT': 2,
                  'QCAD_CELL_FIXED': 3}

CELL_MODES = {'QCAD_CELL_MODE_NORMAL': 0,
              'QCAD_CELL_MODE_CROSSOVER': 1,
              'QCAD_CELL_MODE_VERTICAL': 2,
              'QCAD_CELL_MODE_CLUSTER': 3}

## general global parameters

R_MAX = 2.5     # max cell-cell interaction range (relative to grid spacing)


def build_hierarchy(fn):
    '''Build a dict hierarchy containing all objects, their parameters, and
    childen.'''

    fp = open(fn, 'r')

    # general re expression. may need to change if future format changes
    re_start = re.compile('^\[.+\]$')
    re_term = re.compile('^\[#.+\]$')

    hier = {'label': 'Hierarchy', 'children': [], 'vars': {}}

    key_stack = ['Hierarchy']  # stack of active keys, pop of top of stack
    dict_stack = [hier]   # stack of corresponding dict objects.

    line_cnt = 0
    for line in fp:
        line_cnt += 1
        line = line.strip()  # remove endline and possible whitespace

        # must check object termination first
        if re_term.match(line):
            key = line[2:-1]
            if key_stack[-1] == key:
                d = dict_stack.pop()
                key_stack.pop()
                try:
                    dict_stack[-1]['children'].append(d)
                except:
                    print('Somehow over-popped dict_stack...')
                    return None
            else:
                print('Start-end mismatch in line {0}'.format(line_cnt))
                return None

        # for a new object, create a new dict template
        elif re_start.match(line):
            key = line[1:-1]
            key_stack.append(key)
            d = {'label': key, 'children': [], 'vars': {}}
            dict_stack.append(d)

        # otherwise check for new variable to add to most recent dict
        else:
            if '=' in line:
                var, val = line.split('=')
                dict_stack[-1]['vars'][var] = val
    fp.close()

    return hier


def proc_hierarchy(hier):
    '''Process the extracted data hierarchy to extract useful information. In
    the current information, we are interested in the overall cell grid spacing
    (for deciding on the range of included cell) and the properties of each
    cell in the circuit'''

    cells = []
    spacing = None

    # hierarchy should onlt have two children: VERSION and TYPE:DESIGN. The
    # former might be useful in later implentations for selecting formatting
    # options but for now all be care about is the DESIGN objects

    hier = [child for child in hier['children']
            if child['label'] == 'TYPE:DESIGN'][0]

    # for now assert that there can be only one cell layer, no vertical x-over
    layers = [child for child in hier['children']
              if child['label'] == 'TYPE:QCADLayer']

    # get grid spacing
    substrate = [layer for layer in layers if layer['vars']['type'] == '0'][0]
    substrate = [child for child in substrate['children']
                 if child['label'] == 'TYPE:QCADSubstrate'][0]
    spacing = float(substrate['vars']['grid_spacing'])

    # isolate cell layers
    cell_layers = [layer for layer in layers if layer['vars']['type'] == '1']

    # merge cell layers, will lead to qdot conflict if vertical x-over
    cell_dicts = [layer['children'] for layer in cell_layers]
    cell_dicts = reduce(lambda x, y: x+y, cell_dicts)

    # create cell objects
    cells = []

    for cd in cell_dicts:
        cell = {}

        # cell type
        cell['cf'] = CELL_FUNCTIONS[cd['vars']['cell_function']]
        cell['cm'] = CELL_MODES[cd['vars']['cell_options.mode']]
        cell['clk'] = int(cd['vars']['cell_options.clock'])

        # position, first child will be the QCADesignObject
        design_object = cd['children'][0]
        cell['x'] = float(design_object['vars']['x'])
        cell['y'] = float(design_object['vars']['y'])

        # quantum dots
        qdot_dicts = [child for child in cd['children']
                      if child['label'] == 'TYPE:CELL_DOT']

        qdots = []
        for d in qdot_dicts:
            dot = {}
            dot['x'] = float(d['vars']['x'])
            dot['y'] = float(d['vars']['y'])
            dot['c'] = float(d['vars']['charge'])
            qdots.append(dot)

        cell['qdots'] = qdots

        # keep track of polarization if cell is fixed: don't rely on labels
        if cell['cf'] == CELL_FUNCTIONS['QCAD_CELL_FIXED']:
            pol = qdots[0]['c']+qdots[2]['c']-qdots[1]['c']-qdots[3]['c']
            pol /= qdots[0]['c']+qdots[2]['c']+qdots[1]['c']+qdots[3]['c']
            cell['pol'] = pol

        cells.append(cell)

    return cells, spacing


def zone_cells(cells, spacing, show=False):
    '''Split cells into clock zones. Distinguishes disjoint zones with the
    same zone index'''

    N = len(cells)  # number of cells

    # construct connectivity matrix
    J = np.zeros([N, N], dtype=float)
    DR = R_MAX*spacing
    for i in xrange(N-1):
        for j in xrange(i+1, N):
            dx = cells[i]['x']-cells[j]['x']
            dy = cells[i]['y']-cells[j]['y']
            r2 = dx*dx+dy*dy
            if r2 <= DR*DR:
                Ek = new_getEk(cells[i], cells[j], DR=DR/spacing)
                J[i, j] = Ek
                J[j, i] = Ek

    # make full cell connectivity Graph
    G = nx.Graph(J)

    if show:
        plt.figure(0)
        plt.clf()
        nx.draw_graphviz(G)
        plt.show()

    # get indices for each clock index
    clk = [cell['clk'] for cell in cells]
    clk_ind = list(set(clk))    # will sort by default
    inds = [[i for i, x in enumerate(clk) if x == ind] for ind in clk_ind]

    # split graph into sub-graphs with the same clock indices
    sub_G = {ind: G.subgraph(inds[ind]) for ind in clk_ind}

    # split disconnected components for each label graph
    sub_ind = {ind: nx.connected_components(sub_G[ind]) for ind in clk_ind}

    ## find zone order

    # create abstract zone connectivity graph
    G = nx.DiGraph()
    # nodes
    for clk in clk_ind:
        for i in xrange(len(sub_ind[clk])):
            key = (clk, i)
            G.add_node(key, inds=sub_ind[clk][i])
    # edges
    for clk in clk_ind:
        adj_clk = 3 if clk == 0 else clk-1
        if not adj_clk in sub_ind:
            continue
        for i in xrange(len(sub_ind[clk])):
            k1 = (clk, i)
            for j in xrange(len(sub_ind[adj_clk])):
                k2 = (adj_clk, j)
                if np.any(J[G.node[k1]['inds'], :][:, G.node[k2]['inds']]):
                    G.add_edge(k2, k1)

    if show:
        plt.figure(1)
        plt.clf()
        nx.draw_graphviz(G)
        plt.show()

    # find input nodes, no predecessors
    predecs = {n: len(G.predecessors(n)) for n in G.nodes_iter()}
    inputs = [ky for ky, val in predecs.iteritems() if val == 0]

    # expand from inputs
    visited = {key: False for key in G.nodes()}

    nodes = inputs
    order = [nodes]
    while nodes:
        new_nodes = set()
        for node in nodes:
            new_nodes.update(G.successors(node))
            visited[node] = True
        # remove already visited nodes from new nodes
        new_nodes = [node for node in new_nodes if not visited[node]]
        nodes = new_nodes
        if nodes:
            order.append(nodes)

    # reformat order list to contain zone indices

    form_func = lambda n: sub_ind[n[0]][n[1]]
    order = [[form_func(zone) for zone in shell] for shell in order]

    return order, J


def reorder_cells(cells, zones, J, flipy=False):
    '''Renumber cells by position rather than the default QCADesigner placement
    order. Cells ordered by the tuple (zone, y, x)'''

    keys = {}

    # assign sortable tuples for each cell
    for iz in xrange(len(zones)):
        zone = zones[iz]
        for iz_sub in xrange(len(zone)):
            for ind in zone[iz_sub]:
                cell = cells[ind]
                y_sign = -1 if flipy else 1
                keys[ind] = (iz, iz_sub, y_sign*cell['y'], cell['x'])

    order = zip(*sorted([(keys[i], i) for i in keys]))[1]

    # relabel cells and reorder the J matrix
    cells = [cells[i] for i in order]
    J = J[order, :][:, order]

    # relabel each of the zones index lists
    inv_map = {order[i]: i for i in order}
    label_func = lambda lst: sorted([inv_map[c] for c in lst])
    zones = [[label_func(zn) for zn in shell] for shell in zones]

    return cells, zones, J


def parse_qca_file(fn, one_zone=False):
    '''Parse a QCADesigner file to extract cell properties. Returns an ordered
    list of cells, the QCADesigner grid spacing in nm, a list structure of the
    indices of each clock zone (propogating from inputs), and a coupling matrix
    J which contains the Ek values for cells within a radius of R_MAX times the
    grid spacing'''

    # build data hierarchy
    hier = build_hierarchy(fn)

    # extract useful information from data hierarchy
    cells, spacing = proc_hierarchy(hier)

    if one_zone:
        for cell in cells:
            cell['clk'] = 0

    # group into clock zones
    zones, J = zone_cells(cells, spacing)

    # reorder cells by zone and position
    cells, zones, J = reorder_cells(cells, zones, J)

    # if only one zone is requested, don't need the zone order structure
    if one_zone:
        zones = zones[0]

    return cells, spacing, zones, J


if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print 'No file input....'
        sys.exit()

    parse_qca_file(fn)
