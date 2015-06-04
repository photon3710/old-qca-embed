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

from pprint import pprint

## mapping for all possible cell functions and modes

CELL_FUNCTIONS = {'QCAD_CELL_NORMAL': 0,
                  'QCAD_CELL_INPUT': 1,
                  'QCAD_CELL_OUPUT': 2,
                  'QCAD_CELL_FIXED': 3}

CELL_MODES = {'QCAD_CELL_MODE_NORMAL': 0,
              'QCAD_CELL_MODE_CROSSOVER': 1,
              'QCAD_CELL_MODE_VERTICAL': 2,
              'QCAD_CELL_MODE_CLUSTER': 3}


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
    layers = [child for child in hier['children']]

    # get grid spacing
    substrate = [layer for layer in layers if layer['vars']['type'] == '0'][0]
    spacing = float(substrate['vars']['type'])

    # isolate cell layers
    cell_layers = [layer for layer in layers if layer['vars']['type'] == '1']

    # merge cell layers
    cell_dicts = [layer['children'] for layer in cell_layers]
    cell_dicts = reduce(lambda x, y: x+y, cell_dicts)

    # create cell objects
    cells = []

    for cd in cell_dicts:
        cell = {}

        # cell type
        cell['cf'] = CELL_FUNCTIONS[cd['vars']['cell_function']]
        cell['cm'] = CELL_MODES[cd['vars']['cell_options.mode']]

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

        pprint(cell)

    return cells, spacing


def parse_qca_file(fn):
    ''' '''

    # build data hierarchy
    hier = build_hierarchy(fn)

    # extract useful information from data hierarchy
    cells, spacing = proc_hierarchy(hier)


if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print 'No file input....'
        sys.exit()

    parse_qca_file(fn)
