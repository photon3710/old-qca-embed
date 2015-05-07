#!usr/bin/python

#---------------------------------------------------------
# Name: dense_test.py
# Purpose: Test code for the dense placement algorithm
# Author:    Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

import sys
#from time import time
from auxil import pinch

ADJ_RADIUS = 1.5

flagchars = ['[', ']']

keywords = {'cell': 'QCADCell',
            'init': 'TYPE',
            'term': '#TYPE',
            'delim': ':'}

header_term = '#VERSION'

TYPEMAP = {'QCAD_CELL_NORMAL': 0,
           'QCAD_CELL_OUTPUT': 1,
           'QCAD_CELL_FIXED': 2,
           'QCAD_CELL_INPUT': 3}

SPACING = 'grid_spacing'
KILL_ON_NO_POL = False


def buildHierarchy(data):
    '''Constructs a hierarchy of dicts describing the objects in a QCAD
    file. Objects contained within others (for example a QCADDesignObject
    within a QCACell) are appended to the children of their container objects.
    Parameters are included within objects with keys given by their names as
    they appear in the file'''

    hierarchy = []
    stack = []
    readflag = False
    cells = []

    spacing = 20.

    for line in data:

        try:
            if SPACING in line:
                spacing = float(line.rpartition('=')[2].strip())
                continue

            if False:
                print '>>' + line

            # ignore header
            if not readflag:
                if header_term in line:
                    readflag = True
                continue

            # object control
            if flagchars[0] in line:  # start or complete object
                line = pinch(line, *flagchars)

                flag, label = line.split(keywords['delim'])

                if flag == keywords['init']:    # start a new object
                    stack.append({})
                    stack[-1]['label'] = label
                    stack[-1]['children'] = []
                    if label == keywords['cell']:
                        cells.append(stack[-1])

                else:    # complete the most recent object
                    obj = stack.pop()
                    if len(stack) == 0:
                        hierarchy.append(obj)
                    else:
                        stack[-1]['children'].append(obj)

            # add object parameters
            else:
                var, val = line.split('=')
                stack[-1][var] = val
        except:
            continue

    return hierarchy, cells, spacing


def reorderCells(cells):
    ''' reorder cells from the bottom left corner '''

    ## bin possible values

    ys = []
    xs = []

    for cell in cells:
        x, y = cell['x'], cell['y']
        if not y in ys:
            ys.append(y)
        if not x in xs:
            xs.append(x)

    # we will use key k = a*y+x with a greater than (max(delta x)/min(delta y))

    xs.sort()
    ys.sort()

    dx_max = xs[-1]-xs[0]

    if len(ys) > 1:
        dy_min = min([ys[i]-ys[i-1] for i in xrange(1, len(ys))])
    else:
        dy_min = 1

    alpha = 1+dx_max/dy_min

    newcells = sorted(cells, key=lambda cell: alpha*cell['y']+cell['x'])

    for i in xrange(len(newcells)):
        newcells[i]['number'] = i

    return newcells


def processCells(cells):
    '''Selects out the position, size and type of each cell and creates
    a list of cell dicts {'x':x, 'y':y, 'cx':width, 'cy':height,
    'type':(NORMAL->0, OUTPUT->1 or FIXED->2)}'''

    output = []

    for i in xrange(len(cells)):
        cell = cells[i]
        obj = cell['children'][0]
        dat = {}
        dat['number'] = i
        dat['x'] = float(obj['x'])
        dat['y'] = float(obj['y'])
        #dat['dot_diam'] = float(cell['cell_options.dot_diameter'])
        dat['cx'] = float(cell['cell_options.cxCell'])
        dat['cy'] = float(cell['cell_options.cyCell'])
        dat['type'] = TYPEMAP[cell['cell_function'].rstrip()]

        if dat['type'] == TYPEMAP['QCAD_CELL_INPUT']:
            dat['pol'] = 0

        if dat['type'] == TYPEMAP['QCAD_CELL_FIXED']:  # read cell polarization
            try:
                dat['pol'] = float(cell['children'][-1]['psz'])
            except KeyError:
                print 'FIXED cell missing polarization label...'
                if KILL_ON_NO_POL:
                    sys.exit()
                else:
                    dat['pol'] = 0
            #print 'Detected FIXED type cell... pol: ' + str(dat['pol'])

        output.append(dat)

    return output


def correctNumbering(cells):

    cn, dn = 0, 0

    for cell in cells:
        if cell['type'] in [0, 1]:
            cell['number'] = cn
            cn += 1
        else:
            cell['number'] = 'D%d' % dn
            dn += 1

    return cells


def parseQCAFile(filename):
    '''Iterate through a QCAD file and construct the hierarchy of objects.
    Selects out and processes the cells for easy use. Returns list of
    cells'''

    fp = open(filename, 'r')
    data = fp.readlines()
    fp.close()

    hierarchy, cells, spacing = buildHierarchy(data)

    print 'counted ' + str(len(cells)) + ' cells'

    proc_cells = processCells(cells)

    # reorder cells

    proc_cells = reorderCells(proc_cells)

    return proc_cells, spacing

if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except:
        print 'No file detected'
        sys.exit()
    cells, spacing = parseQCAFile(fname)
