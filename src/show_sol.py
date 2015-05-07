#!/usr/bin/python

#---------------------------------------------------------
# Name: show_sol.py
# Purpose: Test code for the dense placement algorithm
# Author:	Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

#import pylab as plt
import Tkinter as tk
from parse_qca import parseQCAFile, correctNumbering
import sys


######################################################
## GLOBALS

VERBOSE = False
USE_LABELS = False

# window parameters
WIN_W = 1400
WIN_H = 640

### drawing parameters

## numbering
CELL_LABEL_SCALE = 1.1

## spacing
CONNECT_LENGTH = .2        # percent of tile length
GRID_PAD = .1            # percent of tile length padding
PAD_PERC = .1

W_FACT = 20

## connectors
CONN_THICK = 1            # relative to qubit thickness: factor
CONN_CIRC_RAD = 1.3        # relative to CONN_THICK: factor
CONN_CIRC_COLOR = 'black'
CONN_COLOR = 'black'

## qubits

# [unused, used, used:cell, disabled]
QBIT_COLORS = ['lightgrey', 'red', 'blue', 'magenta']

QCA_COLORS = ['green', 'blue', 'red']    # [normal, output, driver]
ROUTE_MARK = '_Route.txt'
QBIT_SCALE = 1.

QCA_W = .5*WIN_W
QCA_H = WIN_H

QCA_SCALE = 1.5

QBIT_W = WIN_W-QCA_W
QBIT_H = WIN_H

######################################################
## QCA CLASS


class QCA_Frame(tk.Frame):

    def __init__(self, root, filename):

        # load qca file

        cells, spacing = parseQCAFile(filename)

        cells = correctNumbering(cells)

        # fix window size

        bounds, offset = self.findBindingBox(cells)

        region = list(bounds)
        region[0] -= offset[0]
        region[1] -= offset[1]
        region[2] -= offset[0]
        region[3] -= offset[1]

        self.region = [int(region[i]*QCA_SCALE) for i in xrange(4)]

        print 'Scroll region: ' + str(self.region)

        # initialise frame
        tk.Frame.__init__(self, root)
        self.canvas = tk.Canvas(self, width=QCA_W, height=QCA_H)
        self.canvas.configure(background='white')

        #self.xsb.grid(row=1, column=0, sticky='ew')
        #self.ysb.grid(row=0, column=1, sticky='ns')
        self.canvas.grid(row=0, column=0, sticky='nsew')
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.canvas.bind('<ButtonPress-1>', self.scroll_start)
        self.canvas.bind('<B1-Motion>', self.scroll_move)

        # draw circuit

        self.drawCircuit(cells)

        print 'QCA Frame: ' + str(self.winfo_width())

    def scroll_start(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def scroll_move(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=2)

    ## QCA drawing methods

    def findBindingBox(self, cells):
        ''' finds the corners of the smallest box containing all QCA cells.
        Returns the bounds [xmin, xmax, ymin, ymax] and the offset of the
        bottom left corner (for drawing purposes)'''

        bounds = []

        # find the bounding box

        for cell in cells:
            if len(bounds) == 0:
                bounds.append(cell['x'])
                bounds.append(cell['x']+cell['cx'])
                bounds.append(cell['y'])
                bounds.append(cell['y']+cell['cy'])
            else:
                if cell['x'] < bounds[0]:
                    bounds[0] = cell['x']
                if cell['x'] > bounds[1]:
                    bounds[1] = cell['x']
                if cell['y'] < bounds[2]:
                    bounds[2] = cell['y']
                if cell['y'] > bounds[3]:
                    bounds[3] = cell['y']

        # pad bounding box

        dx = PAD_PERC*(bounds[1]-bounds[0])
        dy = PAD_PERC*(bounds[3]-bounds[2])

        bounds[0] -= dx
        bounds[1] += dx
        bounds[2] -= dy
        bounds[3] += dy

        offset = [bounds[0], bounds[2]]

        print bounds
        print offset

        return bounds, offset

    def drawCell(self, cell, scale, offset):
        '''draw a single QCA cell (rectangle with number) in the
        global canvas'''

        canvas = self.canvas

        # shorthand parameters

        n = cell['number']
        x = (cell['x'] - offset[0])*scale    # cell x position
        y = (cell['y'] - offset[1])*scale    # cell y position
        w = cell['cx']*scale                # cell box width
        h = cell['cy']*scale                # cell box height
        #d = cell['dot_diam']*scale            # cell dot diameter

        # correct for tkinter coordinates
        y = self.region[3]-y

        # draw rectangle
        canvas.create_rectangle(x, y, x+w, y-h, fill=QCA_COLORS[cell['type']])

        # draw label: number of QCA cell
        if USE_LABELS:
            canvas.create_text(x+.5*w, y-.5*h, text=str(n), font=('Arial', 14))

    def drawCircuit(self, cells):
        '''draw all QCA cells'''

        canvas = self.canvas

        canvas.delete(tk.ALL)

        # find scaling and offset

        bounds, offset = self.findBindingBox(cells)
        #scale = min(WIN_W/(bounds[1]-bounds[0]),WIN_H/(bounds[3]-bounds[2]))
        scale = QCA_SCALE

        # draw cells

        for cell in cells:
            self.drawCell(cell, scale, offset)


######################################################
## EMBEDDING CLASS


class Qubit_Frame(tk.Frame):

    def __init__(self, root, qbits, paths, dis_qbits, M, N, L):
        ''' '''

        print dis_qbits
        self.qbits = qbits
        self.paths = paths
        self.dis_qbits = dis_qbits
        self.M = M
        self.N = N
        self.L = L

        tk.Frame.__init__(self, root)
        self.canvas = tk.Canvas(self, width=QBIT_W, height=QBIT_H)
        self.canvas.configure(background='white')

        self.canvas.grid(row=0, column=0, sticky='nsew')

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.canvas.bind('<ButtonPress-1>', self.scroll_start)
        self.canvas.bind('<B1-Motion>', self.scroll_move)

        self.drawSolution()

        print 'Qubit Frame: ' + str(self.winfo_width())

    def scroll_start(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def scroll_move(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=2)

    ## Embedding drawing methods

    def drawTile(self, origin, scale, colors):

        canvas = self.canvas

        length = scale
        width = scale//W_FACT
        delta = length//self.L

        outline = ''

        ## draw vertical qubits

        x = origin[0] + .5*delta
        y = origin[1]

        for i in xrange(0, self.L):
            color = QBIT_COLORS[colors[0][i]]
            canvas.create_rectangle(x-.5*width, y, x+.5*width, y+length,
                                    fill=color, outline=outline)
            x += delta

        ## draw horizontal qubits

        x = origin[0]
        y = origin[1] + length - .5*delta

        for i in xrange(0, self.L):
            color = QBIT_COLORS[colors[1][i]]
            canvas.create_rectangle(x, y-.5*width, x+length, y+.5*width,
                                    fill=color, outline=outline)
            y -= delta

    def drawLabels(self, scale):
        ''' '''

        canvas = self.canvas

        offset = scale*GRID_PAD
        tile_sep = scale*(1+CONNECT_LENGTH)

        qb_length = scale
        qb_delta = qb_length//self.L

        font_size = int(CELL_LABEL_SCALE*scale/10)

        for cell in self.qbits:

            if self.qbits[cell] is None:
                continue

            qbit = self.qbits[cell]

            row, col, h, l = qbit

            # top left corner of tile
            x0 = offset+col*tile_sep
            y0 = offset+tile_sep*(self.M-1)-row*tile_sep

            if h == 0:    # vertical qubit
                x = x0+l*qb_delta
                y = y0+qb_length
            else:        # horizontal qubit
                x = x0+qb_length
                y = y0+qb_length - l*qb_delta

            canvas.create_text(x, y, text=str(cell), font=('Arial', font_size))

    def drawConnectors(self, scale):
        ''' '''

        canvas = self.canvas

        length = scale*CONNECT_LENGTH
        thick = (scale//W_FACT)*CONN_THICK
        rad = thick*CONN_CIRC_RAD

        offset = scale*GRID_PAD
        tile_sep = scale*(1+CONNECT_LENGTH)

        qb_length = scale
        qb_delta = qb_length//4

        outline = ''

        for key in self.paths:
            path = self.paths[key]
            for i in xrange(len(path)-1):
                q1, q2 = path[i: i+2]
                r, c = q1[:2]
                x0 = offset + c*tile_sep
                y0 = offset + tile_sep*(self.M-1) - r*tile_sep

                if q1[:2] == q2[:2]:    # internal coupler
                    if q1[2] == 0:   # q1 is the vertical
                        x = x0+qb_delta*(q1[3]+.5)
                        y = y0 + qb_length - qb_delta*(q2[3]+.5)
                    else:
                        x = x0+qb_delta*(q2[3]+.5)
                        y = y0 + qb_length - qb_delta*(q1[3]+.5)
                    canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                       fill=CONN_CIRC_COLOR, outline=outline)

                elif q1[2]:     # horizontal external coupler
                    y = y0 + qb_length - qb_delta*(q1[3]+.5)
                    if q1[1] < q2[1]:   # coupler at right of tile
                        x = x0+qb_length
                    else:
                        x = x0-length
                    canvas.create_rectangle(x, y-.5*thick, x+length,
                                            y+.5*thick, fill=CONN_COLOR,
                                            outline=outline)

                    canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                       fill=CONN_CIRC_COLOR, outline=outline)
                    x += length
                    canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                       fill=CONN_CIRC_COLOR, outline=outline)

                else:       # vertical external coupler
                    x = x0+qb_delta*(q1[3]+0.5)
                    if q1[0] < q2[0]:   # coupler at bottom of tile
                        y = y0-length
                    else:
                        y = y0+qb_length
                    canvas.create_rectangle(x-.5*thick, y, x+.5*thick,
                                            y+length, fill=CONN_COLOR,
                                            outline=outline)

                    canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                       fill=CONN_CIRC_COLOR, outline=outline)
                    y += length
                    canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                       fill=CONN_CIRC_COLOR, outline=outline)

    def drawSolution(self):
        ''' '''

        effective_num_tiles = max(self.M, self.N)*(1+CONNECT_LENGTH)+2*GRID_PAD
        scale = min(QBIT_H, QBIT_W)/effective_num_tiles*QBIT_SCALE

        tile_sep = scale*(1+CONNECT_LENGTH)    # separation between tiles
        tile_offset = scale*GRID_PAD        # offset of tiles from edge

        x0 = tile_offset
        y0 = tile_offset + tile_sep*(self.M-1)

        # generate colors
        qbit_colors = self.setColors()

        for i in xrange(self.M):
            y = y0-i*tile_sep
            for j in xrange(self.N):
                x = x0+j*tile_sep
                colors = qbit_colors[i][j]
                self.drawTile([x, y], scale, colors)

        if USE_LABELS:
            self.drawLabels(scale)
        self.drawConnectors(scale)

    def setColors(self):

        qbit_colors = {}

        # initialise qubits as unassigned
        for i in xrange(self.M):
            qbit_colors[i] = {}
            for j in xrange(self.N):
                d = {}
                for h in [0, 1]:
                    d[h] = {}
                    for l in xrange(self.L):
                        d[h][l] = 0
                qbit_colors[i][j] = d

        for path in self.paths.values():
            for qb in path:
                r, c, h, l = qb
                qbit_colors[r][c][h][l] = 1

        for qbit in filter(None, self.qbits.values()):
            r, c, h, l = qbit
            qbit_colors[r][c][h][l] = 2

        # handle broken qbits
        for qb in self.dis_qbits:
            r, c, h, l = qb
            qbit_colors[r][c][h][l] = 3

        return qbit_colors

######################################################
### HANDLER FUNCTIONS


def pinch(st, l, r):
    return st.partition(l)[2].rpartition(r)[0]


def strToQbit(st):
    if 'None' in st:
        return None
    else:
        try:
            return tuple(map(int, pinch(st, '(', ')').split(',')))
        except:
            return None


def load_sol(fn):

    try:
        fp = open(fn, 'r')
    except IOError:
        print('Failed to open file... check filename/path')
        sys.exit()

    qb_flag = False
    pt_flag = False
    hd_flag = False
    dq_flag = False

    qbits = {}
    paths = {}
    dis_qbits = []

    for line in fp:
        if len(line) < 3:
            continue
        if '<header>' in line:
            hd_flag = True
            continue
        elif '<dis_qbits>' in line:
            dq_flag = True
            continue
        elif '<qubits>' in line:
            qb_flag = True
            continue
        elif '<paths>' in line:
            pt_flag = True
            continue
        elif '</header>' in line:
            hd_flag = False
            continue
        elif '</dis_qbits>' in line:
            dq_flag = False
            continue
        elif '</qubits>' in line:
            qb_flag = False
            continue
        elif '</paths>' in line:
            pt_flag = False
            continue

        # parse header
        if hd_flag:
            if 'M' in line:
                M = int(line.split('=')[1])
            elif 'N' in line:
                N = int(line.split('=')[1])
            elif 'L' in line:
                L = int(line.split('=')[1])

        # disabled qubits
        if dq_flag:
            q = strToQbit(line.strip())
            dis_qbits.append(q)

        # parse qubits
        if qb_flag:
            cell, qbit = line.strip().split(':')
            cell = int(cell)
            qbit = strToQbit(qbit)
            qbits[cell] = qbit

        # parse paths
        elif pt_flag:
            key, path_st = line.strip().split(':')
            key = tuple(map(int, key.split(';')))
            path = map(lambda x: strToQbit(x), path_st.split(';'))
            paths[key] = path

    return qbits, paths, dis_qbits, M, N, L


def drawSol(qbits, paths, dis_qbits, M, N, L, qca_fn=None):
    ''' '''

    root = tk.Tk()
    root.title('Embed Viewer')
    root.configure(background='white')

    frame = tk.Frame(root, width=WIN_W, height=WIN_H)
    frame.configure(background='white')

    if qca_fn:
        QCA_Frame(frame, qca_fn).pack(side='left')
    Qubit_Frame(frame, qbits, paths, dis_qbits, M, N, L).pack(side='left')

    frame.pack()

    root.mainloop()


def main(fn, qca_fn):

    qbits, paths, dis_qbits, M, N, L = load_sol(fn)

    drawSol(qbits, paths, dis_qbits, M, N, L, qca_fn=qca_fn)


if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print('No file entered')
        sys.exit()

    try:
        qca_fn = sys.argv[2]
    except:
        qca_fn = None

    main(fn, qca_fn)
