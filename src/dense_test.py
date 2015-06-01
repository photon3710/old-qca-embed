#---------------------------------------------------------
# Name: dense_test.py
# Purpose: Test code for the dense placement algorithm
# Author:	Jacob Retallick
# Created: 02.08.2014
# Last Modified: 06.05.2015
#---------------------------------------------------------

from dense_placement.embed import denseEmbed, setChimeraSize, \
    getCouplerFlags, setQbitAdj
from dense_placement.convert import convertToModels

from build import createInv, createMAJ, createWire
from auxil import generateAdjDict, adjToCoef, coefToConn, \
    convertToNearestNeighbour

from parse_qca import parseQCAFile
from random import random
import sys

M, N, L = 8, 8, 4

ndis = 0

dis_coup = []
dis_qbits = list(set(map(int, [random()*M*N*L*2 for _ in xrange(ndis)])))

typ = 'inv'
NEAREST = True
NUM_SUCCESS = 5

# initialise chimera parameters
setChimeraSize(M, N, L)

# get coupler flags
CF = getCouplerFlags(dis_coup, dis_qbits)

# set up chimera qbit adjacency dict
setQbitAdj(CF)

# read from file or use default generator to build QCA circuit
try:
    fname = sys.argv[1]
    cells, spacing = parseQCAFile(fname)
except:
    print 'No input file detected... running default generator'
    if typ == 'wire':
        NN = 8
        cells, spacing = createWire(N=[NN], C=[-1], P=[-1, None])
    elif typ == 'maj':
        cells, spacing = createMAJ([2, 2, 2], [1, 1, 1])
    elif typ == 'inv':
        cells, spacing = createInv([1, 2, 6, 4], [1])
    else:
        print 'Unrecognized circuit type'
        sys.exit()


adj, drivers = generateAdjDict(cells, spacing)
if NEAREST:
    adj = convertToNearestNeighbour(adj, drivers)
h, J = adjToCoef(adj)
source = coefToConn(h, J)

if False:
    for i in source:
        print '%d: %s' % (i, source[i])


print '\n\n'
i = 1   # iteration count
j = 0   # number of successful embeddings
while j < NUM_SUCCESS:
    print '\n\nRUNNING TRIAL %d :: %d\n\n' % (i, j)
    i += 1
    try:
        cell_map, paths = denseEmbed(source, write=False)
        n_ex = sum(map(lambda x: len(x)-2, paths.values()))
        print 'Extra qubits needed: %d' % n_ex
        j += 1
    except Exception as e:
        #print e.message
        pass


# run conversion method to minimize maximum model size
models, max_model = convertToModels(paths, cell_map)
