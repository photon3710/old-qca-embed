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
from auxil import coefToConn, convert_to_lim_adjacency

from parse_qca import parse_qca_file
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

    fn = 'test_circuits/feedback'#sys.argv[1]
    cells, spacing, zones, J, feedback = parse_qca_file(fn, show=False)
    J = convert_to_lim_adjacency(cells, spacing, J)

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

source = coefToConn(J, J)
print source

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
        print 'start'
        cell_map, paths = denseEmbed(source, write=True)
        print 'end'
        n_ex = sum(map(lambda x: len(x)-2, paths.values()))
        print 'Extra qubits needed: %d' % n_ex
        j += 1
    except Exception as e:
        print "ERROR: " + e.message


# run conversion method to minimize maximum model size
##models, max_model = convertToModels(paths, cell_map)
