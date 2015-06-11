from new_parse_qca import parse_qca_file
from scipy.optimize import curve_fit

import numpy as np
import sys
import pylab as plt

## Physical Parameters
eps0 = 8.85412e-12  # permittivity of free space
epsr = 12.          # relative permittivity
q0 = 1.602e-19      # elementary charge


def getEk(c1, c2, scale=1.):
    '''Compute the kink energy using the qdot positions (in nm) of two cells.
    If the cell displacement is greater than DR return False'''

    qdots_1 = c1['qdots']
    qdots_2 = c2['qdots']

    # compute displacements

    x1 = [qd['x'] for qd in qdots_1]
    y1 = [qd['y'] for qd in qdots_1]

    x2 = [qd['x'] for qd in qdots_2]
    y2 = [qd['y'] for qd in qdots_2]

    X1 = np.array([x1, y1]).T.reshape([4, 1, 2])
    X2 = np.array([x2, y2]).T.reshape([1, 4, 2])

    R = np.sqrt(np.sum(pow(X1-X2, 2), axis=2))*scale

    if np.min(R) == 0:
        print 'qdot overlap detected'
        sys.exit()

    # QCADesigner orders qdots either CW or CCW so same and diff configurations
    # are always alternating indices.

    Q = np.array([1, -1, 1, -1])    # template for charge arrangement

    Q = np.outer(Q, Q)

    Ek = -1e9*q0*np.sum((Q/R))/(8*np.pi*eps0*epsr)

    return Ek


scale = np.linspace(.01, 10, 1000)

fn = '../dat/test_circ/t2'

cells, spacing, zones, J = parse_qca_file(fn, one_zone=True)

try:
    c1, c2 = cells[:2]
except:
    print('Insufficients cells in QCA circuit... need at least 2')
    sys.exit()

Eks = []
for s in scale:
    Eks.append(getEk(c1, c2, scale=s))
    
#plt.plot(np.log(scale), np.log(Eks))
#plt.xlabel('scale')
#plt.ylabel('Ek')
#plt.plot()

print spacing
print 'Base factor: {0}'.format(getEk(c1, c2, scale=1./spacing)*epsr)

