#!/usr/bin/python

from auxil import new_getEk
import numpy as np
import pylab as plt

rot1 = np.pi/4
rot2 = np.pi/4

#rot1 = rot2 = 0

NUM_R = 100
NUM_TH = NUM_R

Rmin = 2
Rmax = 3

Thmax = np.pi/2

DX = 1e-9


def genPos(r, th, rot=0):
    ''' '''

    # Assume dx = dy = 1 and r = r/dx

    # zero centered positions
    X = .5*np.matrix([[1, 1], [1, -1], [-1, -1], [-1, 1]]).T
    R = np.matrix([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    Y = R*X

    # offset by r and th

    dx = r*np.cos(th)
    dy = r*np.sin(th)

    cell = {'qdots': []}

    for i in xrange(4):
        qd = {}
        qd['x'] = DX*(dx + Y[0, i])
        qd['y'] = DX*(dy + Y[1, i])
        cell['qdots'].append(qd)

    return cell


def main():
    ''' '''
    # generate qdot positions

    c1 = genPos(0, 0, rot1)

    R = np.linspace(Rmin, Rmax, NUM_R)
    TH = np.linspace(-Thmax, Thmax, NUM_TH)

    EKs = np.zeros([NUM_R, NUM_TH], dtype=float)

    for i_r in xrange(NUM_R):
        for i_th in xrange(NUM_TH):
            c2 = genPos(R[i_r], TH[i_th], rot2)
            EKs[i_r, i_th] = new_getEk(c1, c2)
            
    plt.figure(0)
    plt.clf()
    plt.imshow(EKs.T, extent=[Rmin, Rmax, -Thmax, Thmax])
    plt.colorbar()
    plt.show(block=False)
    
    r, th = np.meshgrid(R, TH)
    Z = np.cos(4*th)/pow(r, 5)
    
    plt.figure(1)
    plt.clf()
    plt.imshow(Z, extent=[Rmin, Rmax, -Thmax, Thmax])
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    main()

