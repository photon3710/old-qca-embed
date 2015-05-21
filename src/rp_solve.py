#!/usr/bin/python

from parse_qca import parseQCAFile
from auxil import generateAdjDict, adjToCoef, convertToNearestNeighbour,\
    stateToPol
from solve import solve, solveSparse
from pymetis import part_graph
from math import ceil
import scipy.sparse as sp
from time import time

import numpy as np
import networkx as nx
import sys

from pprint import pprint

## SOLVER PARAMETERS

# threshold values
N_THRESH = 8            # maximum partition size for exact solver
MEM_THRESH = 1e6        # maximum mode count product
STATE_THRESH = 0.05     # required amplitude for state contribution

# flags and counts
N_PARTS = 2             # maximum number of partitions at each iteration
                        # best results seem to be for 2 partitions
USE_NN = False          # Convert to nearest neighbour

# resolution
E_RES = 1e-1            # resolution for energy binning relative to max J


def comp_on_comp(h, J, gam, modes):
    '''Compute the on-comp parameters for the given modes'''

    h = np.matrix(h)
    J = np.matrix(np.triu(J))

    N = h.shape[1]

    sz = h*modes.T
    sz2 = np.diag(modes*J*modes.T)

    if gam == 0.:
        g = None
    else:
        mds = np.asarray(modes)
        diff = np.abs(mds.reshape([-1, 1, N])-mds.reshape([1, -1, N]))
        g = gam*(np.sum(diff, axis=2) == 2)

    return [sz, sz2, g]


def comp_off_comp(J, mds1, mds2):
    '''Compute the comp-comp mode coupling parameters'''

    JJ = mds1*J*mds2.T

    return JJ


def gen_comp_Ham(h, J, gam, ind, modes):
    '''Formulate the system Hamiltonian for a given list of component indices
    and modes.

    inputs: h (list)    :   on-site cell parameters
            J (mat)     :   cell-cell coupling parameters, not upper diag
            gam (float) :   on-site cell tunneling parameter
            ind (list)  :   list of lists of component indices. Lists must be
                            non intersecting and have a union containing all
                            cell indices.
            modes (list):   list of mode matrices. Each row in the mode matrix
                            has elements of +- 1 giving the mode spins.'''

    # check validity of input arguments
    N = len(h)
    all_ind = reduce(lambda x, y: x+y, ind)
    assert len(all_ind) == len(set(all_ind)), 'Repeated index...'
    assert sorted(all_ind) == range(N), 'Incomplete index set...'

    # sort component indices
    ind = map(sorted, ind)

    # create the h, J, gam parameters for each component
    h = np.array(h)
    gam = np.array(gam)
    h_m = [h[ind_m] for ind_m in ind]
    J_m = {}
    for i in xrange(N):
        for j in xrange(i, N):
            J_m[(i, j)] = J[ind[i], :][:, ind[j]]

    # on-comp parameters
    on_comp = []
    for i in xrange(N):
        on_comp.append(comp_on_comp(h_m[i], J_m[(i, i)], gam, modes[i]))

    # coupling parameters
    off_comp = {}
    for i in xrange(N):
        for j in xrange(i+1, N):
            off_comp[(i, j)] = comp_off_comp(J_m[(i, j)], modes[i], modes[j])


def partition(h, J, nparts):
    '''Split graph into nparts partitions: returns the indices, h and J
    parameters for each partition as well as coupling parameters between each
    partition'''

    # networkx graph
    G = nx.Graph(J)

    # connectivity list
    conn = G.adjacency_list()

    # partition in nparts
    ncuts, labels = part_graph(nparts, adjacency=conn)

    # indices of each partition
    parts = []
    for v in xrange(nparts):
        parts.append([i for i, x in enumerate(labels) if x == v])

    # make sure indices are sorted
    parts = map(sorted, parts)

    # make h, J for each partition
    h = np.array(h)
    h_p = [h[parts[i]] for i in xrange(nparts)]
    J_p = [J[parts[i], :][:, parts[i]] for i in xrange(nparts)]

    # coupling matrices
    C_p = {}
    for i in xrange(nparts):
        for j in xrange(i+1, nparts):
            C_p[(i, j)] = J[parts[i], :][:, parts[j]]

    return parts, h_p, J_p, C_p


def get_prod_states(state, comp=False):
    '''Get spin representation of the product states which contribute to the
    given state'''

    # number of spins
    N = int(np.log2(state.shape[0]))

    # isolate contributing product states
    inds = (np.abs(state) > STATE_THRESH).nonzero()[0]

    # sort by contribution magnitude
    inds = sorted(inds, key=lambda x: np.abs(state)[x], reverse=True)

    if comp:
        prod_states = inds
    else:
        prod_states = []
        for ind in inds:
            bstr = format(ind, '#0%db' % (N+2))[2::]
            ps = tuple(np.array(map(int, bstr))*2-1)
            prod_states.append(ps)

    return prod_states


def proc_solve(spec, sols, e_res):
    '''Process output from 'solve' function to match rp_solve'''

    states = sols['eigsh']['vecs']
    n_states = states.shape[1]

    states = [states[:, i] for i in xrange(n_states)]

    # associate energies with product states

    Egaps = list(spec-spec[0])
    prod_states = [get_prod_states(state) for state in states]

    # bin energies and correct
    Ebin = []
    PSbin = []

    # calculate energy resolution
    while Egaps:
        E = Egaps[0]
        try:
            i = next(i for i, x in enumerate(Egaps) if x > E+e_res)
        except:
            i = None
        Ebin.append(E)
        ps = set(reduce(lambda x, y: x+y, prod_states[:i]))
        PSbin.append(ps)
        if i is None:
            Egaps = []
            prod_states = []
        else:
            Egaps = Egaps[i:]
            prod_states = prod_states[i:]

    # redefine Egaps and prod_states
    Egaps = Ebin
    prod_states = PSbin

    return Egaps, states, prod_states


def select_modes(Es, PS):
    '''Based on the energy gaps for the solution of each partition, select
    which product state modes are included. Input modes should be formated
    as vectors of +- 1.'''

#    print Es
#    print PS
    N = len(Es)     # number of partitions

    i_p = [1]*N     # number of included indices for each partition
    m_p = [0]*N     # number of included modes for each partition

    prod = 1.

    # force inclusion of ground state modes
    for p in xrange(N):
        m_p[p] += len(PS[p][0])
        prod *= m_p[p]

    # check fill condition
    if prod > MEM_THRESH:
        print 'Not enough memory capacity to facilitate ground state inclusion'
        return None

    # determine number of new product states for each partition
    K = [[] for i in xrange(N)]
    for i in xrange(N):
        ps = set(PS[i][0])
        for j in xrange(1, len(PS[i])):
            temp = len(ps)
            ps.update(PS[i][j])
            K[i].append((Es[i][j], i, len(ps)-temp))

    # order K values for prioritising mode inclusion
    order = sorted(reduce(lambda x, y: x+y, K))

    # main inclusion loop
    check = [False]*N   # set True if term tejected
    while prod < MEM_THRESH and not all(check) and order:
        E, p, k = order.pop(0)  # candidate
        # shouldn't add higher energy term if lower energy rejected
        if check[p]:
            continue
        # check if adding new components pushes prod over threshold
        temp = prod
        prod += k*prod/m_p[p]
        if prod >= MEM_THRESH:
            prod = temp
            check[p] = True
        else:
            i_p[p] += 1
            m_p[p] += k

    # format the modes for output: list of mode matrices

    modes = []

    for p in xrange(N):
        mds = []
        for i in xrange(i_p[p]):
            mds += PS[p][i]
        mds = sorted(set(mds))
        modes.append(np.matrix(mds))

    print '%d modes selected...' % prod
    print 'mode dist: %s' % str(m_p)

    return modes


def gen_comp_diag(sz, sz2, J_m):
    '''Generate the diagonal terms in the component formalism Hamiltonian'''

    N = len(sz)     # number of components

    print 'Constructing diagonal component formalism parameters'
    # generate size counters
    M = [sz[p].shape[1] for p in xrange(N)]  # number of modes per partition
    C = np.cumprod(M)   # cumulative product of mode sizes
    C = np.insert(C, 0, 1)  # useful for avoiding special case

    # on comp energy terms
    print '\tComputing on-comp energy terms...'
    on_site = np.zeros([1, C[-1]], dtype=float)
    for p in xrange(N):
        a = sz[p] + sz2[p]
        temp = np.tile(np.repeat(a, C[p]), C[-1]/C[p+1])
        on_site += temp

    # external comp-comp mode coupling
    print '\tComputing comp-comp mode coupling terms...'
    ext_coup = np.zeros([1, C[-1]], dtype=float)
    for p1 in xrange(N):
        for p2 in xrange(p1+1, N):
            # handle p1 behaviour
            j = J_m[(p1, p2)].T
            if not j.any():
                continue
            #j = np.arange(j.size).reshape(j.shape)
            a = np.tile(np.repeat(j, C[p1], axis=1), C[p2]/C[p1+1])
            # handle p2 behaviour
            b = np.tile(a.flatten(), C[-1]/C[p2+1])
            ext_coup += b
    print '\t...done'
    return on_site + ext_coup


def gen_comp_tunn(G):
    '''Generate the off-diagonal tunneling terms in the component formalism
    Hamiltonian'''

    print 'Constructing off-diagonal component formalism parameters'

    N = len(G)  # number of partitions
    Nm = [x.shape[0] for x in G]    # number of modes per partition
    Cm = np.insert(np.cumprod(Nm), 0, 1)     # size of each partition sub-block

    # for each mode update the 'diagonal' submatrix
    mat = sp.coo_matrix((1, 1), dtype=float)
    for p in xrange(N):
        # construct the sub-block container
        A = []
        for m in xrange(Nm[p]):
            a = [None]*m + [mat]
            for n in xrange(m+1, Nm[p]):
                if G[p][n, m] == 0:
                    a.append(None)
                else:
                    dat = np.repeat(G[p][n, m], Cm[p])
                    a.append(sp.diags([dat], [0]))
            A.append(a)
        mat = sp.bmat(A)

    print '\t...done'
    return mat.T


def general_decomp(n, c):
    '''Decompose a number into a general basis cumprod c'''

    rep = []
    for i in xrange(len(c)):
        t, n = divmod(n, c[i])
        rep.append(t)

    return rep


def correct_prod_state(pstates, modes, inds):
    '''Correct product state to account for mode space representation'''

    t = time()
    print 'Correcting product states...',

    # prep work

    Nps = len(pstates)  # number of product state lists
    N = len(modes)      # number of partitions

    inds = reduce(lambda x, y: x+y, inds)

    nmodes = [mds.shape[0] for mds in modes]    # number of modes per part
    C = np.cumprod([1]+nmodes[:-1])[::-1]     # cumprod of mode counts

    ps_modes = []
    # for each product state list
    for i in xrange(Nps):
        ps_list = pstates[i]
        ps_mds = []
        for ps in ps_list:
            rep = general_decomp(ps, C)[::-1]
            ps_m = [modes[j][rep[j]].tolist()[0] for j in xrange(N)]
            ps_m = reduce(lambda x, y: x+y, ps_m)
            # reorder using indices
            ps_mds.append(tuple([x[1] for x in sorted(zip(inds, ps_m))]))
        ps_modes.append(ps_mds)

    print 'done'
    print 'correct prod time: %.5f s' % (time()-t)
    return ps_modes


def proc_comp_solve(sols, modes, inds, e_res):
    '''Process the output of the sparse solver for the component formalism
    Hamiltonian'''

    t = time()
    print '\nProcessing component solver...'

    spec = sols['eigsh']['vals']

    states = sols['eigsh']['vecs']
    n_states = states.shape[1]

    print '\n\n%d states solved...' % n_states
    states = [states[:, i] for i in xrange(n_states)]

    Egaps = list(spec-spec[0])
    prod_states = [get_prod_states(state, comp=True)
                   for state in states]

    # bin energies and correct
    Ebin = []
    PSbin = []
    while Egaps:
        E = Egaps[0]
        try:
            i = next(i for i, x in enumerate(Egaps) if x > E+e_res)
        except:
            i = None
        Ebin.append(E)
        ps = set(reduce(lambda x, y: x+y, prod_states[:i]))
        PSbin.append(ps)
        if i is None:
            Egaps = []
            prod_states = []
        else:
            Egaps = Egaps[i:]
            prod_states = prod_states[i:]

    # redefine Egaps and prod_states
    Egaps = Ebin
    prod_states = PSbin

    # correct states to account for mode space
    prod_states = correct_prod_state(prod_states, modes, inds)

    print '\t...done'
    print 'Proc comp time: %.5f s' % (time()-t)

    return Egaps, states, prod_states


def solve_comp(h_p, J_p, C_p, gam, modes, inds, e_res):
    '''Solve component formalism'''

    print '\nRunning component solver...'
    t = time()

    Egaps = []
    states = []
    prod_states = []

    N = len(h_p)    # number of partitions

    # get on component parameter
    outputs = [comp_on_comp(h_p[p], J_p[p], gam, modes[p]) for p in xrange(N)]
    sz, sz2, g = map(list, zip(*outputs))

    # get component mode coupling parameters
    J_m = {}
    for i in xrange(N):
        for j in xrange(i+1, N):
            J_m[(i, j)] = comp_off_comp(C_p[(i, j)], modes[i], modes[j])

    # construct diagonal elements
    diag = gen_comp_diag(sz, sz2, J_m)
    diag = sp.diags([diag[0]], [0])

    # construct tunneling elements as sparse matrix
    if not g[0] is None:
        off_diag = gen_comp_tunn(g)
        H = diag + off_diag
    else:
        H = diag

    # run sparse matrix solver
    print 'H matrix size: %s' % str(H.shape)
    print 'Running sparse solver...',
    t1 = time()
    sols = solveSparse(H)
    print 'done'
    print 'sparse solver time: %.5f s' % (time()-t1)

    # process output
    Egaps, states, prod_states = proc_comp_solve(sols, modes, inds, e_res)

    #pprint(prod_states[0:2])
    #print Egaps
    print '...done'
    print 'Component solver time: %.5f s' % (time()-t)
    return Egaps, states, prod_states


def rp_solve(h, J, gam):
    '''Solve ising spin-glass configuration using recursive partitioning
    with low-energy spectrum mode composition'''

    print 'Detected problem size: %d...' % len(h)

    e_res = np.max(np.abs(J))*E_RES
    print e_res

    if len(h) <= N_THRESH:
        print 'Running exact solver...'
        gs, es, spectrum, sols = solve(h, J, gamma=gam, full_output=True)
        Egaps, states, prod_states = proc_solve(spectrum, sols, e_res)
    else:
        print 'Running recursive partition...'
        # ensure J not triu
        J = np.triu(J)+np.triu(J, 1).T

        # partition into some number of subgraphs
        nparts = int(min(ceil(len(h)*1./N_THRESH), N_PARTS))
        parts, h_p, J_p, C_p = partition(h, J, nparts)

        # solve each partition
        Es = []
        PS = []
        for i in xrange(nparts):
            Egaps, states, prod_states = rp_solve(h_p[i], J_p[i], gam)
            Es.append(Egaps)
            PS.append(prod_states)

        # select modes to include
        modes = select_modes(Es, PS)

        # solve component system
        Egaps, states, prod_states = solve_comp(h_p, J_p, C_p, gam,
                                                modes, parts, e_res)

    return Egaps, states, prod_states


def echo_ps(ps):
    '''Nice output format for product states'''
    s = ''.join(['+' if p < 0 else '-' for p in ps])
    return s

if __name__ == '__main__':

    try:
        fn = sys.argv[1]
    except:
        print('No input file...')
        sys.exit()

    cells, spacing = parseQCAFile(fn)
    adj, drivers = generateAdjDict(cells, spacing)
    if USE_NN:
        adj = convertToNearestNeighbour(adj)
    h, J = adjToCoef(adj)

    gam = np.max(np.abs(J))*.1
    #gam = 0
    Egaps, states, prod_states = rp_solve(h, J, gam)

    #print np.round(stateToPol(states[0]),1)
    #print np.round(stateToPol(states[1]),1)
    print '\n'.join(map(echo_ps, prod_states[0]))
    print '\n\n'
    print '\n'.join(map(echo_ps, prod_states[1]))
    print '\n\n'
    #pprint(prod_states[1])
    print Egaps
