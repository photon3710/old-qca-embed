#!/usr/bin/python

from __future__ import division
import pylab as plt
import os
import sys
from math import ceil
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq, curve_fit   # Brent method for root finding

M = 6
N = M*M*8

FILE_ROOT = '../%dgen/0/' % N
IMG_ROOT = '../../img/'

SHOW_MIN = False

FS = 18

SAVE = False

BW = 20
MAX_RANGE = N
MAX_RANGE = 600
MIN_RANGE = 10

REL = False
LOGY = False
SHOW_FIT = True

LINEAR_DENSE = True


def load_outfile(fname):
    ''' '''

    print 'reading file: %s' % fname,

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open file: %s' % fname
        return None

    data = {'lim': [], 'full': []}

    for line in fp:
        if len(line) < 3 or '#' in line:
            continue
        elif '<lim\>' in line:
            dat = data['lim']
        elif '<full\>' in line:
            dat = data['full']
        else:
            d = line.strip().split()
            n = int(d[0])
            av = float(d[1])
            mi = int(d[2])
            try:
                t = float(d[3])
                p = float(d[4])
                dat.append([n, av, mi, t, p])
            except KeyError:
                dat.append([n, av, mi])

    #print '\t(%d,%d)' % (len(data['lim']), len(data['full']))

    print max([x[3] for x in data['full']])
    return data


def plot_qubit_usage(all_data, img_dir, typ):
    ''' '''
    plt.figure('qubit-usage')
    plt.clf()
    ind = 2 if SHOW_MIN else 1

    X = {'full': [], 'lim': []}
    Y = {'full': [], 'lim': []}

    for key in ['full', 'lim']:
        c = ['g', 'r'] if key == 'full' else ['w', 'm']
        s = 'x' if key == 'full' else 'o'
        x = X[key]
        y = Y[key]
        for d in all_data[key]:
            #print '%s: \t %s' % (key, str(d))
            if d[ind] == -1:
                continue
                plt.plot(d[0], 0, c[1]+s, markersize=5, markeredgewidth=2)
            else:
                x.append(d[0])
                y.append(d[ind])
        if REL:
            y = [y[i]*1./x[i] for i in xrange(len(Y))]
        if key == 'full':
            plt.plot(x, y, c[0]+'x', markersize=5, markeredgewidth=2)
        else:
            plt.plot(x, y, c[0]+'o', markersize=5, markeredgewidth=2,
                     markeredgecolor='blue')
    x_max = max([1, max(X['full']), max(X['lim'])])
    xx = np.linspace(0, x_max, 20)
    plt.plot(xx, xx, 'k', ls='--', linewidth=2, alpha=.5)
    plt.xlim(0, x_max)
    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc='best', fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    ylab = ('Minimum' if SHOW_MIN else 'Average') + 'Qubit Usage'
    plt.ylabel(ylab, fontsize=FS)
    plt.tick_params(axis='both', labelsize=FS)
    #plt.title('Average Qubit Usage vs. Number of Cells')
    if SAVE:
        plt.savefig(img_dir+'-gen-qbits.eps', bbox_inches='tight')
    plt.show(block=False)

    return X, Y


def plot_success_prob(all_data, bw, img_dir, typ):
    ''' '''
    plt.figure('success-prob')
    plt.clf()

    X = {'full': [], 'lim': []}
    Y = {'full': [], 'lim': []}

    for key in ['full', 'lim']:

        c = ['g', 'r'] if key == 'full' else ['w', 'm']
        max_bin = max(map(lambda x: x[0], all_data[key]))
        num_bins = int(ceil(max_bin/bw))+1

        success = [0 for _ in xrange(num_bins)]
        fail = [0 for _ in xrange(num_bins)]

        for d in all_data[key]:
            if d[1] > 0:
                success[int(d[0]/bw)] += 1
            else:
                fail[int(d[0]/bw)] += 1

        for i in xrange(num_bins):
            if success[i]+fail[i] > 0:
                Y[key].append(success[i]/(success[i]+fail[i]))
            else:
                Y[key].append(0)
        X[key] = map(lambda x: bw*(.5+x), range(num_bins))

        if key is 'full':
            plt.plot(X[key], Y[key], c[0]+'x', markersize=8, markeredgewidth=3)
        else:
            plt.plot(X[key], Y[key], c[0]+'o', markersize=8, markeredgewidth=3,
                     markeredgecolor='blue')

    # generate splines
    for key in X:
        x, y = X[key], Y[key]
        f = interp1d(x, y, kind='slinear')
        xx = np.linspace(min(x), max(x), 100)
        yy = f(xx)
        if key == 'full':
            plt.plot(xx, yy, 'g--', linewidth=2)
        else:
            plt.plot(xx, yy, 'b--', linewidth=2)

    loc = 'best'
    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc=loc, fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    plt.ylabel('Success probability', fontsize=FS)
    plt.xlim([MIN_RANGE, MAX_RANGE])
    plt.tick_params(axis='both', labelsize=FS)
    if SAVE:
        plt.savefig(img_dir+'-gen-prob.eps', bbox_inches='tight')
    plt.show(block=False)

    return X, Y


def plot_runtimes(all_data, img_dir, typ):
    ''' '''
    plt.figure('run-times')
    plt.clf()

    X = {'full': [], 'lim': []}
    Y = {'full': [], 'lim': []}

    for key in ['full', 'lim']:

        x = X[key]
        y = Y[key]
        for d in all_data[key]:
            try:
                n, av, mi, t, p = d
            except:
                continue
            if t < 0:
                continue
            x.append(n)
            y.append(t)
        if LOGY:
            y = np.log(y)
        if key is 'full':
            plt.plot(x, y, 'gx', markersize=5, markeredgewidth=2)
        else:
            plt.plot(x, y, 'wo', markersize=5, markeredgewidth=2,
                     markeredgecolor=[0, 0, 1., .5])

    loc = 'best'
    plt.legend(['Full Adjacency', 'Limited Adjacency'],
               numpoints=1, loc=loc, fontsize=FS)
    plt.xlabel('Number of Cells', fontsize=FS)
    if LOGY:
        plt.ylabel('log(Run-time) (s)', fontsize=FS)
    else:
        plt.ylabel('Run-time (s)', fontsize=FS)
    plt.tick_params(axis='both', labelsize=FS)
    if SAVE:
        plt.savefig(img_dir+'-gen-runtimes.eps', bbox_inches='tight')
    plt.show(block=False)

    return X, Y


def zero_cross(x, y, offset):
    '''Find a zero crossing of (x,y-offset)'''

    # get data interpolation
    f = interp1d(x, [yy-offset for yy in y], kind='linear')
    a = x[1]
    b = x[-1]
    root, r = brentq(f, a, b, full_output=True)
    if r.converged:
        return root
    else:
        return None


def fit_runtimes(x, y, typ):
    '''Fit the runtime data to the algorithm specific model'''

    x = np.array(x)
    y = np.array(y)

    if typ == 'dense':
        if LINEAR_DENSE:
            fit_func = lambda x, a: a*x
            func = fit_func
        else:
            fit_func = lambda x, a, b: a+b*np.log10(x)
            func = lambda x, a, b: np.power(10, a)*np.power(x, b)
            y = np.log10(y)
    else:
        y = np.log10(y)
        fit_func = lambda x, a, b: a*x+b
        func = lambda x, a, b: np.power(10, a*x+b)

    popt, pcov = curve_fit(fit_func, x, y)

    # plot fit
    if SHOW_FIT:
        plt.figure('run-times')
        X = np.linspace(min(x), max(x), 100)
        Y = func(X, *popt)
        plt.plot(X, Y, 'k')

    return popt


def fit_qubits(x, y, typ):
    ''' '''

    x = np.array(x)
    y = np.array(y)
    
    print x
    print y

    fit_func = lambda x, a, b: a*x
    func = fit_func

    popt, pcov = curve_fit(fit_func, x, y)

    if SHOW_FIT:
        plt.figure('qubit-usage')
        X = np.linspace(min(x), max(x), 100)
        Y = func(X, *popt)
        plt.plot(X, Y, 'k')

    return popt


def analyse(data, typ):
    ''' '''

    metrics = {'full': {}, 'lim': {}}

    # largest embedded circuit
    for key in metrics:
        metrics[key]['x_max'] = max(data['qb_x'][key])
        metrics[key]['y_max'] = max(data['qb_y'][key])

    # 50% success prob
    for key in metrics:
        x, y = data['sp_x'][key], data['sp_y'][key]
        root = zero_cross(x, y, 0.5)
        metrics[key]['50p'] = root

    # qubit usage fit
    for key in metrics:
        popt = fit_qubits(data['qb_x'][key], data['qb_y'][key], typ)
        metrics[key]['qb_popt'] = popt

    # runtime fit
    for key in metrics:
        popt = fit_runtimes(data['rt_x'][key], data['rt_y'][key], typ)
        metrics[key]['rt_popt'] = popt

    # echo metrics
    for key in metrics:
        print('\n**** {0} ****'.format(key.upper()))
        print('Largest embedded circuit: {0}'.format(metrics[key]['x_max']))
        print('Most qubits used: {0}'.format(metrics[key]['y_max']))
        print('50 \% threshold: {0}'.format(metrics[key]['50p']))
        print('Qubit Fit Model: ax -> {0}'.format(
            *metrics[key]['qb_popt']))
        if typ == 'dense':
            if LINEAR_DENSE:
                print('Runtime Fit model: ax+b -> {0}'.format(
                    *metrics[key]['rt_popt']))
            else:
                a, b = metrics[key]['rt_popt']
                a = np.power(10, a)
                print('Runtime Fit model: ax^b -> {0}, {1}'.format(a, b))
        else:
            print('Runtime Fit model: 10^(ax+b) -> {0}, {1}'.format(
                *metrics[key]['rt_popt']))

    return metrics


def main(file_dir, img_dir, bw, typ='dense'):
    ''' '''

    fnames = os.listdir(file_dir)
    fnames = map(lambda s: file_dir+s, fnames)

    all_data = {'lim': [], 'full': []}

    for fname in fnames:
        data = load_outfile(fname)
        if not data is None:
            all_data['lim'] += data['lim']
            all_data['full'] += data['full']

    data = {}

    # plot qubit usage
    x, y = plot_qubit_usage(all_data, img_dir, typ)
    data['qb_x'] = x
    data['qb_y'] = y

    for key in all_data:
        print('{0}: {1} trials'.format(key, len(all_data[key])))

    # plot success probability
    x, y = plot_success_prob(all_data, bw, img_dir, typ)
    data['sp_x'] = x
    data['sp_y'] = y

    # plot run-times
    x, y = plot_runtimes(all_data, img_dir, typ)
    data['rt_x'] = x
    data['rt_y'] = y

    # obtain metrics
    analyse(data, typ)

    # delay pyplot
    plt.show(block=True)


if __name__ == '__main__':

    try:
        typ = sys.argv[1].strip().lower()
        assert typ in ['dense', 'heur'], 'Invalid algorithm type'
    except:
        typ = 'dense'

    try:
        bw = int(sys.argv[2])
    except:
        bw = BW

    if typ == 'dense':
        file_dir = FILE_ROOT + 'dense/'
        img_dir = IMG_ROOT + 'dense'
    else:
        file_dir = FILE_ROOT + 'heur/'
        img_dir = IMG_ROOT + 'heur'

    main(file_dir, img_dir, bw, typ=typ)
