#!/usr/bin/python

from convert import convertToModels, writeSol
#from math import ceil
import os
import re
import sys


def strToQbit(s):
    ''' convert a formatted string to a qubit tuple '''
    r, c, h, l = map(int, s.strip()[1:-1].split(','))
    return (r, c, h, l)


def readSub(data):
    ''' '''

    read_flag = False
    rx = re.compile('<.*>')
    dat = []
    while data:
        line = data.pop(0)
        if len(line) < 3 or line.strip()[0] == '#':
            continue
        if not read_flag:
            m = rx.search(line)
            if not m is None:
                name = m.group(0)[1:-1]
                read_flag = True
                continue
        if read_flag:
            if '</%s>' % name in line:
                break
            dat.append(line.strip('\n'))

    return dat


def loadSol(fname):
    '''Load a Dense Placement solutions file'''

    try:
        fp = open(fname, 'r')
    except:
        print 'Failed to open filename: %s' % fname
        return None

    data = fp.readlines()
    fp.close()

    header_dat = readSub(data)
    dis_qbit_dat = readSub(data)
    qbit_dat = readSub(data)
    path_dat = readSub(data)

    # parse header data
    M = int(header_dat[0].split()[-1])
    N = int(header_dat[1].split()[-1])
    L = int(header_dat[2].split()[-1])

    # parse dis_qbit data
    dis_qbits = dis_qbit_dat
    dis_qbits = []

    # parse qbit data
    qbits = {}
    for line in qbit_dat:
        cell, qbit = line.split(':')
        cell = int(cell)
        qbit = strToQbit(qbit)
        qbits[cell] = qbit

    # parse paths data
    paths = {}
    for line in path_dat:
        key, path = line.split(':')
        key = tuple(map(int, key.strip().split(';')))
        path = map(strToQbit, path.split(';'))
        paths[key] = path

    sol_dict = {}
    sol_dict['M'] = M
    sol_dict['N'] = N
    sol_dict['L'] = L
    sol_dict['dis_qbits'] = dis_qbits
    sol_dict['qbits'] = qbits
    sol_dict['paths'] = paths

    return sol_dict


def main(fname):
    ''' '''

    sol_dict = loadSol(fname)
    models, max_models = convertToModels(sol_dict['paths'], sol_dict['qbits'])
    return models, max_models


def runDirectory(dname):
    ''' '''

    regex = re.compile('^sol[0-9]+$')

    if not dname[-1] == '/':
        dname += '/'

    try:
        fnames = filter(regex.match, os.listdir(dname))
    except:
        print 'Invalid directory path: %s' % dname
        return None

    # create write directory
    write_dir = dname+'conv_dir/'
    if not os.path.isdir(write_dir):
        os.mkdir(write_dir)

    for fname in fnames:
        models, max_model = main(dname + fname)
        if not models is None:
            # write model to file
            fn = write_dir + fname
            if os.path.isfile(fn) or os.path.isdir(fn):
                print 'File already exists'
                continue
            if not writeSol(models, max_model, fn):
                print 'Failed to write to file: %s' % fn


if __name__ == '__main__':

    try:
        dname = sys.argv[1]
        runDirectory(dname)
    except Exception as e:
        print e.message
        print 'No file entered'
