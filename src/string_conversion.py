#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Stephane
#
# Created:     29-06-2015
# Copyright:   (c) Stephane 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import numpy as np

def str_to_list(string, integer=True):
    '''inverse of str(list[])'''
    l = []
    for i in string.strip('[]').split(','):
        if i != '':
            if integer:
                l.append(int(i.strip(' ')))
            else:
                 l.append(float(i.strip(' ')))
    return l

def str_to_nparr(string, integer=True):
    '''inverse of str(ndarray)'''
    l = []
    for i in string.strip('[]').split(' '):
        if i != '':
            if integer:
                l.append(int(i))
            else:
                l.append(float(i))
    return np.asarray(l)

def str_to_list_of_tup(string, integer=True):
    '''inverse of str([(,),(,)]'''
    arr = []
    if string == '[]':
        return arr
    for i in string.strip('[]').split('),'):
        l = []
        if i.strip('( )') == '':
            l.append(tuple())
        else:
            for j in i.strip('( )').split(','):
                if integer:
                    l.append(int(j))
                else:
                    l.append(float(j))
        arr.append(tuple(l))
    return arr

def str_to_2d_list_of_tup(string, integer=True):
    '''inverse of str([[(,),(,)],[(,)...]])'''
    arr = []
    for i in string.split('],'):
        l = []
        if i.strip('[ ]') == '':
            l.append(tuple())
        else:
            for j in i.strip(',[ ]').split('),'):
                l.append(str_to_tuple(j.strip('[ ]'), integer=integer))
        arr.append(l)
    return arr

def str_to_2d_np_arr(string, integer=True):
    '''inverse function of str(np.array)'''
    arr = []
    for i in string.strip('[ ]').split('['):
        l = []
        if i.strip(' ') != '':
            for j in i.strip(',[ ]\n').split(' '):
                if j.strip('[ ]') != '':
                    if integer:
                        l.append(int(j.strip('[ ]')))
                    else:
                        l.append(float(j.strip('[ ]')))
            arr.append(np.asarray(l))
    return np.asarray(arr)

def str_to_tuple(string, integer=True):
    '''inverse operation of str(tuple)'''
    arr = string.strip('()').split(',')
    l = []
    for i in arr:
        if integer:
            l.append(int(i))
        else:
            l.append(float(i))
    return tuple(l)

def str_to_2d_tuple(string, integer=True):
    '''inverse operation of str(tuple of tuples)'''
    arr = string.strip('),').split('),')
    l = []
    for i in arr:
        if i != '':
            i = i.strip('( )').split(',')
            il = []
            for j in i:
                if j != '':
                    if integer:
                        il.append(int(j.strip('( )')))
                    else:
                        il.append(float(j.strip('( )')))
            l.append(tuple(il))

    return tuple(l)

def str_to_list_of_nparr(string, integer=True):
    '''inverse of str([nd.array])'''
    arr = string.strip('[]').split('array([')
    l = []
    for a in arr:
        a = a.strip(']), ')
        npl = []
        if a != '':
            a = a.split(',')
            for n in a:
                if n != '':
                    if integer:
                        npl.append(int(n))
                    else:
                        npl.append(float(n))
        if npl:
            l.append(np.asarray(npl))
    return l