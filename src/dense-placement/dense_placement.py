#---------------------------------------------------------
# Name:		dense_placement.py
# Purpose:	
# Author:	Jacob Retallick
# Created:	2014-08-06
# Last Modified: 2014-09-14
#---------------------------------------------------------

from __future__ import division
from random import randint, random, shuffle
from bisect import bisect
import numpy as np
import math,sys
from copy import deepcopy as cp

import routing as Routing


## Globals

# constants

M=8	# Number of rows
N=8	# Number of columns
L=4	# Number of qubits per half tile

# WORKING VARIABLES
	
_numAdj={} 		# number of unplaced adjacent cells, source indexed
_source={}		# adjacency list for cell connectivity

_doNow=[]		# ordered list of cells to be placed 
_doNext=set() 	# set of cells to be placed after doNow

_qubits={}		# source keyed dict of each cell's qubit
_cells={}		# 4-tup keyed dict of each qbit's cell
_qbit_paths={}	# 4-tup-keyed dict of paths containing each qbit
_qbitAdj={}		# 4-tup keyed dict of adjacent qubit lists

_cell_flags={}	# source keyed flag dicts for each cell
_qbit_flags={}	# 4-tuple keyed flag dicts for each qubit
_reserved={}	# source keyed dict of sets of taken adjacent qubits

_paths=[]		# list of all paths in the embedding

_vacancy[0,0,0,0]	# number of free columns/rows of tiles [L,R,D,U]
_tile_occ={}		# number of qbits in each tile row/col

## Handles/Knobs

FIRST_PROB_POW=3		# power for firstCell probabilistic method
FIRST_QBIT_SIG=1.		# gaussian dist. tile select standard dev.
FIRST_QBIT_ATTEMPTS=5	# number of attempts to find starting qbit

WRITE=True
ROUTE_PATH='routing' if WRITE else ''

	
	
	
	
	
# unchecked
def indexToLinear(tup,index0=False):
	'''convert a 4-tuple index to a linear qubit index. Tuple format
	tup=(row,col,horiz?,index) with row,col,index starting from 
	bottom left tile/qubit'''
	
	global M,N,L
	
	qpr=2*N*L 	# qbits per row
	qpt=2*L		# qbits per tile
	
	return (0 if index0 else 1)+qpr*tup[0]+qpt*tup[1]+L*tup[2]+tup[3]
	
# unchecked
def indexToTuple(index,index0=False):
	'''converts a linear qubit index to a 4-tuple format'''
	
	global M,N,L
	qpr=2*N*L 	# qbits per row
	qpt=2*L		# qbits per tile
	
	if not index0:
		index-=1
			
	row,rem=divmod(index,qpr)
	col,rem=divmod(rem,qpt)
	horiz,ind=divmod(rem,L)
	
	return (row,col,horiz,ind)

	
	
	
	
	
# unchecked
def setChimeraSize(m,n,l):
	'''
	updates the Chimera graph size.
	
	inputs: m (int)	: number of tile rows
			n (int)	: number of tile columns
			l (int)	: number of horizontal or vertical qubits per tile
			
	outputs: none
	'''
	
	global M,N,L
	M,N,L=m,n,l

# unchecked
def getCouplerFlags(dis_coup=[], dis_qbits=[]):
	'''Returns a dictionary of coupler flags with 4-tuple pair
	indexing. Inputs should have linear indices starting with 1.'''
	
	global M,N,L
	
	coupler_flags={}
	
	# generate all couplers for M,N,L: flag default True
	# keys of type (q1,q2) where q1<q2
	for row in xrange(M):
		for col in xrange(N):	# for each tile
			# vertical qubits
			for v in xrange(L):
				q1=(row,col,0,v)
				# internal couplers
				for h in xrange(L):
					q2=(row,col,1,h)
					coupler_flags[(q1,q2)]=True
				# external couplers
				if row<(M-1):
					q2=(row+1,col,0,v)
					coupler_flags[(q1,q2)]=True
			# horizontal qubits, vertical handles internal couplers
			for h in xrange(L):
				q1=(row,col,1,h)
				# external couplers
				if col<(N-1):
					q2=(row,col+1,1,h)
					coupler_flags[(q1,q2)]=True
					
	# set disabled couplers as False
	
	for i1,i2 in dis_coup:
		if i1==i2:
			print 'Self directed coupler detected ...%d' %i1
			sys.exit()
		if i2<i1:
			i1,i2=i2,i1
		q1,q2 = indexToTuple(i1),indexToTuple(i2)
		## catch check
		# print i1, q1
		# print i2, q2
		try:
			coupler_flags[(q1,q2)]=False
		except KeyError:
			print 'Invalid coupler: %s -> %s, the grid size is likely incorrect' %(i1,i2)
			sys.exit()
	
	# deactive couplers which connect to disabled qubits
	for qbit in dis_qbits:
		q=indexToTuple(qbit)
		if q[2]==0: # vertical qubit
			for h in xrange(L):
				q2=(q[0],q[1],1,h)
				coupler_flags[(q,q2)]=False
			if q[0]>0:
				q1=(q[0]-1,q[1],q[2],q[3])
				coupler_flags[(q1,q)]=False
			if q[0]<(M-1):
				q2=(q[0]+1,q[1],q[2],q[3])
				coupler_flags[(q,q2)]=False
		else: # horizontal qubit
			for v in xrange(L):
				q1=(q[0],q[1],0,v)
				coupler_flags[(q1,q)]=False
			if q[1]>0:
				q1=(q[0],q[1]-1,1,q[3])
				coupler_flags[(q1,q)]=False
			if q[1]<(N-1):
				q2=(q[0],q[1]+1,1,q[3])
				coupler_flags[(q,q2)]=False
				
	return coupler_flags



# unchecked
def initialize(coupler_flags):
	''' '''
	global _qbitAdj
	pass

# unchecked
def reset():
	''' '''
	pass
	
# unchecked
def setQbitAdj(coupler_flags, qbitAdj):
	'''reset qubit adjacency by reference'''
	
	# reset qubit adjacenecy
	for n in xrange(N):
		for m in xrange(M):
			for h in xrange(2):
				for l in xrange(L):
					qbitAdj[(n,m,h,l)]=[]
	
	# generate qubit adjacenecy				
	for coupler in coupler_flags:
		if coupler_flags[coupler]:
			q1,q2=coupler
			qbitAdj[q1].append(q2)
			qbitAdj[q2].append(q1)
			
	# sort each keyed list
	for key in qbitAdj:
		qbitAdj[key].sort()
	
# unchecked
def setFlags():
	'''Initialise *_flags and reserved dicts'''
	global _numAdj, _qbitAdj, _cell_flags, _qbit_flags, _reserved
	
	# initialise cell_flags
	
	for key in _numAdj:
		d={}
		d['placed']=True
		_cell_flags[key]=d
		_reserved[key]=set()
		
	# initialise qbit_flags
	for key in _qbitAdj:
		d={}
		d['taken']=False
		d['reserved']=False
		d['assigned']=False
		_qbit_flags[key]=d
	
	
	
	
	
	
# unchecked
def setVancancy():
	''' '''
	global _tile_occ, _vacancy
	
	# compute left/right vacancy
	occupied=map(lambda x:x>0,_tile_occ['c'])
	left=occupied.index(True)
	right=occupied[::-1].index(True)
	
	# compute bottom/top vacancy
	occupied=map(lambda x:x>0,_tile_occ['r'])
	bot=occupied.index(True)
	top=occupied[::-1].index(True)
	
	_vancancy=[left,right,bot,top]
		
# unchecked
# affects tile_occ
def reserveQubits(cells):
	'''for each cell in cells, check if adjacent qubits should be
	reserved. Reserve if appropriate'''
	global _reserved, _qbit_flags, _qubits, _qbitAdj, _numAdj
	
	for cell in cells:
		
		# get qbit assigned to cell
		try:
			cell_qb = _qubits[cell]
		except KeyError:
			print 'Cell <%s> has no assigned qubit...' %str(cell)
			continue
		
		# get list of all adjacent unplaced qubits
		qbs=[q for q in _qbitAdj[cell_qb] if not _qbit_flags[q]['taken']]
		
		# wipe old reservations
		for qb in _reserved[cell]:
			_qbit_flags[qb]['reserved']=False
		_reserved[cell].clear()
		
		if _numAdj[cell]==len(qbs):
			# reserve all adjacent cells
			for qb in qbs:
				_qbit_flags[qb]['reserved']=True
				_reserved[cell].add(qb)
			
		elif _numAdj[cell]>len(qbs):
			# insufficient adjacent qubits... reserve error
			print 'Reserve error: insufficient free qubits for cell %s' %str(cell)
			return False
	
	
# unchecked
# affects tile_occ
def releaseQubits(cells):
	'''release reservations on adjacent cells for each cell in cells'''
	global _reserved, _qbit_flags
	
	for cell in cells:
		# for each reserved adjacent qbit unset reserved flag
		for qb in _reserved[cell]:
			_qbit_flags[qb]['reserved']=False
		# forget reserved set
		_reserved[cell].clear()
		
# unchecked
# affects tile_occ
def assignQubit(cell,qbit):
	'''assign a qubit to a cell and update flags and the number of
	unplaced adjacent cells. Qubit reservations not handled'''
	
	global _qubits, _numAdj, _source
	global _cell_flags, _qbit_flags, _tile_occ
	
	_qubits[cell]=qbit
	_cell_flags[cell]['placed']=True
	_qbit_flags[qbit]['taken']=True
	_qbit_flags[qbit]['reserved']=True # should not matter
	_qbit_flags[qbit]['assigned']=True
	
	# update tile_occ and vacancy
	_tile_occ['r'][qbit[0]]+=1
	_tile_occ['c'][qbit[1]]+=1
	setVancancy()
	
# unchecked, unfinished
# affects tile_occ
def assignPaths(paths):
	''' '''
	
	global _qbit_flags, _reserved, _qbitAdj, _qubits, _paths, _vacancy
	
	# we want to flag as taken each qbit in each chain
	# and update the reservations of any qbits adjacent to the chain
	
	for path in paths:
		for qbit in path[1:-1]:		# end qbits are already taken
			# flag chain qbit as taken
			_qbit_flags[qbit]['taken']=True
			# for each adjacent qubit with an assigned cell, update
			# reservations
			proxim_cells
			for adj_qb in _qbitAdj[qbit]:
				if _qbit_flags[adj_qb]['assigned']:
					
					
					
			
		
		
# unchecked	
def firstCell(numAdj2, M1=False):
	'''returns the first cell to be placed'''
	global _numAdj
	
	# create adjacency worths for each cell
	worth={key:(_numAdj[key],numAdj2[key]) for key in _numAdj}
	# sort by decreasing worth
	order=sorted(_numAdj.keys(),key=lambda x: worth[x])[::-1]
	
	### method 1: max adj
	
	if M1:
		# determine how many cells have the maximum worth
		num_max=worth.values().count(worth[order[0]])
		# randomly select one of these cells
		i=int(random()*num_max)
		return order[i]
		
	### method 2: fully probabilistic
	# probability is ~ numAdj**POW for some power

	else:
		# give a probability score for each cell
		probs={key:pow(worth[key][0],FIRST_PROB_POW) for key in worth}
		# normalise and compute comparison values
		total_prob=sum(probs.values())
		comps=[0]
		for key in order:
			probs[key]/=total_prob
			comps.append(comps[-1]+probs[key])
		# randomly select starting key
		i=max(bisect(comps,random()),1)
		
		return order[i-1]

# unchecked
def firstQubit(adj, M1=False):
	'''Selects the qubit corresponding to the first cell'''
	global _qbitAdj
	
	### method 1: middle cell
	
	if M1:
		
		# select candidate tile(s)
		n,m=[N//2],[M//2]
		if N%2==0: n.append(N//2-1)
		if M%2==0: m.append(M//2-1)
		tiles=[(_n,_m) for _n in n for _m in m]
		
		# shuffle tiles
		shuffle(tiles)
		
		for tile in tiles:
			r,c=tile
			# try to find suitable qubit
			order=[(h,i) for h in xrange(2) for i in xrange(L)]
			shuffle(order)
			for h,i in order:
				qbit=(r,c,h,i)
				if len(_qbitAdj[qbit])>=adj:
					return qbit
		return None
		
	
	### method 2: Gaussian dist
		
	else:
		
		if N%2: # if odd rows
			Y=np.arange(-(N//2),N//2+1)
		else:
			Y=.5+np.arange(-(N//2),N//2)
			
		if M%2: # if odd rows
			X=np.arange(-(M//2),M//2+1)
		else:
			X=.5+np.arange(-(M//2),M//2)
		
		# generate probabilities
		CDF=[]
		for ax in [X,Y]:
			Z=np.exp(-ax*ax/(2*FIRST_QBIT_SIG))
			Z/=np.sum(Z)
			cdf=[0.]
			for z in Z:
				cdf.append(cdf[-1]+z)
			CDF.append(cdf)
			
		
		# attempt to find qubit
		attempt=0
		while attempt<FIRST_QBIT_ATTEMPTS:
			attempt+=1
			# pick tile
			r=max(bisect(CDF[0],random()),1)-1
			c=max(bisect(CDF[1],random()),1)-1
			# pick qubit
			order=[(h,i) for h in xrange(2) for i in xrange(L)]
			shuffle(order)
			for h,i in order:
				qbit=(r,c,h,i)
				if len(qbitAdj[qbit])>=adj:
					return qbit
		return -1
	
	
	
			
	
# unchecked, unfinished
def placeCell(cell):
	'''Attempt to find a suitable qbit to place input cell on.
	
	inputs: cell(int)		: source index of cell to place
			
	output: qbit(tuple)		: 4-tup for qbit to assign to cell
			paths(list)		: list of paths from placed cell qubits to 
							  qbit
	'''
	
	global _qbitAdj, _source, _qubits
	global _cell_flags, _qbit_flags, _reserved, _vacancy
	
	### Initialise
	
	seam_flag=False
	qbit=None
	
	# find qubits for adjacent cells
	adj_qbits=[_qubits[c] for c in _source[cell] if _cell_flags[c]['placed']]
	
	# find required availability of target qbit
	avb=len(_source[cell])
	
	# every time a seam is opened, we should check to see if there is a
	# better qubit to consider
	
	while qbit is None:
		
	
		### Open Seam
		
		if seam_flag:
			seam_flag=False
			
			# check for vacancies
			if not any(vacancies):
				break
				
			# find available seams
			seams=set()
			for qb in adj_qbits:
				s=availableSeams(qb)
				seams=seams.union(s)
				
			# rank seams
			seam_ranks=map(rankSeam, seams)
			
			# select seam to open
			seam=selectSeam(seam_ranks)
			
			# open seam
			success=openSeam(seam)
			
			if not success:
				print 'Failed to open seam'
				return None,[]
			
			
		### Pick qubit to assign
		
		# run multisource search method
		qbit = multiSourceSearch(adj_qbits,avb)
		
		# check if found
		if qbit is None:
			seam_flag=True
			continue
		
		### Find paths
		
		routes=[[qb,qbit] for qb in adj_qbits]
		
		cost=Routing.Routing(routes,writePath=ROUTE_PATH)
		
		# check successfull routing
		
		if cost<Routing.COST_BREAK:
			qbit=None
			continue
		
	# get paths
	paths=cp(Routing.getPaths())
	
	return qbit,paths
		
		
# unchecked, unfinished
def forgetCell(cell):
	''' '''
	pass
	




### Multisource Search

# unchecked
def extend_Dijkstra(src):
	''' '''
	global _qbitAdj, _qbit_flags
	
	BIG_VAL=2*len(_qbitAdj)	# large value for initial node costs
	
	visited={qbit:False for qbit in _qbitAdj}
	costs={qbit:BIG_VAL for qbit in _qbitAdj}
	next_qb=set()
	next_qb.add(
	costs[src]=0
	
	while next_qb:
		
		# pick lowest cost qbit and yield
		qbit = sorted(next_qb,key=lambda x: costs[x])[0]
		next_qb.remove(qbit)
		yield qbit
		
		# mark as visited
		visited[qbit]=True
		
		# update costs of all unvisited adjacent nodes
		for qb in _qbitAdj[qbit]:
			if not visited[qb]:
				costs[qb]=min(costs[qb],costs[qbit]+1)
				next_qb.add(qb)
				
# unchecked, NOT IMPLEMENTED
def extend_Astar(src):
	''' '''
	global _qbitAdj, _qbit_flags
	
	print 'A* extension not yet implemented'
	sys.exit()
	
# unchecked
def multiSourceSearch(srcs,adj,typ='Dijkstra'):
	''' '''
	global _qbitAdj
	
	# create path extension generators
	
	if typ.upper()=='DIJKSTRA':
		extend_func=extend_Dijkstra
	else:
		extend_func=extend_Astar
		
	######## temporary override
	extend_func=extend_Dijkstra
	
	
	extend={src:extend_func(src) for src in srcs}
	
	visits={qbit:0 for qbit in _qbitAdj}
	
	# search loop
	
	while True:
		
		cands=[]	# candidate nodes
		
		for src in srcs:
			
			# extend 
			try:
				node=next(extend[src])
			except:
				print 'multiSourceSearch ERROR: No possible path'
				return None
			
			# increment visited node count
			visits[node]+=1
			
			# if visited from all sources add as candidate
			if visits[node]==len(srcs):
				cands.append(node)
				
		## check for suitable candidate
		
		# sort by suitability
		cands=sorted(map(lambda x: [suitability(x),x],cands))[::-1]
		
		# filter out unsuitable qbit
		cands=filter(lambda x: x[0]>=adj,cands)
		
		# select qbit
		if cands:
			node=cands[0][1]	# simple version: change later
			return node
			
	
	
### Seam Splitting

# seam indexing format [h,n]
	# ex. seam between columns 3 and 4 -> 	[0,3]
	#	  seam between rows 2 and 3 ->		[1,2]
	
# unchecked, unfinished
def availableSeams(qbit):
	'''Returns a set of possible seams to split: [seam_index, direct]'''
	global _vacancy
	
	seams=[]
	
	# left/right seams
	
	if _vacancy[0]>0 and qbit[1]>0: 	# open left
		seams.append([ [0,qbit[1]-1] ,False])
		seams.append([ [0,qbit[1] ,False)
		
	if _vacancy[1]>0 and qbit[1]<(N-1):	# open right
		seams.append([ [0,qbit[1]-1] ,True])
		seams.append([ [0,qbit[1]] ,True])
	
	# down/up seams
	
	if _vacancy[2]>0 and qbit[0]>0: 	# open down
		seams.append([ [1,qbit[0]-1] ,False])
		seams.append([ [1,qbit[0]] ,False])
		
	if _vacancy[3]>0 and qbit[0]<(M-1):	# open up
		seams.append([ [1,qbit[0]-1] ,True])
		seams.append([ [1,qbit[0]] ,True])
		
	return seams
	

# unchecked, unfinished	
def rankSeam(seam):
	'''Give seam a score based on the number of conflicts which arise
	from potentially opening each seam (filters bad seams)'''
	
	global _qbit_flags, _paths
	
	sm,dr=seam	# seam and open direction
	
	## get list of assigned qbits and connectors to be moved
	
	qbits=[]
	connects=[]
	
	def check_qb(qb):
		if dr:
			return qb[sm[0]]>sm[1]
		else:
			return qb[sm[0]]<=sm[1]
		
	def target_qb(qb):
		q=cp(qb)
		if dr:
			q[sm[0]]+=1
		else:
			q[sm[0]]-=1
		return q 
		
	# list of qubits to be moved
	for qb in filter(check,_qbit_flags):
		if qb['assigned']:
			qbits.append(qb)
	
	# list of connectors to be 'moved'
	for path in _paths:	# for each connector on each path
		for i in xrange(1,len(path)):	# if either qubit is to be moved
			connect=sorted([path[i],path[i-1]])
			if not any(map(check_qb,connect)):
				connects.append(connect)
			
	## detect qbit conflicts
	
	qbit_conflicts=0
	
	# for each qbit to be moved, check target qbit
	for qb in map(target_qb,qbits):
		try: 
			t=_qbit_flags[qb]
		except KeyError:
			qbit_conflicts+=1
	
	## detect connector conflicts
	
# unchecked, unfinished
def openSeam(seam):
	''' '''
	global _paths
	
	damage=[]	# list of embedded 
	
	## Open seam
	
	
	
	
	# all flags, etc. must be fully updated to allow recursion
	
	## Repair damage
	
	# select cells to forget
	damage=damage	# modify to select subset of cells needed to replace
	
	# forget cells
	for cell in damage:
		forgetCell(cell)
	
	# attempt to find new placements
	for cell in damage:
		
		qbit,paths=placeCell(cell)
		
		if qbit is None:
			print 'Recursive seam failure...'
			return False
			
		# assign qubit and paths
		assignQubit()
		assignPaths(paths) # handles reservations
		
		_paths+=paths
		
		
		
		
		
	
	
#######################################################################
#######################################################################
### MAIN ###

# unchecked			
def denseEmbed(source):
	'''
	Attempts to find an embedding of the source graph in the target
	graph.
	
	inputs:	source (dict)	: adjacency dict: source graph
			coupler_flags(dict): from getCouplerFlags.
			
	outputs: qubits (dict) 	: source node indexed mapping of assigned
							 qubits
			 routes	(dict) 	: (node1,node2) indexed dictionary of qubit
			 				 routes; node1 < node2
			 info (dict)	: describe later...
	'''
	global _numAdj, _qbitAdj, _source, _cell_flags
	
	### INITIALIZE ###
	
	_source=source
	
	# generate numAdj from source
	_numAdj={key:len(source[key]) for key in source}
	
	# generate numAdj2
	numAdj2={}	# sum of numAdj over each cells adjacent to key
	for key in source:
		numAdj2[key]=sum([_numAdj[adj] for adj in source[key]])
	
	# initialize flags and reserved dicts
	setFlags()
	
	# configure routing algorithm
	routing.initialize(_qbitAdj)
	
	### INITIAL SEED ###
	
	# select first cell
	cell=firstCell(_numAdj,numAdj2)
	
	# select first qubit
	qbit=firstQubit(_numAdj[cell],_qbitAdj)
	
	if qbit is None:
		print 'No suitable first qubit found...'
		sys.exit()
	
	# take first qubit
	assignQubit(cell,qbit)
	
	# handle reservations
	reserveQubits(cell,qbit)
	
	# update do*
	doNow=sorted(source[qbit])
	doNext.clear()
	
	
	### GENERAL PLACEMENT ###
	
	while doNow:
		
		# place all cells in doNow
		
		for cell in doNow:
			
			# find qbit and paths to qbit from placed cells
			qbit,paths=placeCell(cell)
			
			# abort on failed placement
			if qbit is None:
				print 'No placement of cell %s found' %str(cell)
				return -1
				
			# assign qubit and paths
			assignQubit()
			assignPaths(paths) # handles reservations
			
			_paths+=paths
			
			# add unplaced adjacent cells to doNext
			for c2 in source[cell]:
				if not _cell_flags[c2]['placed']:
					doNext.add(c2)
			
		
		# update do*
		
		doNow=sorted(doNext,key=lambda x: numAdj[x])
		doNext.clear()
		
	
if __name__ == '__main__':
	pass
