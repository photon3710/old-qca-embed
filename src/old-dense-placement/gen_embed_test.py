#!usr/bin/python


from __future__ import division

import time,sys
sys.path.append('./dense-placement-old/')
import run_dense_embed as demb
import numpy as np
from dwave_sapi import find_embedding, get_chimera_adjacency
from generator import generateCircuit, generateCircuit2, GtoCoef
import pylab as plt

RUN_DENSE=True

NUM_RUNS=200
NUM_TRIALS=5
TIMEOUT=10	# seconds


def procDenseEmbed(embed):
	''' '''
	
	qbits,routes=embed
	chain_lengths=[len(route) for route in routes]
	
	used_qbits=set()
	for route in routes:
		for qb in route:
			used_qbits.add(qb)
	
	qubits=len(used_qbits)
	
	return qubits,chain_lengths
	
	
	
def procHeurEmbed(embed):
	'''Process an embedding from the heuristic algorithm'''
	
	num_qubits=0
	chain_lengths=[]
	flags=[False]*512
	fail=False
	
	for chain in embed:
		if fail:
			break
		for qbit in chain:
			if flags[qbit]:
				fail=True
				break
			flags[qbit]=True
		chain_lengths.append(len(chain))
		
	if fail:
		return -1, []
	else:
		return sum(chain_lengths), chain_lengths


def formatCoef(h,J):
	'''Format coefficient matrices for input to heuristic algorithm'''
	
	S={}
	N=len(h)	
	
	for i in xrange(N):
		for j in xrange(N):
			if abs(J[i,j]) > 0.01:
				S[(i,j)]=1
			else:
				S[(i,j)]=0
	
	return S, N
	
		
def runDense(h,J, trials):
	'''run the dense placement method given the  '''
	
	# compute adjacency
	adjacency=[]
	N=len(h)
	for i in xrange(N):
		adjacency.append([j for j in xrange(i+1,N) if abs(J[i][j])>0.01])
	
	if False:
		for i in xrange(N):
			print '%d: %s' %(i,adjacency[i])
		
	solutions=demb.run(adjacency,trials,tofile=False)
	
	good_embeds=[]
	times=[]
	nums=[len(solutions), trials]
	
	for sol in solutions:
		good_embeds.append([sol['qubits'],sol['routes']])
		times.append(sol['time'])

	return good_embeds,times,nums
	
		
	
def runHeuristic(S,S_size, max_count, flagSol=False):
	'''run the heuristic method'''
	
	M = 8
	N = M
	L = 4
	
	A = get_chimera_adjacency(M,N,L)
	A_size = M * N * L * 2
	
	trial_num = 0
	success_num = 0
	count = 0
	
	good_embeds = []
	times = []
	
	while count < max_count:
		
		t1 = time.clock() 
		embeddings = find_embedding(S, S_size, A, A_size, verbose=0, tries=1, timeout=TIMEOUT)
		t2 = time.clock()
		trial_num += 1
		
		if len(embeddings) == S_size:	# successful embedding
			success_num += 1
			good_embeds.append(embeddings)
			print 'solution ' + str(success_num) + ' found...'
			
			times.append(t2-t1)
		
		count = success_num if flagSol else trial_num
	
	print '\n'*2+'Embedded ' + str(success_num) + ' of ' + str(trial_num) + ' attempts' + '\n'*2
	
	return good_embeds, times, [success_num,trial_num]
	
	
def main():
	
	
	for full_adj in [False]:	
		D = []
		print ('\n'+'*'*50)*2 + '\n** ',
		print 'FULL ADJACENCY' if full_adj else 'LIM ADJACENCY'
		
		for run in xrange(NUM_RUNS):
		
			print '\n'+'*'*40 + "\nRun %d\n\n" %(run+1)
			G=generateCircuit2(full_adj=full_adj)
			h,J=GtoCoef(G)
			if RUN_DENSE:
				good_embeds, times, nums = runDense(h,J,trials=NUM_TRIALS)
			else:
				S,S_Size=formatCoef(h,J)
				good_embeds, times, nums = runHeuristic(S,S_Size,max_count=NUM_TRIALS,flagSol=False)
			prob = nums[0]/nums[1]
		
			QUBITS=[]
		
			if good_embeds:
				for embed in good_embeds:
					if RUN_DENSE:
						qubits,chain_lengths=procDenseEmbed(embed)
					else:
						qubits,chain_lengths=procHeurEmbed(embed)
			else:
				qubits,chain_lengths=-1,[]
				
			QUBITS.append(qubits)
			print 'Success Rate: %f' %prob
		
			if not QUBITS:
				QUBITS=[-1]
			
			if QUBITS:
				mean_qubits = np.mean(QUBITS)
				min_qubits = np.min(QUBITS)
		
				D.append([len(h),mean_qubits,min_qubits])
			
	
		# mean qubits
	
		c=['g','r'] if full_adj else ['b','m']
		X,Y=[],[]
		for d in D:
			if d[1]==-1:
				continue
				plt.plot(d[0],0,c[1]+'x',markersize=5, markeredgewidth=2)
			else:
				X.append(d[0])
				Y.append(d[1])
		plt.plot(X,Y,c[0]+'x',markersize=5, markeredgewidth=2)
	
	plt.legend(['Limited Adjacency'],numpoints=1,loc='upper left')
	plt.xlabel('Number of Cells')
	plt.ylabel('Average Qubit Usage')
	#plt.title('Average Qubit Usage vs. Number of Cells')
	plt.show()
	
	# min qubits
	
	#for d in D:
	#	if d[2]==-1:
	#		plt.plot(d[0],0,'rx')
	#	else:
	#		plt.plot(d[0],d[2],'kx')
	
	#plt.xlabel('Number of Cells')
	#plt.ylabel('Minimum Qubit Usage')
	#plt.title('Minimum Qubit Usage vs. Number of Cells')
	#plt.show()
	
	
if __name__ == '__main__':
	main()
