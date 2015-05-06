#!usr/bin/python


import sys, os, time

sys.path.append('../')

from parse_qca import parseQCAFile
from aux import generateAdjDict, convertToNearestNeighbour

import DenseInitialPlacement
import ExpandFunctions
import RoutingFunctionv6


### GLOBALS

outdir='./bin/'

numberOfQubits = 512
qubitsPerCell = 8
triesPerSize = 10
minGridSize = 8
maxGridSize = 8
solutionsPerCircuit=10

FULL_ADJ=True
TO_FILE=False

def formatAdjacency(adjDict):
	''' '''
	
	cells=sorted(adjDict.keys())
	adjacency=[]
	cell_map={cells[i]:i for i in xrange(len(cells))}
	
	for key in cells:
		adj=map(lambda x:cell_map[x[0]],adjDict[key][1::])
		adjacency.append(adj)
		
	return adjacency
	
def setDir(filename):
	'''Erase old data, confirm directory structure'''
	
	global FULL_ADJ
	
	root=filename.rpartition('/')[2].replace(' ','-')
	
	os.chdir('./bin')
	ls=os.listdir('.')
	if not root in ls:
		os.mkdir(root)
	os.chdir(root)
	ls=os.listdir('.')
	d='Full' if FULL_ADJ else 'Limited'
	if not d in ls:
		os.mkdir(d)
	os.chdir(d)
	for f in os.listdir('.'):
		os.remove(f)
	
	
	

def run(adjacency, num_trials=None, tofile=True):
	''' '''
	
	global FULL_ADJ
	
	SOLUTIONS=[]
	
	if False:
		for cell in xrange(len(adjacency)):
			print '%d: %s' %(cell, str(adjacency[cell]))
	
	if tofile:
		routefile=open('routes.txt','w')
		file2 = open('Route.txt','w')
	#file2.write('Routing input circuit:' + str(new_pol) + ' from file ' + filename + '\r\n')

	#Establish Connections
	connections=DenseInitialPlacement.connections(adjacency)

	#Reorder Adjacency
	newAdjacency,order,numberAdjacent=DenseInitialPlacement.reorderVer2(adjacency,connections)

	successes=0
	
	if tofile:
		solution=open('Solution'+str(successes)+'.txt','w')
		specs=open('Specs'+'.txt','w')
		specs.write('Tries Per Size: '+str(triesPerSize)+'\r\n')
		specs.write('Maximum Grid Size: '+str(maxGridSize)+'\r\n\r\n\r\n\r\n')
		
	trial=0
	for i in range(minGridSize,maxGridSize+1):
		#print 'Using', str(i)+'x'+str(i), 'qubit grid from now on'
		tStart=time.clock()
		count=0
		done=False
		numberOfQubits=i*i*qubitsPerCell
		DenseInitialPlacement.setNumQubits(numberOfQubits)
		ExpandFunctions.genCouplers(i)
		couplers='couplers'+str(i)+'.txt'
		RoutingFunctionv6.setNumberQubits(numberOfQubits)
		#Initial Placement
		
		while True:
			if num_trials and trial >= num_trials:
				break
			t0=time.clock()
			qubits,tQ,rows,columns,routeSolution = DenseInitialPlacement.placeQubits(newAdjacency,order,numberAdjacent,connections)
			t1=time.clock()
			vacancy=DenseInitialPlacement.perimeterVacancy(tQ,rows,columns)
			vacancyCopy=list(vacancy)
			takenQubits=[]
			for j in range(len(tQ)):
				if tQ[j]==1:
					takenQubits.append(j)
			tiles=0
			for j in range(i*i):
				for k in range(qubitsPerCell):
					if tQ[j*qubitsPerCell+k]==1:
						tiles+=1
						break
			couplersUsed=0
			for j in range(len(routeSolution)):
				couplersUsed+=len(routeSolution[j])-1
			gridHeight=i-vacancyCopy[0]-vacancyCopy[3]
			gridWidth=i-vacancyCopy[1]-vacancyCopy[2]
			if qubits != []:
				if i != 0:
					adjusted=False
					rowDifference=i*qubitsPerCell
					gridLess=0
					while not adjusted:
						while True:
							if vacancy[0]>0:
								if vacancy[1]>0:
									#print 'shift up,left'
									for j in range(1,i):
										for k in range(len(qubits)):
											if qubits[k] in rows[j]:
												qubits[k]-=rowDifference
										for k in range(len(routeSolution)):
											for l in range(len(routeSolution[k])):
												if routeSolution[k][l] in rows[j]:
													routeSolution[k][l]-=rowDifference

									for j in range(1,i):
										for k in range(len(qubits)):
											if qubits[k] in columns[j]:
												qubits[k]-=qubitsPerCell
										for k in range(len(routeSolution)):
											for l in range(len(routeSolution[k])):
												if routeSolution[k][l] in rows[j]:
													routeSolution[k][l]-=qubitsPerCell

									gridLess+=1
									vacancy[0]-=1
									vacancy[1]-=1
									break
								if vacancy[2]>0:
									#print 'shift up'
									for j in range(1,i):
										for k in range(len(qubits)):
											if qubits[k] in rows[j]:
												qubits[k]-=rowDifference
										for k in range(len(routeSolution)):
											for l in range(len(routeSolution[k])):
												if routeSolution[k][l] in rows[j]:
													routeSolution[k][l]-=rowDifference
									gridLess+=1
									vacancy[0]-=1
									vacancy[2]-=1
									break
							if vacancy[3]>0:
								if vacancy[1]>0:
									#print 'shift left'
									for j in range(1,i):
										for k in range(len(qubits)):
											if qubits[k] in columns[j]:
												qubits[k]-=qubitsPerCell
										for k in range(len(routeSolution)):
											for l in range(len(routeSolution[k])):
												if routeSolution[k][l] in rows[j]:
													routeSolution[k][l]-=qubitsPerCell
									gridLess+=1
									vacancy[3]-=1
									vacancy[1]-=1
									break
								if vacancy[2]>0:
									#print ''
									gridLess+=1
									vacancy[3]-=1
									vacancy[2]-=1
									break
							adjusted=True
							break
					if gridLess > 0:
						numberOfQubits=(i-gridLess)*(i-gridLess)*qubitsPerCell
						for j in range(len(qubits)):
							qubitGroup=int(qubits[j]/qubitsPerCell)
							qubitIndex=qubits[j]%qubitsPerCell

							row=int(qubitGroup/i)
							col=int(qubitGroup%i)

							newQubitGroup=((i-gridLess)*row+col)*qubitsPerCell
							qubits[j]= int(newQubitGroup+qubitIndex)
						for j in range(len(routeSolution)):
							for k in range(len(routeSolution[j])):
								qubitGroup=int(routeSolution[j][k]/qubitsPerCell)
								qubitIndex=routeSolution[j][k]%qubitsPerCell

								row=int(qubitGroup/i)
								col=int(qubitGroup%i)

								newQubitGroup=((i-gridLess)*row+col)*qubitsPerCell
								routeSolution[j][k]= int(newQubitGroup+qubitIndex)

					RoutingFunctionv6.setRoute(routeSolution,couplers)
					if tofile:
						solution=open('Solution'+str(successes)+'.txt','w')
						solution.write('new:\r\n')
						for j in range(len(routeSolution)):
							for k in range(len(routeSolution[j][:-1])):
								solution.write(str(routeSolution[j][k])+' ')
							solution.write(str(routeSolution[j][-1])+'\r\n')
						solution.write('shared:\r\n')
						solution.close()
						specs.write('Solution'+str(successes)+' Specifications:\r\n')
						specs.write('Starting Grid Size: '+str(i)+'\r\n')
						specs.write('Grid Size: '+str(i-gridLess)+'\r\n')
						specs.write('Used Height: '+str(gridHeight)+'\r\n')
						specs.write('Used Width: '+str(gridWidth)+'\r\n')
						specs.write('Qubits Used: '+str(len(takenQubits))+'\r\n')
						specs.write('Couplers Used: '+str(couplersUsed)+'\r\n')
						specs.write('Tiles Used: '+str(tiles)+'\r\n')
						specs.write('Time (seconds): '+str(t1-t0)+'\r\n')
						specs.write('Fails Before Solution was Found: '+str(count)+'\r\n\r\n\r\n\r\n')
					successes+=1
					#print successes, 'success'
					tElapsed=time.clock()-tStart
					seconds=int(tElapsed%60)
					minutes=int((tElapsed/60)%60)
					hours=int(tElapsed/216000)
					#print 'Elapsed Time ',hours,':',minutes,':',seconds
					count=0
					
					temp={}
					temp['qubits']=qubits
					temp['routes']=routeSolution
					temp['time']=t1-t0
					SOLUTIONS.append(temp)
					trial+=1
				if successes >= solutionsPerCircuit:
					done=True
					break
			else:
				#print count,'consecutive fails'
				if i==maxGridSize-1:
					trial+=1
				adjusted=False
				gridLess=0
				while not adjusted:
					while True:
						if vacancy[0]>0:
							if vacancy[1]>0:
								gridLess+=1
								vacancy[0]-=1
								vacancy[1]-=1
								break
							if vacancy[2]>0:
								gridLess+=1
								vacancy[0]-=1
								vacancy[2]-=1
								break
						if vacancy[3]>0:
							if vacancy[1]>0:
								gridLess+=1
								vacancy[3]-=1
								vacancy[1]-=1
								break
							if vacancy[2]>0:
								#print ''
								gridLess+=1
								vacancy[3]-=1
								vacancy[2]-=1
								break
						adjusted=True
						break
				if tofile:
					specs.write('Fail:\r\n')
					specs.write('Starting Grid Size: '+str(i)+'\r\n')
					specs.write('Grid Size: '+str(i-gridLess)+'\r\n')
					specs.write('Used Height: '+str(gridHeight)+'\r\n')
					specs.write('Used Width: '+str(gridWidth)+'\r\n')
					specs.write('Qubits Used: '+str(len(takenQubits))+'\r\n')
					specs.write('Couplers Used: '+str(couplersUsed)+'\r\n')
					specs.write('Tiles Used: '+str(tiles)+'\r\n')
					specs.write('Time (seconds): '+str(t1-t0)+'\r\n\r\n\r\n\r\n')
				tElapsed=time.clock()-tStart
				seconds=int(tElapsed%60)
				minutes=int((tElapsed/60)%60)
				hours=int(tElapsed/216000)
				#print 'Elapsed Time ',hours,':',minutes,':',seconds
				#print 'Grid Size=',i
				if count >= triesPerSize:
					#print '\n\n\n\n\n\nWON\'T WORK!!!\n\n\n\n\n'
					#if i == maxGridSize-1:
					#    while True:
					#        stop=True
					break
				count+=1
		if done:
			#print 'SUCCESS!!!\nRouted on '+str(numberOfQubits)+' qubit grid'
			#print 'Failed',count, 'times'
			break

	if tofile:
		specs.close()
		#file2.write('\n----------------------------------------\n')
		file2.write('\r\n----------------------------------------\r\n')
		##print '\n ---------------------------------------------'
		#file2.write('\nFinal Placement:\n')
		file2.write('\r\nFinal Placement:\r\n')
		##print 'Final Placement:'
		for i in range(len(qubits)):
			#file2.write('Cell ' + str(number[i]) + ' placed on qubit ' + str(qubits[i]) + '\n')
			file2.write('Cell ' + str(i) + ' placed on qubit ' + str(qubits[i]) + '\r\n')
		##    print 'Cell ' + str(i) + ' placed on qubit ' + str(qubits[i])
			#file3.write(str(qubits[i]) + ' ' + str(number[i]) + '\n')
			#file3.write(str(qubits[i]) + ' ' + str(number[i]) + '\r\n')

		"""IMPORTANT"""
		routefile.close()
		"""IMPORTANT"""

	
		#file2.write('\nRoutes:\n')
		file2.write('\r\nRoutes:\r\n')
		##print '\nRoutes:'
		#file2.write(str(RoutingFunctionv6.getRoute())+'\n')
		file2.write(str(RoutingFunctionv6.getRoute())+'\r\n')
		##print RoutingFunctionv4.getRoute()
		#file2.write('\n--------------- JOB DONE ---------------\n')
		file2.write('\r\n--------------- JOB DONE ---------------\r\n')
		file2.close()
		#file3.close()
		
	return SOLUTIONS

if __name__ =='__main__':
	
	fname=sys.argv[1]
	
	if len(sys.argv)>2:
			FULL_ADJ=int(sys.argv[2])==0
	
	# load circuit
	cells,spacing = parseQCAFile(fname)
	
	setDir(fname)
	
	# generate adjacency
	adjDict,driver_index=generateAdjDict(cells,spacing)
	if not FULL_ADJ:
		adjDict=convertToNearestNeighbour(adjDict)
	
	adjacency=formatAdjacency(adjDict)
	
	# run embedding
	solutions=run(adjacency,10,tofile=TO_FILE)
	
	for sol in solutions:
		print '*'*40+'\n\n'
		print sol['time']
		print sol['qubits']
		print sol['routes']
