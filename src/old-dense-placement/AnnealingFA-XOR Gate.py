import random
import math
import RoutingFunctionv6
import DenseInitialPlacement
import ExpandFunctions
import sys
import time
import os

numberOfQubits = 512
qubitsPerCell = 8
triesPerSize = 10
minGridSize = 7
maxGridSize = 20
solutionsPerCircuit=10

filedir='./circuits/'
outdir='./bin/'

couplers='couplers.txt'
#filename = 'S_R_Flip_Flop_bench'
#filename = 'selectable_oscillator_bench'
filename = 'xor_gate_bench'
#filename = 'Full-Adder3' #The one we use
#filename = 'threeMajorityGate'
#filename = 'Full-Adder2'
#filename = 'Full-Adder'
#filename = 'majorityGateInverter'
#filename = 'majorityGate'
#filename = 'Inverter'
#filename = 'serial_adder_bench'
#filename = 'walus_loop_memory_cell_bench'
#filename = '4-bit_2-1_mux_bench'
#filename = '4_bit_accumulator_bench'

file1 = open(filedir+filename, 'r')
data =file1.readlines()
file1.close()

os.chdir(outdir)
routefile=open('routes.txt','w')

#print filename

#if numberOfQubits!=512:
#    RoutingFunctionv6.setNumberQubits(numberOfQubits)
#    couplers='couplers'+str(math.pow(numberOfQubits/8,0.5))+'.txt'

# Tuning Knobs. Higher numbers will yield a better annealing process
tempRange = 3
changesPerTemp = 2
numberPossibleMoves = 3
maxAttempts = 10


###### Lists ######

cells = []
cells_ind = []
cells_x = []
cells_y = []
newCells = []
position_x = []
position_y = []
newPosition_x = []
newPosition_y = []
distance = []
distance_temp = []
angle = []
angle_temp = []
angle_modified = []
angle_modified_temp = []
driver = []
initial_pol = []
new_pol = []
adjacency = []
temp_adj = []
qubits = []
temp_qubits = []
history = []
route = []
number = []

grid=[]
grid_x=[]
grid_y=[]
grid_position_x=[]
grid_position_y=[]
##### Lists #####

##### Distances and Angles #####

for i in range(0,len(data)):
    #if data[i] == '[TYPE:QCADCell]\n':
    if data[i] == '[TYPE:QCADCell]\r\n':
        cells.append(i)

    ### Finding Divers ###

    #if data[i] == 'cell_function=QCAD_CELL_NORMAL\n':
    if data[i] == 'cell_function=QCAD_CELL_NORMAL\r\n':
        driver.append(0)
    #elif data[i] == 'cell_function=QCAD_CELL_OUTPUT\n':
    elif data[i] == 'cell_function=QCAD_CELL_OUTPUT\r\n':
        driver.append(0)
    #elif data[i] == 'cell_function=QCAD_CELL_FIXED\n':
    elif data[i] == 'cell_function=QCAD_CELL_FIXED\r\n':
        driver.append(1)

    ### Finding Drivers ###

##### Finding Grid Spacing #####

spacing = data[28]
spacing = float(spacing[13:len(spacing)-2])
#print(spacing)

##### Finding Grid Spacing #####

for i in range(0,len(cells)):
    cells_x.append(data[cells[i]+2])
    cells_y.append(data[cells[i]+3])
    temp_x = cells_x[i]
    temp_y = cells_y[i]
    position_x.append(float(temp_x[2:len(temp_x)-2]))
    position_y.append(float(temp_y[2:len(temp_y)-2]))

sizex = (max(position_x)-min(position_x))/spacing + 1
sizey = (max(position_y)-min(position_y))/spacing + 1

 ###FINDING GRID PLACEMENT###

file3 = open(filename +'_Relations.txt','w')
#file3.write(str(int(sizex)) + ' '  + str(int(sizey)) + '\n')
file3.write(str(int(sizex)) + ' '  + str(int(sizey)) + '\r\n')



a = min(position_y)

while a <= max(position_y):
    b = min(position_x)
    temp_newCells =[]
    temp_position_x = []
    for i in range(len(cells)):
        if position_y[i] == a:
            temp_newCells.append(i)
            temp_position_x.append(position_x[i])
    while b <= max(position_x):
        for j in range(len(temp_newCells)):
            if temp_position_x[j] == b:
                newCells.append(temp_newCells[j])
        b = b + spacing
##    newCells.append(temp_newCells)
    a = a + spacing


for i in range(len(cells)):
    newPosition_x.append(position_x[newCells[i]])
    newPosition_y.append(position_y[newCells[i]])


position_x = []
position_y = []

for i in range(len(cells)):
    position_x.append(newPosition_x[i])
    position_y.append(newPosition_y[i])

for i in range(len(cells)):
    number.append(int((position_x[i]-min(position_x))/spacing + ((position_y[i] - min(position_y))/spacing)*sizex))

#print number

for i in range(0,len(cells)):
    distance_temp = []
    angle_temp = []
    for j in range(0,len(cells)):
        temp_x = position_x[i] - position_x[j]
        temp_y = position_y[i] - position_y[j]
        distance_temp.append(math.sqrt(math.pow(temp_x,2) + math.pow(temp_y,2)))

        if j != i:
            angle_temp.append(math.acos(temp_x /math.sqrt(math.pow(temp_x,2) + math.pow(temp_y,2))))
        else:
            angle_temp.append(0)

        if j == len(cells)-1:
            distance.append(distance_temp)
            angle.append(angle_temp)


for i in range(len(angle)):
    angle_modified_temp = []
    angle_temp = angle[i]

    for j in range(len(angle_temp)):
        angle_modified_temp.append(math.cos(4 * angle_temp[j]))

        if j == len(angle_temp) - 1:
            angle_modified.append(angle_modified_temp)

##### Distances and Angles #####
### Initalize Grid ####

for i in range(int(sizey)):
    grid.append([0]*int(sizex))
    grid_x.append([None]*int(sizex))
    grid_y.append([None]*int(sizex))

### Initalize Grids ####

#uses spacing to find x and y value of each cell
for i in range(len(cells)):
    grid_position_x.append(int((position_x[i] - min(position_x))/spacing))
    grid_position_y.append(int((position_y[i] - min(position_y))/spacing))

#print grid_position_x
#print grid_position_y


###place in appropriate grids
for i in range(len(cells)):
    temp_grid = grid[grid_position_y[i]]
    temp_grid[grid_position_x[i]] = i+1
    temp_x = grid_x[grid_position_y[i]]
    temp_x[grid_position_x[i]] = position_x[i]
    temp_y = grid_y[grid_position_y[i]]
    temp_y[grid_position_x[i]] = position_y[i]

### Initialize Grid ###
##### Drivers ######

for i in range(len(cells)):


    if driver[i] == 0:
        initial_pol.append(0)

    else:
        if i == len(cells)-1:
            for j in range(len(data)-cells[i]):
                try:
                    if float(data[cells[i]+j][7:12]) > 1.55 and float(data[cells[i]+j][7:12]) < 1.65:
                        initial_pol.append(1)
                        break
                    elif float(data[cells[i]+j][7:12]) == 0:
                        initial_pol.append(-1)
                        break
                except ValueError:
                    pass
        else:
            for j in range(cells[i+1]-cells[i]):
                try:
                    if float(data[cells[i]+j][7:12]) > 1.55 and float(data[cells[i]+j][7:12]) < 1.65:
                        initial_pol.append(1)
                        break
                    elif float(data[cells[i]+j][7:12]) == 0:
                        initial_pol.append(-1)
                        break
                except ValueError:
                    pass

for i in range(len(cells)):
    new_pol.append(initial_pol[newCells[i]])

file2 = open(filename + '_Route.txt','w')
##print 'Simulating input circuit:',initial_pol, 'from file',filename,'...'
#file2.write('Routing input circuit:' + str(new_pol) + ' from file ' + filename + '\n')
file2.write('Routing input circuit:' + str(new_pol) + ' from file ' + filename + '\r\n')




for i in range(len(cells)):
    temp_adj = []
    distance_temp = distance [i]

    #Full Adjacency
    for j in range(i,len(cells)):
        if distance_temp[j] == spacing or distance_temp[j] == spacing*pow(2,0.5):
            temp_adj.append(j)

    #Adjacency with diagonal only at inverters
##    for j in range(len(cells)):
##        if distance_temp[j] == spacing:
##            try:
##                if i in adjacency[j]:
##                    continue
##                else:
##                    temp_adj.append(j)
##            except IndexError:
##                temp_adj.append(j)

##    x=grid_position_x[i]
##    y=grid_position_y[i]
##    for k in [-1,1]:
##        for m in [-1,1]:
##            try:
##                if grid[y+k][x+m]>0 and grid[y+k][x]==0 and grid[y][x+m]==0:
##                    print "diag at:",i,"k: ",k,"m: ",m
##                    temp_diag1=(grid[y+k][x+m]-1) #because grid incremented by 1
##                    if temp_diag1<i and not(i in adjacency[temp_diag1]):
##                        adjacency[temp_diag1].append(i)
##                    elif temp_diag1>i:
##                        temp_adj.append(temp_diag1)
##            except IndexError:
##                continue

    adjacency.append(temp_adj)

#print adjacency

def select(temp,state1,old,new):

    cost2 = []
    prob = []
    for i in range(0,len(old)):
        cost2.append(RoutingFunctionv6.modifyRoute(old[i],new[i],0))
        if cost2[i] > 2*state1 :
            prob.append(0)
        else:
            if cost2[i] != 9999:
                prob.append(math.exp(round((state1-cost2[i])/(temp),5)))
            else:
                prob.append(0)
    if sum(prob) == 0:
        return -1
    for i in range(len(prob)):
        prob[i] = prob[i]/sum(prob)
    while True:
        choice = random.randint(0,len(old)-1)
        if prob[choice] != 0:
            break
    #file2.write('Cost = ' + str(cost2) + '\n')
    #file2.write('Probablity = ' + str(prob) + '\n')

    if cost2[choice] < state1:
        RoutingFunctionv6.modifyRoute(old[choice],new[choice],1)
        return choice
    elif random.uniform(0,1) < prob[choice]:
        RoutingFunctionv6.modifyRoute(old[choice],new[choice],1)
        return choice
    else:
        return -1

def check(totalList,newlist):
    # Returns true if a list in in another list, false if it is not.
     for i in range(len(totalList)):
        temp_hist = totalList[i]
        for j in range(len(temp_hist)):
            if temp_hist[j] != newlist[j]:
                break
            elif j == len(temp_hist)-1 and temp_hist[j] == newlist[j]:
                return True
     return False

if False:
	for cell in xrange(len(adjacency)):
		print '%d: %s' %(cell, str(adjacency[cell]))
		
#Establish Connections
connections=DenseInitialPlacement.connections(adjacency)

#Reorder Adjacency
newAdjacency,order,numberAdjacent=DenseInitialPlacement.reorderVer2(adjacency,connections)

successes=0
solution=open(str(filename)+'Solution'+str(successes)+'.txt','w')
specs=open(str(filename)+'Specs'+'.txt','w')
specs.write('Tries Per Size: '+str(triesPerSize)+'\r\n')
specs.write('Maximum Grid Size: '+str(maxGridSize)+'\r\n\r\n\r\n\r\n')
for i in range(minGridSize,maxGridSize):
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
                                print ''
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
                solution=open(str(filename)+'Solution'+str(successes)+'.txt','w')
                solution.write('new:\r\n')
                for j in range(len(routeSolution)):
                    for k in range(len(routeSolution[j][:-1])):
                        solution.write(str(routeSolution[j][k])+' ')
                    solution.write(str(routeSolution[j][-1])+'\r\n')
                solution.write('shared:\r\n')
                solution.close()
                specs.write(str(filename)+'Solution'+str(successes)+' Specifications:\r\n')
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
            if successes >= solutionsPerCircuit:
                done=True
                break
        else:
            #print count,'consecutive fails'
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
                            print ''
                            gridLess+=1
                            vacancy[3]-=1
                            vacancy[2]-=1
                            break
                    adjusted=True
                    break
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

specs.close()
#file2.write('\n----------------------------------------\n')
file2.write('\r\n----------------------------------------\r\n')
##print '\n ---------------------------------------------'
#file2.write('\nFinal Placement:\n')
file2.write('\r\nFinal Placement:\r\n')
##print 'Final Placement:'
for i in range(len(qubits)):
    #file2.write('Cell ' + str(number[i]) + ' placed on qubit ' + str(qubits[i]) + '\n')
    file2.write('Cell ' + str(number[i]) + ' placed on qubit ' + str(qubits[i]) + '\r\n')
##    print 'Cell ' + str(i) + ' placed on qubit ' + str(qubits[i])
    #file3.write(str(qubits[i]) + ' ' + str(number[i]) + '\n')
    file3.write(str(qubits[i]) + ' ' + str(number[i]) + '\r\n')

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
file3.close()
#execfile("RoutingViewerSolutionTest.py",{'numberOfQubits':numberOfQubits})
