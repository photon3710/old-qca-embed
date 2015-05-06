#!usr/bin/python

from random import randint, random
import math
import RoutingFunctionv6

maxPaths=50
breaker=100
seamBreaker=10
breakCost=9999
numberOfQubits=512
qubitsPerCell=8
halfQPC=qubitsPerCell/2
gridSize=int(math.sqrt(numberOfQubits/qubitsPerCell))
rowDifference=qubitsPerCell*gridSize
cellLimit=2*gridSize
couplers='couplers'+str(gridSize)+'.txt'



def setNumQubits(noq):
    global numberOfQubits, gridSize, rowDifference, cellLimit, couplers
    numberOfQubits = noq
    gridSize=int(math.sqrt(numberOfQubits/qubitsPerCell))
    rowDifference=qubitsPerCell*gridSize
    cellLimit=2*gridSize
    couplers='couplers'+str(gridSize)+'.txt'
    return

def placeQubits(adjacency,order,numAdj,connect):
    rows=defineRows()
    columns=defineColumns()
    index=[order[0]] #Cells whose adjacent cells are to be placed
    toDo = [] #Cells to be imported into index on the next loop iteration
    #toDoNext=[]
    placed = [0]*(len(adjacency)) #1 if cell has been placed, 0 if not
    placed[order[0]]=1
    allPlaced = [1]*len(adjacency) #Used as value to compare placed against
    qsMod=int((gridSize-1)/2) #Ensures that qubitStart is placed on the middlemost qubit
    qubitStart = qubitsPerCell*qsMod*(gridSize+1) #Qubit that cell order[0] is to be placed on
    takenQubits =[0]*(numberOfQubits) #0 if qubit is available, 1 if not
    takenQubits[qubitStart]=1 #Because cell 0 will occupy qubit 'qubitStart'
    allTaken=[1]*(numberOfQubits) #Used as value to compare takenQubits against
    qubits=[-1]*len(adjacency) #The nth entry is the qubit that cell n lies on
    qubits[order[0]]=qubitStart
    full=False
    disable=[qubitStart]
    adjacentCells=list(numAdj)
    flag=[0]*len(adjacency)
    flaggedQubits=[[] for i in range(len(adjacency))]
    toConnect=[[] for i in range(len(connect))]
    for i in range(len(connect)):
        toConnect[i]=list(connect[i])
    placedAdjacent=[[] for i in range(len(adjacency))]
    for k in range(len(toConnect[order[0]])):
        if order[0] not in placedAdjacent[toConnect[order[0]][k]]:
            placedAdjacent[toConnect[order[0]][k]].append(order[0])
    outsideCouplers=[]

    route=[]
    routefile=open('routes.txt','w')
    trace=open('trace.txt','w')
    #solution=open('solution.txt','w')
    #solution.write('new:\n')

    while placed != allPlaced: #For as long as every cell needs placement
        #If there are no remaining qubits
        #if takenQubits==allTaken:
        #    full=True
        #    return [qubits,full]
        while index != []: #For as long as

            #Need to keep track of the next list of cells to place
            #For every cell in 'index', we check that...
            ACCopy=list(adjacentCells)
            for i in range(len(index)):

                for j in range(len(adjacency[order.index(index[i])])):
                    currentCell=adjacency[order.index(index[i])][j]
                    #print 'Placing cell',currentCell
                    #If the cell to be placed is the last cell adjacent to cell 'index[i]', then it can occupy the last qubit
                    #That is adjacent to the qubit occupied by 'index[i]'
                    #T/F if it is/is not the last item in the adjacency list

                    #That has not been placed...
                    if placed[currentCell] != 1:
                        #Call placing function
                        if len(placedAdjacent[currentCell]) == 1:
                            #LOWER FLAGS
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(index[i],flag,flaggedQubits,takenQubits,disable)
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(currentCell,flag,flaggedQubits,takenQubits,disable)

                            if index[i] in toConnect[currentCell] and currentCell in toConnect[index[i]]:
                                toConnect[index[i]].pop(toConnect[index[i]].index(currentCell))
                                toConnect[currentCell].pop(toConnect[currentCell].index(index[i]))
                                adjacentCells[index[i]]-=1
                                adjacentCells[currentCell]-=1
                                #placedAdjacent[index[i]].append(currentCell)
                                #placedAdjacent[currentCell].append(index[i])


                            #Find the shortest path to a qubit with sufficient adjacent qubits
                            adjQub=[[] for k in range(maxPaths)]

                            #adjQubOld=denseAdjacentQubits(qubits[index[i]],takenQubits)
                            adjQubOld=orderQubits(innerAdjacentQubits(qubits[index[i]],takenQubits),takenQubits)+orderQubits(outerAdjacentQubits(qubits[index[i]],takenQubits),takenQubits)

                            for k in range(len(adjQubOld)):
                                adjQub[k]=list([adjQubOld[k]])

                            noPlace=True
                            count=0
                            while noPlace:
                                count+=1
                                if count>=breaker:
                                    count=0
                                    done=False
                                    vacancy=perimeterVacancy(takenQubits,rows,columns)
                                    if vacancy==[0,0,0,0]:
                                        return [[],takenQubits,rows,columns,route]
                                    seams=availableSeams(qubits[index[i]],takenQubits,outsideCouplers)
                                    noSeam=True
                                    if seams != []:
                                        for k in range(len(seams)):
                                            if (seams[k] <= gridSize and (vacancy[1]+vacancy[2]) > 0) or (seams[k] > gridSize and (vacancy[0]+vacancy[3]) > 0):
                                                noSeam=False
                                                #print '\n\n\n\n\n4\n\n\n\n\n'
                                                [[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers,toAppend]=openSeam(seams[k],vacancy,rows,columns,[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers)
                                                if toAppend==-1:
                                                    return [[],takenQubits,rows,columns,route]
                                                for l in range(len(toAppend)):
                                                    takenQubits[toAppend[l]]=1
                                                    disable.append(toAppend[l])

                                                adjQub=[[] for l in range(maxPaths)]
                                                adjQubOld=orderQubits(innerAdjacentQubits(qubits[index[i]],takenQubits),takenQubits)+orderQubits(outerAdjacentQubits(qubits[index[i]],takenQubits),takenQubits)
                                                for l in range(len(adjQubOld)):
                                                    adjQub[l]=list([adjQubOld[l]])
                                                break
                                    if noSeam:
                                        return [[],takenQubits,rows,columns,route]

                                for k in range(len(adjQub)):
                                    if len(adjQub[k])==0:
                                        break
                                    for l in range(len(adjQub[k])-1):
                                        takenQubits[adjQub[k][l]]=1
                                    #If this qubit is sufficient
                                    if numAdjacency(currentCell,adjQub[k][-1],takenQubits,qubits,toConnect) >= adjacentCells[currentCell]:
                                        #print 'Place on qubit '+str(adjQub[k][-1])
                                        #print 'It has '+str(numAdjacency(currentCell,adjQub[k][-1],takenQubits,qubits,toConnect))+' adjacent qubits'
                                        qubits[currentCell]=adjQub[k][-1]
                                        if adjQub[k][-1] in flaggedQubits[index[i]]:
                                            flaggedQubits[index[i]].pop(flaggedQubits[index[i]].index(adjQub[k][-1]))
                                        if adjQub[k][-1] in flaggedQubits[currentCell]:
                                            flaggedQubits[currentCell].pop(flaggedQubits[currentCell].index(adjQub[k][l]))
                                        takenQubits[adjQub[k][-1]]=1
                                        disable.append(adjQub[k][-1])
                                        for l in range(len(adjQub[k])-1):
                                            takenQubits[adjQub[k][l]]=0
                                        noPlace=False
                                        break
                                    else:
                                        #print 'qubit '+str(adjQub[k][-1])+' only has '+str(numAdjacency(currentCell,adjQub[k][-1],takenQubits,qubits,toConnect))+' adjacent qubits'
                                        for l in range(len(adjQub[k])-1):
                                            takenQubits[adjQub[k][l]]=0
                                if not noPlace:
                                    break
                                #Copy adjQub to adjQubOld
                                adjQubOld=[]
                                for k in range(len(adjQub)):
                                    if adjQub[k] != []:
                                        adjQubOld.append(adjQub[k])
                                    else:
                                        break

                                #Erase adjQub
                                for k in range(len(adjQub)):
                                    if adjQub[k]!=[]:
                                        adjQub[k]=[]
                                    else:
                                        break

                                #For each route...
                                for k in range(len(adjQubOld)):
                                    #Mark the leading qubits as taken...
                                    for l in range(len(adjQubOld[k])-1):
                                        takenQubits[adjQubOld[k][l]]=1

                                    #And find the qubits adjacent to the last on this list.
                                    #tempAdjQub=denseAdjacentQubits(adjQubOld[k][-1],takenQubits)
                                    tempAdjQub=orderQubits(innerAdjacentQubits(adjQubOld[k][-1],takenQubits),takenQubits)+orderQubits(outerAdjacentQubits(adjQubOld[k][-1],takenQubits),takenQubits)

                                    #Unmark the leading qubits
                                    for l in range(len(adjQubOld[k])-1):
                                        takenQubits[adjQubOld[k][l]]=0
                                    #Place this new list in the first empty spot in adjQub
                                    for l in range(len(tempAdjQub)):
                                        for m in range(len(adjQub)):
                                            if adjQub[m]==[]:
                                                adjQub[m]=list(adjQubOld[k])
                                                adjQub[m].append(tempAdjQub[l])
                                                break

                            #print "Routing "+str(qubits[index[i]])+" to "+str(qubits[currentCell])
                            RoutingFunctionv6.disableQubits(disable)
                            cost=RoutingFunctionv6.Routing([[qubits[index[i]],qubits[currentCell]]],couplers,routefile)
                            if cost == breakCost:
                                return [[],takenQubits,rows,columns,route]
                            a= RoutingFunctionv6.getRoute()
                            #if a[0]==breakCOST:
                            for k in range(len(a[0][:-1])):
                                if abs(a[0][k+1]-a[0][k])==rowDifference or abs(a[0][k+1]-a[0][k])==qubitsPerCell:
                                    outsideCouplers.append([min(a[0][k+1],a[0][k]),max(a[0][k+1],a[0][k])])
                            route.append(a[0])
                            for k in range(1,len(a[0])-1):
                                if a[0][k] in flaggedQubits[index[i]]:
                                    flaggedQubits[index[i]].pop(flaggedQubits[index[i]].index(a[0][k]))
                                if a[0][k] in flaggedQubits[currentCell]:
                                    flaggedQubits[currentCell].pop(flaggedQubits[currentCell].index(a[0][k]))
                                takenQubits[a[0][k]]=1
                                if a[0][k] not in disable:
                                    disable.append(a[0][k])

                            #trace.write('Placing cell '+str(currentCell)+' adjacent to cell '+str(index[i])+'\n')
                            #trace.write('adjacentCells= '+str(adjacentCells)+'\n')
                            trace.write('Placing cell '+str(currentCell)+' adjacent to cell '+str(index[i])+'\r\n')
                            trace.write('adjacentCells= '+str(adjacentCells)+'\r\n')

                            #Remembers the just placed cell for 'index' on the next iteration and notes that it has been placed
                            toDo.append(currentCell)
                            placed[currentCell]=1
                            for k in range(len(toConnect[currentCell])):
                                if currentCell not in placedAdjacent[toConnect[currentCell][k]]:
                                    placedAdjacent[toConnect[currentCell][k]].append(currentCell)


                        elif len(placedAdjacent[currentCell]) == 2:
                            #LOWER FLAGS
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(placedAdjacent[currentCell][0],flag,flaggedQubits,takenQubits,disable)
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(placedAdjacent[currentCell][1],flag,flaggedQubits,takenQubits,disable)

                            for k in range(2):
                                if placedAdjacent[currentCell][k] in toConnect[currentCell] and currentCell in toConnect[placedAdjacent[currentCell][k]]:
                                    toConnect[placedAdjacent[currentCell][k]].pop(toConnect[placedAdjacent[currentCell][k]].index(currentCell))
                                    toConnect[currentCell].pop(toConnect[currentCell].index(placedAdjacent[currentCell][k]))
                                    adjacentCells[placedAdjacent[currentCell][k]]-=1
                                    adjacentCells[currentCell]-=1
                                    #placedAdjacent[placedAdjacent[currentCell][k]].append(currentCell)
                                    #placedAdjacent[currentCell].append(placedAdjacent[currentCell][k])
                            tempDisable=[]
                            while True:

                                #route placedAdjacent[currentCell][0] to index[i]
                                RoutingFunctionv6.disableQubits(disable+tempDisable)
                                cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][0]],qubits[placedAdjacent[currentCell][1]]]],couplers,routefile)

                                noPath=True
                                while noPath:
                                    while cost == breakCost:
                                        vacancy=perimeterVacancy(takenQubits,rows,columns)
                                        seam=findSeam(qubits[placedAdjacent[currentCell][0]],qubits[placedAdjacent[currentCell][1]],takenQubits,outsideCouplers,vacancy)
                                        if seam == -1:
                                            return [[],takenQubits,rows,columns,route]
                                        #print '\n\n\n\n\n3\n\n\n\n\n'
                                        [[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers,toAppend]=openSeam(seam,vacancy,rows,columns,[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers)
                                        if toAppend==-1:
                                            return [[],takenQubits,rows,columns,route]
                                        for l in range(len(toAppend)):
                                            takenQubits[toAppend[l]]=1
                                            disable.append(toAppend[l])
                                        RoutingFunctionv6.disableQubits(disable+tempDisable)
                                        RoutingFunctionv6.resetCouplers(couplers)
                                        cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][0]],qubits[placedAdjacent[currentCell][1]]]],couplers,routefile)

                                    a= RoutingFunctionv6.getRoute()
                                    noPath=False
                                    if len(a[0])==2:
                                        noPath=True
                                        RoutingFunctionv6.disableCouplers([[min(a[0]),max(a[0])]],couplers)
                                        cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][0]],qubits[placedAdjacent[currentCell][1]]]],couplers,routefile)

                                if len(a[0])==3:
                                    #IF SUFFICIENT
                                    if numAdjacency(currentCell,a[0][1],takenQubits,qubits,toConnect) >= adjacentCells[currentCell]:
                                        route.append([a[0][0],a[0][1]])
                                        route.append([a[0][1],a[0][2]])
                                        for k in range(len(a[0][:-1])):
                                            if abs(a[0][k+1]-a[0][k])==rowDifference or abs(a[0][k+1]-a[0][k])==qubitsPerCell:
                                                outsideCouplers.append([min(a[0][k+1],a[0][k]),max(a[0][k+1],a[0][k])])
                                        if a[0][1] in flaggedQubits[placedAdjacent[currentCell][0]]:
                                            flaggedQubits[placedAdjacent[currentCell][0]].pop(flaggedQubits[placedAdjacent[currentCell][0]].index(a[0][1]))
                                        if a[0][1] in flaggedQubits[placedAdjacent[currentCell][1]]:
                                            flaggedQubits[placedAdjacent[currentCell][1]].pop(flaggedQubits[placedAdjacent[currentCell][1]].index(a[0][1]))
                                        if a[0][1] not in disable:
                                            disable.append(a[0][1])
                                        else:
                                            #print 'This should never happen...'
                                            routefile.close()
                                            return [[],takenQubits,rows,columns,route]

                                        takenQubits[a[0][1]]=1
                                        disable.append(a[0][1])
                                        qubits[currentCell]=a[0][1]
                                        placed[currentCell]=1
                                        for k in range(len(toConnect[currentCell])):
                                            if currentCell not in placedAdjacent[toConnect[currentCell][k]]:
                                                placedAdjacent[toConnect[currentCell][k]].append(currentCell)
                                        toDo.append(currentCell)
                                        break
                                    else:
                                        tempDisable.append(a[0][1])
                                else:
                                    #Find the most sufficient for now
                                    for k in range(1,len(a[0])-1):
                                        takenQubits[a[0][k]]=1

                                    sufficiency=[]
                                    for k in range(1,len(a[0])-1):
                                        sufficiency.append(numAdjacency(currentCell,a[0][k],takenQubits,qubits,toConnect))
                                    maxSufficiency=max(sufficiency)
                                    if maxSufficiency >= adjacentCells[currentCell]:
                                        selectIndex=sufficiency.index(maxSufficiency)+1
                                        selectQubit=a[0][selectIndex]
                                        #ADD THESE ROUTES FOR routes.txt
                                        route.append(a[0][:selectIndex+1])
                                        route.append(a[0][selectIndex:])
                                        for k in range(len(a[0][:-1])):
                                            if abs(a[0][k+1]-a[0][k])==rowDifference or abs(a[0][k+1]-a[0][k])==qubitsPerCell:
                                                outsideCouplers.append([min(a[0][k+1],a[0][k]),max(a[0][k+1],a[0][k])])
                                        for k in range(1,len(a[0][:-1])):
                                            for l in range(len(placedAdjacent[currentCell])):
                                                if a[0][k] in flaggedQubits[placedAdjacent[currentCell][l]]:
                                                    flaggedQubits[placedAdjacent[currentCell][l]].pop(flaggedQubits[placedAdjacent[currentCell][l]].index(a[0][k]))
                                            if a[0][k] in flaggedQubits[index[i]]:
                                                flaggedQubits[index[i]].pop(flaggedQubits[index[i]].index(a[0][k]))

                                        qubits[currentCell]=selectQubit
                                        placed[currentCell]=1
                                        for k in range(len(toConnect[currentCell])):
                                            if currentCell not in placedAdjacent[toConnect[currentCell][k]]:
                                                placedAdjacent[toConnect[currentCell][k]].append(currentCell)
                                        toDo.append(currentCell)
                                        for k in range(1,len(a[0])-1):
                                            disable.append(a[0][k])
                                        break
                                    else:
                                        for k in range(1,len(a[0])-1):
                                            takenQubits[a[0][k]]=0
                                            tempDisable.append(a[0][k])

                                    #else: MODIFY PATH!!!!

                        else:
                            #LOWER FLAGS
                            length=len(placedAdjacent[currentCell])
                            for k in range(length):
                                flag,flaggedQubits,takenQubits,disable=lowerFlags(placedAdjacent[currentCell][k],flag,flaggedQubits,takenQubits,disable)

                            for k in range(length):
                                if placedAdjacent[currentCell][k] in toConnect[currentCell] and currentCell in toConnect[placedAdjacent[currentCell][k]]:
                                    toConnect[placedAdjacent[currentCell][k]].pop(toConnect[placedAdjacent[currentCell][k]].index(currentCell))
                                    toConnect[currentCell].pop(toConnect[currentCell].index(placedAdjacent[currentCell][k]))
                                    adjacentCells[placedAdjacent[currentCell][k]]-=1
                                    adjacentCells[currentCell]-=1
                                    #placedAdjacent[placedAdjacent[currentCell][k]].append(currentCell)
                                    #placedAdjacent[currentCell].append(placedAdjacent[currentCell][k])

                            noPath=True
                            #Length found using triangular numbers: like factorial but with addition
                            triLength=length*(length-1)/2
                            tempPath=[[] for k in range(triLength)]
                            tempPathLength=[[] for k in range(triLength)]
                            sufficient=[[] for k in range(triLength)]
                            tempDisable=[[] for k in range(triLength)]
                            potentialPaths=[[] for k in range(triLength)]
                            potentialPathsLength=[[] for k in range(triLength)]

                            while noPath:
                                noPath=False
                                kIndex=-length
                                for k in range(length-1):
                                    kIndex+=(length-k)
                                    for l in range(k+1,length):
                                        lIndex=l-1-k
                                        varIndex=kIndex+lIndex
                                        #Route each pair
                                        RoutingFunctionv6.disableQubits(disable)
                                        cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][k]],qubits[placedAdjacent[currentCell][l]]]],couplers,routefile)

                                        #Route each sufficient qubit on that path to remaining qubits
                                        done=False
                                        while not done:
                                            while cost == breakCost:
                                                tempPath[varIndex]=[1]*breakCost
                                                done=False
                                                vacancy=perimeterVacancy(takenQubits,rows,columns)
                                                seam=findSeam(qubits[placedAdjacent[currentCell][k]],qubits[placedAdjacent[currentCell][l]],takenQubits,outsideCouplers,vacancy)
                                                if seam == -1:
                                                    return [[],takenQubits,rows,columns,route]
                                                #print '\n\n\n\n\n1\n\n\n\n\n'

                                                [[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers,toAppend]=openSeam(seam,vacancy,rows,columns,[disable,qubits,flaggedQubits],[route],[takenQubits],outsideCouplers)
                                                if toAppend==-1:
                                                    return [[],takenQubits,rows,columns,route]
                                                for m in range(len(toAppend)):
                                                    takenQubits[toAppend[m]]=1
                                                    disable.append(toAppend[m])
                                                RoutingFunctionv6.disableQubits(disable)
                                                RoutingFunctionv6.resetCouplers(couplers)
                                                cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][k]],qubits[placedAdjacent[currentCell][l]]]],couplers,routefile)

                                            tempPath[varIndex]= RoutingFunctionv6.getRoute()
                                            done=True
                                            if len(tempPath[varIndex][0])==2:
                                                RoutingFunctionv6.disableCouplers([[min(tempPath[varIndex][0]),max(tempPath[varIndex][0])]],couplers)
                                                cost=RoutingFunctionv6.Routing([[qubits[placedAdjacent[currentCell][k]],qubits[placedAdjacent[currentCell][l]]]],couplers,routefile)
                                                done=False

                                        for m in range(1,len(tempPath[varIndex][0])-1):
                                            tempDisable[varIndex].append(tempPath[varIndex][0][m])
                                            takenQubits[tempPath[varIndex][0][m]]=1

                                        for m in range(len(tempDisable[varIndex])):
                                            if numAdjacency(currentCell,tempPath[varIndex][0][m],takenQubits,qubits,toConnect) >= adjacentCells[currentCell]:
                                                sufficient[varIndex].append(1)
                                            else:
                                                sufficient[varIndex].append(0)

                                        RoutingFunctionv6.disableQubits(tempDisable[varIndex]+disable)
                                        for m in range(len(tempDisable[varIndex])):
                                            if sufficient[varIndex][m] == 1:
                                                tempPP=[]
                                                for n in range(length):
                                                    if n!=k and n!=l:
                                                        cost=RoutingFunctionv6.Routing([[tempPath[varIndex][0][m+1],qubits[placedAdjacent[currentCell][n]]]],couplers,routefile)
                                                        if cost != breakCost:
                                                            b=RoutingFunctionv6.getRoute()
                                                            tempPP.append(b[0])
                                                        else:
                                                            tempPP.append([])
                                                potentialPaths[varIndex].append(tempPP)
                                            else:
                                                potentialPaths[varIndex]=[[]]


                                        for m in range(len(tempDisable[varIndex])):
                                            takenQubits[tempDisable[varIndex][m]]=0

                                    if noPath:
                                        break

                                if not noPath:
                                    for k in range(len(tempPath)):
                                        tempPathLength[k]=len(tempPath[k][0])-2

                                    for k in range(len(tempPath)):
                                        potentialPathsLength[k]=[[] for l in range(len(tempPath[k][0])-2)]

                                    for k in range(len(potentialPaths)):
                                        for l in range(len(potentialPaths[k])):
                                            for m in range(len(potentialPaths[k][l])):
                                                if len(potentialPaths[k][l][m]) == 0:
                                                    potentialPathsLength[k][l].append(breakCost)
                                                else:
                                                    potentialPathsLength[k][l].append(len(potentialPaths[k][l][m])-2)

                                    potentialPathsLengthTotal=[[] for k in range(len(potentialPathsLength))]
                                    for k in range(len(potentialPathsLength)):
                                        for l in range(len(potentialPathsLength[k])):
                                            potentialPathsLengthTotal[k].append(sum(potentialPathsLength[k][l]))

                                    minPPL=[]
                                    for k in range(len(potentialPathsLengthTotal)):
                                        minPPL.append(min(potentialPathsLengthTotal[k]))
                                    total=[]
                                    for k in range(len(tempPathLength)):
                                        total.append(tempPathLength[k]+minPPL[k])
                                    if min(total)>=breakCost:
                                        noPath=True
                                        #return [[],takenQubits,rows,columns,route]
                                    tempPathIndex=total.index(min(total))
                                    potentialPathsIndex=potentialPathsLengthTotal[tempPathIndex].index(min(potentialPathsLengthTotal[tempPathIndex]))

                                    qubits[currentCell]=tempPath[tempPathIndex][0][potentialPathsIndex+1]
                                    route.append(tempPath[tempPathIndex][0][:potentialPathsIndex+2])
                                    route.append(tempPath[tempPathIndex][0][potentialPathsIndex+1:])

                                    for k in range(len(potentialPaths[tempPathIndex][potentialPathsIndex])):
                                        route.append(potentialPaths[tempPathIndex][potentialPathsIndex][k])

                                    for m in range(len(tempPath[tempPathIndex][0][:-1])):
                                        if abs(tempPath[tempPathIndex][0][m+1]-tempPath[tempPathIndex][0][m])==rowDifference or abs(tempPath[tempPathIndex][0][m+1]-tempPath[tempPathIndex][0][m])==qubitsPerCell:
                                            outsideCouplers.append([min(tempPath[tempPathIndex][0][m+1],tempPath[tempPathIndex][0][m]),max(tempPath[tempPathIndex][0][m+1],tempPath[tempPathIndex][0][m])])

                                    for m in range(len(potentialPaths[tempPathIndex][potentialPathsIndex])):
                                        for n in range(len(potentialPaths[tempPathIndex][potentialPathsIndex][m][:-1])):
                                            if abs(potentialPaths[tempPathIndex][potentialPathsIndex][m][n+1]-potentialPaths[tempPathIndex][potentialPathsIndex][m][n])==rowDifference or abs(potentialPaths[tempPathIndex][potentialPathsIndex][m][n+1]-potentialPaths[tempPathIndex][potentialPathsIndex][m][n])==qubitsPerCell:
                                                outsideCouplers.append([min(potentialPaths[tempPathIndex][potentialPathsIndex][m][n+1],potentialPaths[tempPathIndex][potentialPathsIndex][m][n]),max(potentialPaths[tempPathIndex][potentialPathsIndex][m][n+1],potentialPaths[tempPathIndex][potentialPathsIndex][m][n])])

                                    for m in range(1,len(tempPath[tempPathIndex][0])-1):
                                        takenQubits[tempPath[tempPathIndex][0][m]]=1
                                        disable.append(tempPath[tempPathIndex][0][m])
                                        for n in range(len(placedAdjacent[currentCell])):
                                            if tempPath[tempPathIndex][0][m] in flaggedQubits[placedAdjacent[currentCell][n]]:
                                                flaggedQubits[placedAdjacent[currentCell][n]].pop(flaggedQubits[placedAdjacent[currentCell][n]].index(tempPath[tempPathIndex][0][m]))
                                        if tempPath[tempPathIndex][0][m] in flaggedQubits[index[i]]:
                                            flaggedQubits[index[i]].pop(flaggedQubits[index[i]].index(tempPath[tempPathIndex][0][m]))

                                    for m in range(len(potentialPaths[tempPathIndex][potentialPathsIndex])):
                                        for n in range(1,len(potentialPaths[tempPathIndex][potentialPathsIndex][m])-1):
                                            takenQubits[potentialPaths[tempPathIndex][potentialPathsIndex][m][n]]=1
                                            disable.append(potentialPaths[tempPathIndex][potentialPathsIndex][m][n])
                                            for o in range(len(placedAdjacent[currentCell])):
                                                if potentialPaths[tempPathIndex][potentialPathsIndex][m][n] in flaggedQubits[placedAdjacent[currentCell][o]]:
                                                    flaggedQubits[placedAdjacent[currentCell][o]].pop(flaggedQubits[placedAdjacent[currentCell][o]].index(potentialPaths[tempPathIndex][potentialPathsIndex][m][n]))
                                            if potentialPaths[tempPathIndex][potentialPathsIndex][m][n] in flaggedQubits[index[i]]:
                                                flaggedQubits[index[i]].pop(flaggedQubits[index[i]].index(potentialPaths[tempPathIndex][potentialPathsIndex][m][n]))

                                    placed[currentCell]=1

                                    for m in range(len(toConnect[currentCell])):
                                        if currentCell not in placedAdjacent[toConnect[currentCell][m]]:
                                            placedAdjacent[toConnect[currentCell][m]].append(currentCell)

                                    toDo.append(currentCell)


                    #RAISE FLAGS AGAIN
                    flag,flaggedQubits,takenQubits,disable=raiseFlags(flag,flaggedQubits,takenQubits,disable,placed,qubits,order,toConnect,adjacentCells)

                    #trace.write('qubits= '+str(qubits)+'\n')
                    trace.write('qubits= '+str(qubits)+'\r\n')
                    t=[]
                    for k in range(len(takenQubits)):
                        if takenQubits[k]==1:
                            t.append(k)
                    #trace.write('TakenQubits= '+str(t)+'\n')
                    #trace.write('Disabled= '+str(disable)+'\n')
                    #trace.write('toConnect= '+str(toConnect)+'\n')
                    trace.write('TakenQubits= '+str(t)+'\r\n')
                    trace.write('Disabled= '+str(disable)+'\r\n')
                    trace.write('toConnect= '+str(toConnect)+'\r\n')
                    trace.write('Flagged Qubits: ')
                    for k in range(len(flaggedQubits)):
                        if flaggedQubits[k] != []:
                            #trace.write('Cell '+str(k)+': '+str(flaggedQubits[k])+'\n')
                    #trace.write('\n\n')
                            trace.write('Cell '+str(k)+': '+str(flaggedQubits[k])+'\r\n')
                    trace.write('\r\n\r\n')

                    #Route to everything that has been placed
                    toRemove1=[]
                    toRemove2=[]
                    for k in range(len(toConnect[currentCell])):
                        if placed[toConnect[currentCell][k]]==1:
                            #print 'This should never happen'
                            return [[],takenQubits,rows,columns,route]

                            #LOWER FLAGS
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(currentCell,flag,flaggedQubits,takenQubits,disable)
                            flag,flaggedQubits,takenQubits,disable=lowerFlags(toConnect[currentCell][k],flag,flaggedQubits,takenQubits,disable)

                            RoutingFunctionv6.disableQubits(disable)
                            RoutingFunctionv6.Routing([[qubits[currentCell],qubits[toConnect[currentCell][k]]]],couplers,routefile)
                            a= RoutingFunctionv6.getRoute()
                            for l in range(len(a[0][:-1])):
                                if abs(a[0][l+1]-a[0][l])==rowDifference or abs(a[0][l+1]-a[0][l])==qubitsPerCell:
                                    outsideCouplers.append([min(a[0][l+1],a[0][l]),max(a[0][l+1],a[0][l])])
                            route.append(a[0])
                            for l in range(1,len(a[0])-1):
                                if a[0][l] in flaggedQubits[currentCell]:
                                    flaggedQubits[currentCell].pop(flaggedQubits[currentCell].index(a[0][l]))
                                if a[0][l] in flaggedQubits[toConnect[currentCell][k]]:
                                    flaggedQubits[toConnect[currentCell][k]].pop(flaggedQubits[toConnect[currentCell][k]].index(a[0][l]))
                                takenQubits[a[0][l]]=1
                                if a[0][l] not in disable:
                                    disable.append(a[0][l])

                            #trace.write('Routing cell '+str(toConnect[currentCell][k])+' to cell '+str(currentCell)+'\n')
                            #trace.write('adjacentCells= '+str(adjacentCells)+'\n')
                            #trace.write('qubits= '+str(qubits)+'\n')
                            #t=[]
                            #for l in range(len(takenQubits)):
                            #    if takenQubits[l]==1:
                            #        t.append(l)
                            #trace.write('TakenQubits= '+str(t)+'\n')
                            #trace.write('Disabled= '+str(disable)+'\n')
                            #trace.write('toConnect= '+str(toConnect)+'\n')
                            #trace.write('Flagged Qubits: ')
                            #for l in range(len(flaggedQubits)):
                            #    if flaggedQubits[l] != []:
                            #        trace.write('Cell '+str(l)+': '+str(flaggedQubits[l])+'\n')
                            #trace.write('\n\n')

                            adjacentCells[currentCell]-=1
                            adjacentCells[toConnect[currentCell][k]]-=1
                            #placedAdjacent[currentCell].append(toConnect[currentCell][k])
                            #placedAdjacent[toConnect[currentCell][k]].append(currentCell)
                            toRemove1.append(toConnect[currentCell][k])
                            toRemove2.append(currentCell)

                            #RAISE FLAGS AGAIN
                            flag,flaggedQubits,takenQubits,disable=raiseFlags(flag,flaggedQubits,takenQubits,disable,placed,qubits,order,toConnect,adjacentCells)

                    for k in range(len(toRemove1)):
                        toConnect[toRemove1[k]].pop(toConnect[toRemove1[k]].index(toRemove2[k]))
                        toConnect[toRemove2[k]].pop(toConnect[toRemove2[k]].index(toRemove1[k]))


            #Once every adjacent cell to every cell in 'index' has been placed, 'index' is given the cells who have just been placed
            #So that cells adjacent to those can be placed
            toDo=reorder(toDo,ACCopy)
            index=toDo
            toDo=[]

        #In the case that index = [] without every cell being placed, the cell with the lowest index number is appended to 'index'...
        for i in range(len(order)):
            if placed[order[i]]==0:
                index.append(order[i])
                #And marked as placed
                placed[order[i]]=1
                for k in range(len(toConnect[order[i]])):
                    if order[i] not in placedAdjacent[toConnect[order[i]][k]]:
                        placedAdjacent[toConnect[order[i]][k]].append(order[i])

                #Place this cell on the middlemost (y-axis only) available qubit with adequate adjacent qubits
                for j in range(numberOfQubits-qubitStart):
                    if takenQubits[qubitStart+j]==0:
                        if len(denseAdjacentQubits(qubitStart+j,takenQubits)) >= len(adjacency[order.index(order[i])]):
                            disable.append(qubitStart+j)
                            takenQubits[qubitStart+j]=1
                            qubits[order[i]]=qubitStart+j
                            break
                #If that doesn't work we'll try going down from qubitStart
                    if j < qubitStart:
                        if takenQubits[qubitStart-1-j]==0:
                            if len(denseAdjacentQubits(qubitStart-1-j,takenQubits)) >= len(adjacency[order.index(order[i])]):
                                disable.append(qubitStart-1-j)
                                takenQubits[qubitStart-1-j]=1
                                qubits[order[i]]=qubitStart-1-j
                                break

                #Error, no available qubits!
                #This break should really be an error message
                break
    #print route
    #print '\n\noutsideCouplers= '+str(outsideCouplers)
    #for i in range(len(route)):
    #    for j in range(len(route[i][:-1])):
    #        solution.write(str(route[i][j])+' ')
    #    solution.write(str(route[i][-1])+'\n')
    trace.close()
    routefile.close()
    #solution.write('shared:\n')
    #solution.close()
    #RoutingFunctionv6.setRoute(route,couplers)
    return [qubits,takenQubits,rows,columns,route]
    #return [qubits,disable]
    #return qubits

def reorderVer2(adjacency,connections):
    adjacentCells=[0]*len(adjacency)
    newAdjacency=[[]]*len(adjacency)
    order=[]
    toRemove=-1

    #Count the adjacencies per cell
    for i in range(len(adjacency)):
        for j in range(len(adjacency[i])):
            adjacentCells[i]+=1
            adjacentCells[adjacency[i][j]]+=1
    nACopy=list(adjacentCells)

    #Order ranks the cells from most adjacencies to least
    for i in range(len(adjacentCells)):
        order.append(adjacentCells.index(max(adjacentCells)))
        adjacentCells[adjacentCells.index(max(adjacentCells))]=-1

    #Create adjacency list that lists all adjacencies
    for i in range(len(order)):
        temp=[]
        for j in range(len(adjacency)):
            if j==order[i]:
                for k in range(len(adjacency[j])):
                    temp.append(adjacency[j][k])
            for k in range(len(adjacency[j])):
                if adjacency[j][k]==order[i]:
                    temp.append(j)

        newAdjacency[order[i]]=temp

    #Orders the adjacency list starting with most popular cells
    for i in range(len(newAdjacency)):
        tempOrder=[]
        tempNA=[]
        tempAdjacency=[]
        for j in range(len(newAdjacency[i])):
            tempNA.append(nACopy[newAdjacency[i][j]])

        for j in range(len(tempNA)):
            tempOrder.append(tempNA.index(max(tempNA)))
            tempNA[tempNA.index(max(tempNA))]=-1

        for j in range(len(tempNA)):
            tempAdjacency.append(newAdjacency[i][tempOrder[j]])
        newAdjacency[i]=tempAdjacency

    #Orders the adjacency list in the order that cells will be placed
    newOrder=[order[0]]
    index=[order[0]]
    toDo=[]
    currentAC=list(nACopy)
    placed=[0]*len(order)
    placed[order[0]]=1
    allPlaced=[1]*len(order)
    toConnect=[[] for i in range(len(connections))]
    for i in range(len(connections)):
        toConnect[i]=list(connections[i])
    while placed != allPlaced:
        while index != []:
            adjacentToRemove=[]
            for i in range(len(index)):
                #For every cell adjacent to cell 'index[i]'...
                cell_i = index[i]
                for j in range(len(newAdjacency[index[i]])):
					adjacentCell = newAdjacency[index[i]][j]
                    if placed[adjacentCell]==0:
                        placed[adjacentCell]=1
                        #currentAC[index[i]]-=1
                        adjacentToRemove.append(index[i])
                        #currentAC[adjacentCell]-=1
                        adjacentToRemove.append(adjacentCell)
                        toConnect[index[i]].pop(toConnect[index[i]].index(adjacentCell))
                        toConnect[adjacentCell].pop(toConnect[adjacentCell].index(index[i]))
                        #newOrder.append(adjacentCell)
                        toDo.append(adjacentCell)
                        toRemove1=[]
                        toRemove2=[]
                        for k in range(len(toConnect[adjacentCell])):
                            if placed[toConnect[adjacentCell][k]]==1:
                                #currentAC[adjacentCell]-=1
                                adjacentToRemove.append(adjacentCell)
                                #currentAC[toConnect[adjacentCell][k]]-=1
                                adjacentToRemove.append(toConnect[adjacentCell][k])
                                toRemove1.append(adjacentCell)
                                toRemove2.append(toConnect[adjacentCell][k])
                        for k in range(len(toRemove1)):
                            toConnect[toRemove1[k]].pop(toConnect[toRemove1[k]].index(toRemove2[k]))
                            toConnect[toRemove2[k]].pop(toConnect[toRemove2[k]].index(toRemove1[k]))

            toDo=reorder(toDo,currentAC)
            for i in range(len(adjacentToRemove)):
                currentAC[adjacentToRemove[i]]-=1
            for i in range(len(toDo)):
                newOrder.append(toDo[i])
            index=toDo
            toDo=[]
        for i in range(len(order)):
            if placed[order[i]]==0:
                #print 'Ran out of adjacency list!'
                index.append(order[i])
                #And marked as placed
                placed[order[i]]=1
                newOrder.append(order[i])
                break

    temp=[]
    for i in range(len(newOrder)):
        temp.append(newAdjacency[newOrder[i]])
    newAdjacency=temp

    #Reduces redundancies
    for i in range(len(newAdjacency)):
        for j in range(i,len(newAdjacency)):
            for k in range(len(newAdjacency[j])):
                if newAdjacency[j][k] == newOrder[i]:
                    toRemove=k
            if toRemove != -1:
                newAdjacency[j].pop(toRemove)
                toRemove=-1

    return [newAdjacency,newOrder,nACopy]

def reorder(toDo,adjacentCells):
    AC=[]
    newToDo=[]
    for i in range(len(toDo)):
        AC.append(adjacentCells[toDo[i]])
    for i in range(len(toDo)):
        m=max(AC)
        j=AC.index(m)
        newToDo.append(toDo[j])
        AC[j]=-1
    return newToDo

#Orders qubits from most available to least available
def orderQubits(qubitList,takenQubits):
    qubitAvailability=[]
    newQubitList=[]
    for i in range(len(qubitList)):
        qubitAvailability.append(availability(qubitList[i],takenQubits))
    for i in range(len(qubitAvailability)):
        maxIndex=qubitAvailability.index(max(qubitAvailability))
        newQubitList.append(qubitList[maxIndex])
        qubitAvailability[maxIndex]=-1
    return newQubitList

#Returns a list of adjacent qubits to qubitIndex, with the first entry having the most
#adjacent qubits, and ties favoring a dense placement.
def denseAdjacentQubits(qubitIndex,takenQubits):
    adjQub=[]
    orderedQubits=[]
    numAdj=[]

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)
        if random.randint(0,1)==0:
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
        else:
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)

    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)
        if random.randint(0,1)==0:
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
        else:
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)

    #for i in range(len(adjQub)):
    #    numAdj.append(numAdjacency(adjQub[i],takenQubits,qubits,toConnect))

    #for i in range(len(numAdj)):
    #    orderedQubits.append(adjQub[numAdj.index(max(numAdj))])
    #    numAdj[numAdj.index(max(numAdj))]=-1

    #Depending on if you want ordered or not
    #return orderedQubits
    return adjQub

#Returns a list of adjacent qubits to qubitIndex, with the first entry having the most
#adjacent qubits, and ties favoring a dense placement.
def sparseAdjacentQubits(qubitIndex,takenQubits):
    adjQub=[]
    orderedQubits=[]
    numAdj=[]

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        if random.randint(0,1)==0:
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
        else:
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)

    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        if random.randint(0,1)==0:
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
        else:
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)

    #for i in range(len(adjQub)):
    #    numAdj.append(numAdjacency(adjQub[i],takenQubits,qubits,toConnect))

    #for i in range(len(numAdj)):
    #    orderedQubits.append(adjQub[numAdj.index(max(numAdj))])
    #    numAdj[numAdj.index(max(numAdj))]=-1

    #Depending on if you want ordered or not
    #return orderedQubits
    return adjQub

def outerAdjacentQubits(qubitIndex,takenQubits):
    adjQub=[]

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        if random.randint(0,1)==0:
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
        else:
            #Try the qubit in the cell below
            if qubitGroup < numberOfQubits/halfQPC-cellLimit:
                if takenQubits[qubitIndex+rowDifference]==0:
                    adjQub.append(qubitIndex+rowDifference)
            #Try the qubit in the cell above
            if qubitGroup >=cellLimit:
                if takenQubits[qubitIndex-rowDifference]==0:
                    adjQub.append(qubitIndex-rowDifference)

    else: #Qubit runs horizontally
        if random.randint(0,1)==0:
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
        else:
            #Try the qubit in the cell to the right
            if (qubitGroup+1)%cellLimit != 0:
                if takenQubits[qubitIndex+qubitsPerCell]==0:
                    adjQub.append(qubitIndex+qubitsPerCell)
            #Try the qubit in the cell to the left
            if (qubitGroup-1)%cellLimit != 0:
                if takenQubits[qubitIndex-qubitsPerCell]==0:
                    adjQub.append(qubitIndex-qubitsPerCell)

    return adjQub

def innerAdjacentQubits(qubitIndex,takenQubits):
    adjQub=[]

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)

    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                adjQub.append(qubitBase+i)

    return adjQub

def numAdjacency(cellIndex,qubitIndex,takenQubits,qubits,toConnect):
    number=0

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                number+=1
            else:
                if qubitBase+i in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitBase+i:
                            number+=1
        #Try the qubit in the cell above
        if qubitGroup >=cellLimit:
            if takenQubits[qubitIndex-rowDifference]==0:
                number+=1
            else:
                if qubitIndex-rowDifference in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex-rowDifference:
                            number+=1
        #Try the qubit in the cell below
        if qubitGroup < numberOfQubits/halfQPC-cellLimit:
            if takenQubits[qubitIndex+rowDifference]==0:
                number+=1
            else:
                if qubitIndex+rowDifference in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex+rowDifference:
                            number+=1
    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                number+=1
            else:
                if qubitBase+i in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitBase+i:
                            number+=1
        #Try the qubit in the cell to the left
        if (qubitGroup-1)%cellLimit != 0:
            if takenQubits[qubitIndex-qubitsPerCell]==0:
                number+=1
            else:
                if qubitIndex-qubitsPerCell in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex-qubitsPerCell:
                            number+=1
        #Try the qubit in the cell to the right
        if (qubitGroup+1)%cellLimit != 0:
            if takenQubits[qubitIndex+qubitsPerCell]==0:
                number+=1
            else:
                if qubitIndex+qubitsPerCell in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex+qubitsPerCell:
                            number+=1
    #Can only arrive here if there is one or less adjacent cells
    return number

def takenAdjacency(cellIndex,qubitIndex,takenQubits,qubits,toConnect):
    number=0

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==1:
                if qubitBase+i in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitBase+i:
                            number+=1
        #Try the qubit in the cell above
        if qubitGroup >=cellLimit:
            if takenQubits[qubitIndex-rowDifference]==1:
                if qubitIndex-rowDifference in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex-rowDifference:
                            number+=1
        #Try the qubit in the cell below
        if qubitGroup < numberOfQubits/halfQPC-cellLimit:
            if takenQubits[qubitIndex+rowDifference]==1:
                if qubitIndex+rowDifference in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex+rowDifference:
                            number+=1
    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==1:
                if qubitBase+i in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitBase+i:
                            number+=1
        #Try the qubit in the cell to the left
        if (qubitGroup-1)%cellLimit != 0:
            if takenQubits[qubitIndex-qubitsPerCell]==1:
                if qubitIndex-qubitsPerCell in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex-qubitsPerCell:
                            number+=1
        #Try the qubit in the cell to the right
        if (qubitGroup+1)%cellLimit != 0:
            if takenQubits[qubitIndex+qubitsPerCell]==1:
                if qubitIndex+qubitsPerCell in qubits:
                    for j in range(len(toConnect[cellIndex])):
                        if qubits[toConnect[cellIndex][j]] == qubitIndex+qubitsPerCell:
                            number+=1

    return number

def availability(qubitIndex,takenQubits):
    number=0

    qubitGroup=int(qubitIndex/halfQPC)
    if qubitGroup%2 == 0: #Qubit runs vertically
        qubitBase=(qubitGroup+1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                number+=1
        #Try the qubit in the cell above
        if qubitGroup >=cellLimit:
            if takenQubits[qubitIndex-rowDifference]==0:
                number+=1
        #Try the qubit in the cell below
        if qubitGroup < numberOfQubits/halfQPC-cellLimit:
            if takenQubits[qubitIndex+rowDifference]==0:
                number+=1
    else: #Qubit runs horizontally
        qubitBase=(qubitGroup-1)*halfQPC #The lowest horizontal qubit number within the qubit cell
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                number+=1
        #Try the qubit in the cell to the left
        if (qubitGroup-1)%cellLimit != 0:
            if takenQubits[qubitIndex-qubitsPerCell]==0:
                number+=1
        #Try the qubit in the cell to the right
        if (qubitGroup+1)%cellLimit != 0:
            if takenQubits[qubitIndex+qubitsPerCell]==0:
                number+=1
    #Can only arrive here if there is one or less adjacent cells
    return number

def connections(adjacency):
    connections=[[] for i in range(len(adjacency))]
    for i in range(len(adjacency)):
        for j in range(len(adjacency[i])):
            connections[i].append(adjacency[i][j])
            connections[adjacency[i][j]].append(i)
    return connections

def lowerFlags(cell,flag,flaggedQubits,takenQubits,disable):
    if flag[cell]==1:
        flag[cell]=0
        for k in range(len(flaggedQubits[cell])):
            takenQubits[flaggedQubits[cell][k]]=0
            while True:
                if flaggedQubits[cell][k] not in disable:
                    break
                disable.pop(disable.index(flaggedQubits[cell][k]))
    return [flag,flaggedQubits,takenQubits,disable]

def raiseFlags(flag,flaggedQubits,takenQubits,disable,placed,qubits,order,toConnect,adjacentCells):
    for k in range(len(placed)):
        for l in range(len(flaggedQubits[k])):
            takenQubits[flaggedQubits[k][l]]=0
            while True:
                if flaggedQubits[k][l] not in disable:
                    break
                disable.pop(disable.index(flaggedQubits[k][l]))
        flaggedQubits[k]=[]
        flag[k]=0

    for k in range(len(order)):
        if placed[order[k]]==1:
            tempFQ=denseAdjacentQubits(qubits[order[k]],takenQubits)
            if adjacentCells[order[k]]>=numAdjacency(order[k],qubits[order[k]],takenQubits,qubits,toConnect) and adjacentCells[order[k]] !=0:
                flag[order[k]]=1
                flaggedQubits[order[k]]=tempFQ
                for l in range(len(flaggedQubits[order[k]])):
                    takenQubits[flaggedQubits[order[k]][l]]=1
                    disable.append(flaggedQubits[order[k]][l])
    return [flag,flaggedQubits,takenQubits,disable]

def defineRows():
    rows=[[] for i in range(gridSize)]
    for i in range(gridSize):
        qubits=[]
        cellIndex=gridSize*qubitsPerCell*i
        for j in range(gridSize):
            for k in range(qubitsPerCell):
                qubits.append(cellIndex+k)
            cellIndex+=qubitsPerCell
        rows[i]=qubits
    return rows

def defineColumns():
    columns=[[] for i in range(gridSize)]
    for i in range(gridSize):
        qubits=[]
        cellIndex=qubitsPerCell*i
        for j in range(gridSize):
            for k in range(8):
                qubits.append(cellIndex+k)
            cellIndex+=gridSize*qubitsPerCell
        columns[i]=qubits
    return columns

def conduitPath(fromQubit,toQubit,takenQubits,outsideCouplers):
    fromQubitGroup=int(fromQubit/halfQPC)
    fromQubitCell=int(fromQubit/qubitsPerCell)
    if fromQubitGroup%2 == 0: #Qubit runs vertically
        print 'Qubit '+str(fromQubit)+' runs vertically'
        if fromQubitGroup >=cellLimit and fromQubitGroup < numberOfQubits/halfQPC-cellLimit:
            print 'Qubit '+str(fromQubit)+' is not in the top row or bottom row'
            if [fromQubit-rowDifference,fromQubit] not in outsideCouplers:
                print 'Seam '+str(int(fromQubitCell/gridSize)+gridSize+1)+' is available'
            if [fromQubit,fromQubit+rowDifference] not in outsideCouplers:
                print 'Seam '+str(int(fromQubitCell/gridSize)+gridSize+2)+' is available'
        else:
            if fromQubitGroup < cellLimit:
                print 'Qubit '+str(fromQubit)+' is in the top row'
                if [fromQubit,fromQubit+rowDifference] in outsideCouplers:
                    print 'Seam '+str(gridSize+1)+' is available'
                else:
                    print 'Seam '+str(gridSize+1)+' is available'
                    print 'Seam '+str(gridSize+2)+' is available'
            if qubitGroup >= numberOfQubits/halfQPC-cellLimit:
                print 'Qubit '+str(fromQubit)+' is in the bottow row'
                if [fromQubit-rowDifference,fromQubit] in outsideCouplers:
                    print 'Seam '+str(2*gridSize+1)+' is available'
                else:
                    print 'Seam '+str(2*gridSize)+' is available'
                    print 'Seam '+str(2*gridSize+1)+' is available'
        fromQubitBase=(fromQubitGroup+1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[fromQubitBase+i]==0:
                print 'Qubit '+str(fromQubit)+' can exit to the left or right'
                print 'Seam '+str(fromQubitCell%gridSize)+' is available'
                print 'Seam '+str(fromQubitCell%gridSize+1)+' is available'
                break
    else: #Qubit runs horizontally
        print 'Qubit '+str(fromQubit)+' runs horizontally'
        if (fromQubitGroup-1)%cellLimit != 0 and (fromQubitGroup+1)%cellLimit != 0:
            print 'Qubit '+str(fromQubit)+' is not in the left or right column'
            if [fromQubit-qubitsPerCell,fromQubit] not in outsideCouplers:
                print 'Seam '+str(fromQubitCell%gridSize)+' is available'
            if [fromQubit,fromQubit+qubitsPerCell] not in outsideCouplers:
                print 'Seam '+str(fromQubitCell%gridSize+1)+' is available'
        else:
            if (fromQubitGroup-1)%cellLimit == 0:
                print 'Qubit '+str(fromQubit)+' is in the left column'
                if [fromQubit,fromQubit+qubitsPerCell] in outsideCouplers:
                    print 'Seam 0 is available'
                else:
                    print 'Seam 0 is available'
                    print 'Seam 1 is available'
            if (fromQubitGroup+1)%cellLimit == 0:
                print 'Qubit '+str(fromQubit)+' is in the right column'
                if [fromQubit-qubitsPerCell,fromQubit] in outsideCouplers:
                    print 'Seam 8 is available'
                else:
                    print 'Seam 7 is available'
                    print 'Seam 8 is available'
        fromQubitBase=(fromQubitGroup-1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[fromQubitBase+i]==0:
                print 'Qubit '+str(fromQubit)+' can exit from the top or bottom'
                print 'Seam '+str(int(fromQubitCell/gridSize)+9)+' is available'
                print 'Seam '+str(int(fromQubitCell/gridSize)+10)+' is available'
                break

    toQubitGroup=int(toQubit/halfQPC)
    toQubitCell=int(toQubit/qubitsPerCell)
    if toQubitGroup%2 == 0: #Qubit runs vertically
        print 'Qubit '+str(toQubit)+' runs vertically'
        if toQubitGroup >=cellLimit and toQubitGroup < numberOfQubits/halfQPC-cellLimit:
            print 'Qubit '+str(toQubit)+' is not in the top row or bottom row'
            if [toQubit-rowDifference,toQubit] not in outsideCouplers:
                print 'Seam '+str(int(toQubitCell/gridSize)+gridSize+1)+' is available'
            if [toQubit,toQubit+rowDifference] not in outsideCouplers:
                print 'Seam '+str(int(toQubitCell/gridSize)+gridSize+2)+' is available'
        else:
            if toQubitGroup < cellLimit:
                print 'Qubit '+str(toQubit)+' is in the top row'
                if [toQubit,toQubit+rowDifference] in outsideCouplers:
                    print 'Seam '+str(gridSize+1)+' is available'
                else:
                    print 'Seam '+str(gridSize+1)+' is available'
                    print 'Seam '+str(gridSize+2)+' is available'
            if toQubitGroup >= numberOfQubits/halfQPC-cellLimit:
                print 'Qubit '+str(toQubit)+' is in the bottow row'
                if [toQubit-rowDifference,toQubit] in outsideCouplers:
                    print 'Seam '+str(2*gridSize+1)+' is available'
                else:
                    print 'Seam '+str(2*gridSize)+' is available'
                    print 'Seam '+str(2*gridSize+1)+' is available'
        toQubitBase=(toQubitGroup+1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[toQubitBase+i]==0:
                print 'Qubit '+str(toQubit)+' can exit to the left or right'
                print 'Seam '+str(toQubitCell%gridSize)+' is available'
                print 'Seam '+str(toQubitCell%gridSize+1)+' is available'
                break
    else: #Qubit runs horizontally
        print 'Qubit '+str(toQubit)+' runs horizontally'
        if (toQubitGroup-1)%cellLimit != 0 and (toQubitGroup+1)%cellLimit != 0:
            print 'Qubit '+str(toQubit)+' is not in the left or right column'
            if [toQubit-qubitsPerCell,toQubit] not in outsideCouplers:
                print 'Seam '+str(toQubitCell%gridSize)+' is available'
            if [toQubit,toQubit+qubitsPerCell] not in outsideCouplers:
                print 'Seam '+str(toQubitCell%gridSize+1)+' is available'
        else:
            if (toQubitGroup-1)%cellLimit == 0:
                print 'Qubit '+str(toQubit)+' is in the left column'
                if [toQubit,toQubit+qubitsPerCell] in outsideCouplers:
                    print 'Seam 0 is available'
                else:
                    print 'Seam 0 is available'
                    print 'Seam 1 is available'
            if (toQubitGroup+1)%cellLimit == 0:
                print 'Qubit '+str(toQubit)+' is in the right column'
                if [toQubit-qubitsPerCell,toQubit] in outsideCouplers:
                    print 'Seam '+str(gridSize)+' is available'
                else:
                    print 'Seam '+str(gridSize-1)+' is available'
                    print 'Seam '+str(gridSize)+' is available'
        toQubitBase=(toQubitGroup-1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[toQubitBase+i]==0:
                print 'Qubit '+str(toQubit)+' can exit from the top or bottom'
                print 'Seam '+str(int(toQubitCell/gridSize)+gridSize+1)+' is available'
                print 'Seam '+str(int(toQubitCell/gridSize)+gridSize+2)+' is available'
                break

    return

def availableSeams(qubit,takenQubits,outsideCouplers):
    seams=[]
    qubitGroup=int(qubit/halfQPC)
    qubitCell=int(qubit/qubitsPerCell)
    if qubitGroup%2 == 0: #Qubit runs vertically
        #print 'Qubit '+str(qubit)+' runs vertically'
        if qubitGroup >=cellLimit and qubitGroup < numberOfQubits/halfQPC-cellLimit:
            #print 'Qubit '+str(qubit)+' is not in the top row or bottom row'
            if [qubit-rowDifference,qubit] not in outsideCouplers:
                #print 'Seam '+str(int(qubitCell/gridSize)+gridSize+1)+' is available'
                seams.append(int(qubitCell/gridSize)+gridSize+1)
            if [qubit,qubit+rowDifference] not in outsideCouplers:
                #print 'Seam '+str(int(qubitCell/gridSize)+gridSize+2)+' is available'
                seams.append(int(qubitCell/gridSize)+gridSize+2)
        else:
            if qubitGroup < cellLimit:
                #print 'Qubit '+str(qubit)+' is in the top row'
                if [qubit,qubit+rowDifference] in outsideCouplers:
                    #print 'Seam '+str(gridSize+1)+' is available'
                    seams.append(gridSize+1)
                else:
                    #print'Seam '+str(gridSize+1)+' is available'
                    #print'Seam '+str(gridSize+2)+' is available'
                    seams.append(gridSize+1)
                    seams.append(gridSize+2)
            if qubitGroup >= numberOfQubits/halfQPC-cellLimit:
                #print'Qubit '+str(qubit)+' is in the bottow row'
                if [qubit-rowDifference,qubit] in outsideCouplers:
                    #print'Seam '+str(2*gridSize+1)+' is available'
                    seams.append(2*gridSize+1)
                else:
                    #print'Seam '+str(2*gridSize)+' is available'
                    seams.append(2*gridSize)
                    #print'Seam '+str(2*gridSize+1)+' is available'
                    seams.append(2*gridSize+1)
        qubitBase=(qubitGroup+1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                #print'Qubit '+str(qubit)+' can exit to the left or right'
                #print'Seam '+str(qubitCell%gridSize)+' is available'
                seams.append(qubitCell%gridSize)
                #print'Seam '+str(qubitCell%gridSize+1)+' is available'
                seams.append(qubitCell%gridSize+1)
                break
    else: #Qubit runs horizontally
        #print'Qubit '+str(qubit)+' runs horizontally'
        if (qubitGroup-1)%cellLimit != 0 and (qubitGroup+1)%cellLimit != 0:
            #print'Qubit '+str(qubit)+' is not in the left or right column'
            if [qubit-qubitsPerCell,qubit] not in outsideCouplers:
                #print'Seam '+str(qubitCell%gridSize)+' is available'
                seams.append(qubitCell%gridSize)
            if [qubit,qubit+qubitsPerCell] not in outsideCouplers:
                #print'Seam '+str(qubitCell%gridSize+1)+' is available'
                seams.append(qubitCell%gridSize+1)
        else:
            if (qubitGroup-1)%cellLimit == 0:
                #print'Qubit '+str(qubit)+' is in the left column'
                if [qubit,qubit+qubitsPerCell] in outsideCouplers:
                    #print'Seam 0 is available'
                    seams.append(0)
                else:
                    #print'Seam 0 is available'
                    seams.append(0)
                    #print'Seam 1 is available'
                    seams.append(1)
            if (qubitGroup+1)%cellLimit == 0:
                #print'Qubit '+str(qubit)+' is in the right column'
                if [qubit-qubitsPerCell,qubit] in outsideCouplers:
                    #print'Seam '+str(gridSize)+' is available'
                    seams.append(gridSize)
                else:
                    #print'Seam '+str(gridSize-1)+' is available'
                    seams.append(gridSize-1)
                    #print'Seam '+str(gridSize)+' is available'
                    seams.append(gridSize)
        qubitBase=(qubitGroup-1)*halfQPC
        for i in range(halfQPC):
            if takenQubits[qubitBase+i]==0:
                #print'Qubit '+str(qubit)+' can exit from the top or bottom'
                #print'Seam '+str(int(qubitCell/gridSize)+gridSize+1)+' is available'
                seams.append(int(qubitCell/gridSize)+gridSize+1)
                #print'Seam '+str(int(qubitCell/gridSize)+gridSize+2)+' is available'
                seams.append(int(qubitCell/gridSize)+gridSize+2)
                break
    return seams

#vacancy[North,West,East,South]
def perimeterVacancy(takenQubits,rows,columns):
    vacancy=[[] for i in range(4)]
    count=0
    done=False
    while not done:
        for i in range(len(rows)):
            for k in range(len(rows[i])):
                if takenQubits[rows[i][k]]==1:
                    done=True
                    break
            if done:
                break
            else:
                count+=1
    #print str(count)+' empty rows at the top'
    vacancy[0]=count

    count=0
    done=False
    while not done:
        for i in range(len(rows)):
            for k in range(len(rows[len(rows)-1-i])):
                if takenQubits[rows[len(rows)-1-i][k]]==1:
                    done=True
                    break
            if done:
                break
            else:
                count+=1
    #printstr(count)+' empty rows at the bottom'
    vacancy[3]=count

    count=0
    done=False
    while not done:
        for i in range(len(columns)):
            for k in range(len(columns[i])):
                if takenQubits[columns[i][k]]==1:
                    done=True
                    break
            if done:
                break
            else:
                count+=1
    #printstr(count)+' empty columns on the left'
    vacancy[1]=count

    count=0
    done=False
    while not done:
        for i in range(len(columns)):
            for k in range(len(columns[len(columns)-1-i])):
                if takenQubits[columns[len(columns)-1-i][k]]==1:
                    done=True
                    break
            if done:
                break
            else:
                count+=1
    #printstr(count)+' empty columns on the right'
    vacancy[2]=count

    return vacancy

#Returns the seam index that qubit1 and qubit2 can both find
def findSeam(qubit1,qubit2,takenQubits,outsideCouplers,vacancy):

    takenQubitsCopy=list(takenQubits)
    adjQub1=[[] for i in range(maxPaths)]
    adjQubOld1=orderQubits(innerAdjacentQubits(qubit1,takenQubitsCopy),takenQubitsCopy)+orderQubits(outerAdjacentQubits(qubit1,takenQubitsCopy),takenQubitsCopy)
    seam1=availableSeams(qubit1,takenQubitsCopy,outsideCouplers)
    adjQub2=[[] for i in range(maxPaths)]
    adjQubOld2=orderQubits(innerAdjacentQubits(qubit2,takenQubitsCopy),takenQubitsCopy)+orderQubits(outerAdjacentQubits(qubit2,takenQubitsCopy),takenQubitsCopy)
    seam2=availableSeams(qubit2,takenQubitsCopy,outsideCouplers)

    for i in range(len(adjQubOld1)):
        adjQub1[i]=list([adjQubOld1[i]])
    #for i in range(len(qubitList)):
    #    adjQub1[i]=list([qubitList[i]])
    for i in range(len(adjQubOld2)):
        adjQub2[i]=list([adjQubOld2[i]])

    noDeadEnd1=True
    noDeadEnd2=True
    tries=0
    while noDeadEnd1 or noDeadEnd2:
        tries+=1
        if tries > seamBreaker:
            #print 'qubit1=',qubit1
            #print '\nqubit2=',qubit2
            #print '\ntakenQubitsCopy=',takenQubitsCopy
            #print '\noutsideCouplers=',outsideCouplers
            #print '\nvacancy=',vacancy
            #while True:
            #    stop=True
            return -1
        if noDeadEnd1:
            #Qubit 1
            for i in range(len(adjQub1)):
                if len(adjQub1[i])==0:
                    break
                for j in range(len(adjQub1[i])-1):
                    takenQubitsCopy[adjQub1[i][j]]=1

                temp1=availableSeams(adjQub1[i][-1],takenQubitsCopy,outsideCouplers)
                done=False
                while True:
                    for j in range(len(temp1)):
                        if (temp1[j] <= gridSize and (vacancy[1]+vacancy[2]) == 0) or (temp1[j] > gridSize and (vacancy[0]+vacancy[3]) == 0):
                            temp1.pop(j)
                            break
                        if j == len(temp1)-1:
                            done=True
                    if done or len(temp1)==0:
                        break
                for j in range(len(temp1)):
                    seam1.append(temp1[j])
                for j in seam1:
                    if j in seam2:
                        if j > gridSize and vacancy[0]+vacancy[3] > 0:
                            return j
                        if j <= gridSize and vacancy[1]+vacancy[2] > 0:
                            return j

                for j in range(len(adjQub1[i])-1):
                    takenQubitsCopy[adjQub1[i][j]]=0
        if noDeadEnd2:
            #Qubit 2
            for i in range(len(adjQub2)):
                if len(adjQub2[i])==0:
                    break
                for j in range(len(adjQub2[i])-1):
                    takenQubitsCopy[adjQub2[i][j]]=1

                temp2=availableSeams(adjQub2[i][-1],takenQubitsCopy,outsideCouplers)
                done=False
                while True:
                    for j in range(len(temp2)):
                        if (temp2[j] <= gridSize and (vacancy[1]+vacancy[2]) == 0) or (temp2[j] > gridSize and (vacancy[0]+vacancy[3]) == 0):
                            temp2.pop(j)
                            break
                        if j == len(temp2)-1:
                            done=True
                    if done or len(temp2)==0:
                        break
                for j in range(len(temp2)):
                    seam2.append(temp2[j])
                for j in seam2:
                    if j in seam1:
                        if j > gridSize and vacancy[0]+vacancy[3] > 0:
                            return j
                        if j <= gridSize and vacancy[1]+vacancy[2] > 0:
                            return j

                for j in range(len(adjQub2[i])-1):
                    takenQubitsCopy[adjQub2[i][j]]=0
        if noDeadEnd1:
            #Copy adjQub1 to adjQubOld1
            adjQubOld1=[]
            for i in range(len(adjQub1)):
                if adjQub1[i] != []:
                    adjQubOld1.append(adjQub1[i])
                else:
                    break

            #Erase adjQub1
            for i in range(len(adjQub1)):
                if adjQub1[i]!=[]:
                    adjQub1[i]=[]
                else:
                    break

            #Track a dead end
            noDeadEnd1=False
            #For each route...
            for i in range(len(adjQubOld1)):
                #Mark the leading qubits as taken...
                for j in range(len(adjQubOld1[i])-1):
                    takenQubitsCopy[adjQubOld1[i][j]]=1

                #And find the qubits adjacent to the last on this list.
                #tempAdjQub=denseAdjacentQubits(adjQubOld1[i][-1],takenQubitsCopy)
                tempAdjQub=orderQubits(innerAdjacentQubits(adjQubOld1[i][-1],takenQubitsCopy),takenQubitsCopy)+orderQubits(outerAdjacentQubits(adjQubOld1[i][-1],takenQubitsCopy),takenQubitsCopy)

                #Unmark the leading qubits
                for j in range(len(adjQubOld1[i])-1):
                    takenQubitsCopy[adjQubOld1[i][j]]=0
                #Place this new list in the first empty spot in adjQub1
                for j in range(len(tempAdjQub)):
                    noDeadEnd1=True
                    for k in range(len(adjQub1)):
                        if adjQub1[k]==[]:
                            adjQub1[k]=list(adjQubOld1[i])
                            adjQub1[k].append(tempAdjQub[j])
                            break
        if noDeadEnd2:
            #Copy adjQub2 to adjQubOld2
            adjQubOld2=[]
            for i in range(len(adjQub2)):
                if adjQub2[i] != []:
                    adjQubOld2.append(adjQub2[i])
                else:
                    break

            #Erase adjQub2
            for i in range(len(adjQub2)):
                if adjQub2[i]!=[]:
                    adjQub2[i]=[]
                else:
                    break

            #Track a dead end
            noDeadEnd2=False
            #For each route...
            for i in range(len(adjQubOld2)):
                #Mark the leading qubits as taken...
                for j in range(len(adjQubOld2[i])-1):
                    takenQubitsCopy[adjQubOld2[i][j]]=1

                #And find the qubits adjacent to the last on this list.
                #tempAdjQub=denseAdjacentQubits(adjQubOld2[i][-1],takenQubitsCopy)
                tempAdjQub=orderQubits(innerAdjacentQubits(adjQubOld2[i][-1],takenQubitsCopy),takenQubitsCopy)+orderQubits(outerAdjacentQubits(adjQubOld2[i][-1],takenQubitsCopy),takenQubitsCopy)

                #Unmark the leading qubits
                for j in range(len(adjQubOld2[i])-1):
                    takenQubitsCopy[adjQubOld2[i][j]]=0
                #Place this new list in the first empty spot in adjQub2
                for j in range(len(tempAdjQub)):
                    noDeadEnd2=True
                    for k in range(len(adjQub2)):
                        if adjQub2[k]==[]:
                            adjQub2[k]=list(adjQubOld2[i])
                            adjQub2[k].append(tempAdjQub[j])
                            break
    #print '\n\n\n\n\nBoth Dead Ends!!!!\n\n\n\n\n'
    return -1

#Creates a new column or row in place of seam index, and updates all qubit references.
#Also returns a list of qubits to add to takenQubits and disable
def openSeam(seam,vacancy,rows,columns,qubitList,qubitPath,qubitIndex,qubitCouplers):
    toAppend=[]
    #print 'Opening Seam...'

    qubitListCopy=[[] for i in range(len(qubitList))]
    for i in range(len(qubitList)):
        if type(qubitList[i][0]) is list:
            tempList=[[] for j in range(len(qubitList[i]))]
            for j in range(len(qubitList[i])):
                tempList[j]=list(qubitList[i][j])
            qubitListCopy[i]=list(tempList)
        else:
            qubitListCopy[i]=list(qubitList[i])

    qubitPathCopy=[[] for i in range(len(qubitPath))]
    for i in range(len(qubitPath)):
        if type(qubitPath[i][0]) is list:
            tempPath=[[] for j in range(len(qubitPath[i]))]
            for j in range(len(qubitPath[i])):
                tempPath[j]=list(qubitPath[i][j])
            qubitPathCopy[i]=list(tempPath)
        else:
            qubitPathCopy[i]=list(qubitPath[i])

    qubitIndexCopy=[[] for i in range(len(qubitIndex))]
    for i in range(len(qubitIndex)):
        qubitIndexCopy[i]=list(qubitIndex[i])

    qubitCouplersCopy=[[] for i in range(len(qubitCouplers))]
    for i in range(len(qubitCouplers)):
        qubitCouplersCopy[i]=list(qubitCouplers[i])

    couplers=seamCouplers()
    #Horizontal Seam
    if seam > gridSize:
        if vacancy[0]==0 and vacancy[3]==0:
            #print 'NO SPACE'
            toAppend=-1
            return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
        #Shift down
        if vacancy[0] <= vacancy[3]:
            #seam < 2*gridSize
            if 2*gridSize-seam <= 0:
                #print 'ERROR'
                toAppend=-1
                return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
            #print 'Shift from row '+str(seam-gridSize-1)+' down'
            for i in range(2*gridSize-seam):

                #print 'Alter Row '+str(gridSize-2-i)
                #for j in range(len(rows[gridSize-2-i])):

                #Qubit Lists
                for j in range(len(qubitListCopy)):
                    for k in range(len(qubitListCopy[j])):
                        try:
                            for l in range(len(qubitListCopy[j][k])):
                                if qubitListCopy[j][k][l] in rows[gridSize-2-i]:
                                    qubitListCopy[j][k][l]+=rowDifference
                        except TypeError:
                            if qubitListCopy[j][k] in rows[gridSize-2-i]:
                                qubitListCopy[j][k]+=rowDifference
                    #print 'qubitList= '+str(qubitListCopy[j])

                #Qubit Paths
                for j in range(len(qubitPathCopy)):
                    for k in range(len(qubitPathCopy[j])):
                        try:
                            for l in range(len(qubitPathCopy[j][k])):
                                if qubitPathCopy[j][k][l] in rows[gridSize-2-i]:
                                    qubitPathCopy[j][k][l]+=rowDifference
                        except TypeError:
                            if qubitPathCopy[j][k] in rows[gridSize-2-i]:
                                    qubitPathCopy[j][k]+=rowDifference

                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*rowDifference or abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*qubitsPerCell:
                                            qubitPathCopy[j][k].insert(l+1,(qubitPathCopy[j][k][l]+qubitPathCopy[j][k][l+1])/2)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*rowDifference or abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*qubitsPerCell:
                                    qubitPathCopy[j].insert(k+1,(qubitPathCopy[j][k]+qubitPathCopy[j][k+1])/2)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if qubitPathCopy[j][k][l]==qubitPathCopy[j][k][l+1]:
                                            qubitPathCopy[j][k].pop(l)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if qubitPathCopy[j][k]==qubitPathCopy[j][k+1]:
                                    qubitPathCopy[j].pop(k)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    #print 'qubitPath= '+str(qubitPathCopy[j])

                #Qubit Indices
                for j in range(len(qubitIndexCopy)):
                    for k in rows[gridSize-2-i]:
                        if qubitIndexCopy[j][k]==1:
                            qubitIndexCopy[j][k+rowDifference]=1
                            qubitIndexCopy[j][k]=0
                    #print 'qubitIndex= '+str(t)

                #Qubit Couplers
                for j in range(len(qubitCouplersCopy)):
                    for k in range(len(qubitCouplersCopy[j])):
                        if qubitCouplersCopy[j][k] in rows[gridSize-2-i]:
                            qubitCouplersCopy[j][k]+=rowDifference


        #Shift up
        else:
            #seam > gridSize+2 seam-grid-2
            if seam-gridSize-2 <= 0:
                #print 'ERROR'
                toAppend=-1
                return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
            #print 'Shift up to row '+str(seam-gridSize-2)+' up'
            for i in range(seam-gridSize-2):
            #    for k in range(len(rows[i])):
                #print 'Alter Row '+str(i+1)

                #Qubit Lists
                for j in range(len(qubitListCopy)):
                    for k in range(len(qubitListCopy[j])):
                        try:
                            for l in range(len(qubitListCopy[j][k])):
                                if qubitListCopy[j][k][l] in rows[i+1]:
                                    qubitListCopy[j][k][l]-=rowDifference
                        except TypeError:
                            if qubitListCopy[j][k] in rows[i+1]:
                                qubitListCopy[j][k]-=rowDifference
                    #print 'qubitList= '+str(qubitListCopy[j])

                #Qubit Paths
                for j in range(len(qubitPathCopy)):
                    for k in range(len(qubitPathCopy[j])):
                        try:
                            for l in range(len(qubitPathCopy[j][k])):
                                if qubitPathCopy[j][k][l] in rows[i+1]:
                                    qubitPathCopy[j][k][l]-=rowDifference
                        except TypeError:
                            if qubitPathCopy[j][k] in rows[i+1]:
                                    qubitPathCopy[j][k]-=rowDifference

                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*rowDifference or abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*qubitsPerCell:
                                            qubitPathCopy[j][k].insert(l+1,(qubitPathCopy[j][k][l]+qubitPathCopy[j][k][l+1])/2)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*rowDifference or abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*qubitsPerCell:
                                    qubitPathCopy[j].insert(k+1,(qubitPathCopy[j][k]+qubitPathCopy[j][k+1])/2)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if qubitPathCopy[j][k][l]==qubitPathCopy[j][k][l+1]:
                                            qubitPathCopy[j][k].pop(l)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if qubitPathCopy[j][k]==qubitPathCopy[j][k+1]:
                                    qubitPathCopy[j].pop(k)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    #print 'qubitPath= '+str(qubitPathCopy[j])

                #Qubit Indices
                for j in range(len(qubitIndexCopy)):
                    for k in rows[i+1]:
                        if qubitIndexCopy[j][k]==1:
                            qubitIndexCopy[j][k-rowDifference]=1
                            qubitIndexCopy[j][k]=0
                    #print 'qubitIndex= '+str(t)

                #Qubit Couplers
                for j in range(len(qubitCouplersCopy)):
                    for k in range(len(qubitCouplersCopy[j])):
                        if qubitCouplersCopy[j][k] in rows[i+1]:
                            qubitCouplersCopy[j][k]-=rowDifference

    #Vertical Seam
    else:
        if vacancy[1]==0 and vacancy[2]==0:
            #print 'NO SPACE'
            toAppend=-1
            return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
        #Shift right
        if vacancy[1] <= vacancy[2]:
            #seam < gridSize-1
            if gridSize-seam-1 <= 0:
                #print 'ERROR'
                toAppend=-1
                return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
            #print 'Shift from column '+str(seam)+' right'
            for i in range(gridSize-seam-1):
            #    for k in range(len(columns[i])):
                #print 'Alter Column '+str(gridSize-2-i)

                #Qubit Lists
                for j in range(len(qubitListCopy)):
                    for k in range(len(qubitListCopy[j])):
                        try:
                            for l in range(len(qubitListCopy[j][k])):
                                if qubitListCopy[j][k][l] in rows[gridSize-2-i]:
                                    qubitListCopy[j][k][l]+=qubitsPerCell
                        except TypeError:
                            if qubitListCopy[j][k] in columns[gridSize-2-i]:
                                qubitListCopy[j][k]+=qubitsPerCell
                    #print 'qubitList= '+str(qubitListCopy[j])

                #Qubit Paths
                for j in range(len(qubitPathCopy)):
                    for k in range(len(qubitPathCopy[j])):
                        try:
                            for l in range(len(qubitPathCopy[j][k])):
                                if qubitPathCopy[j][k][l] in columns[gridSize-2-i]:
                                    qubitPathCopy[j][k][l]+=qubitsPerCell
                        except TypeError:
                            if qubitPathCopy[j][k] in columns[gridSize-2-i]:
                                    qubitPathCopy[j][k]+=qubitsPerCell

                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*rowDifference or abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*qubitsPerCell:
                                            qubitPathCopy[j][k].insert(l+1,(qubitPathCopy[j][k][l]+qubitPathCopy[j][k][l+1])/2)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*rowDifference or abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*qubitsPerCell:
                                    qubitPathCopy[j].insert(k+1,(qubitPathCopy[j][k]+qubitPathCopy[j][k+1])/2)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if qubitPathCopy[j][k][l]==qubitPathCopy[j][k][l+1]:
                                            qubitPathCopy[j][k].pop(l)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if qubitPathCopy[j][k]==qubitPathCopy[j][k+1]:
                                    qubitPathCopy[j].pop(k)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    #print 'qubitPath= '+str(qubitPathCopy[j])

                #Qubit Indices
                for j in range(len(qubitIndexCopy)):
                    for k in columns[gridSize-2-i]:
                        if qubitIndexCopy[j][k]==1:
                            qubitIndexCopy[j][k+qubitsPerCell]=1
                            qubitIndexCopy[j][k]=0
                    #print 'qubitIndex= '+str(t)

                #Qubit Couplers
                for j in range(len(qubitCouplersCopy)):
                    for k in range(len(qubitCouplersCopy[j])):
                        if qubitCouplersCopy[j][k] in columns[gridSize-2-i]:
                            qubitCouplersCopy[j][k]+=qubitsPerCell

        #Shift left
        else:
            #seam > 1
            if seam-1 <= 0:
                #print 'ERROR'
                toAppend=-1
                return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]
            #print 'Shift up to column '+str(seam-1)+' left'
            for i in range(seam-1):
            #    for k in range(len(columns[i])):
                #print 'Alter Column '+str(i+1)

                #Qubit Lists
                for j in range(len(qubitListCopy)):
                    for k in range(len(qubitListCopy[j])):
                        try:
                            for l in range(len(qubitListCopy[j][k])):
                                if qubitListCopy[j][k][l] in rows[i+1]:
                                    qubitListCopy[j][k][l]-=qubitsPerCell
                        except TypeError:
                            if qubitListCopy[j][k] in columns[i+1]:
                                qubitListCopy[j][k]-=qubitsPerCell
                    #print 'qubitList= '+str(qubitListCopy[j])

                #Qubit Paths
                for j in range(len(qubitPathCopy)):
                    for k in range(len(qubitPathCopy[j])):
                        try:
                            for l in range(len(qubitPathCopy[j][k])):
                                if qubitPathCopy[j][k][l] in columns[i+1]:
                                    qubitPathCopy[j][k][l]-=qubitsPerCell
                        except TypeError:
                            if qubitPathCopy[j][k] in columns[i+1]:
                                    qubitPathCopy[j][k]-=qubitsPerCell

                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*rowDifference or abs(qubitPathCopy[j][k][l]-qubitPathCopy[j][k][l+1])==2*qubitsPerCell:
                                            qubitPathCopy[j][k].insert(l+1,(qubitPathCopy[j][k][l]+qubitPathCopy[j][k][l+1])/2)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*rowDifference or abs(qubitPathCopy[j][k]-qubitPathCopy[j][k+1])==2*qubitsPerCell:
                                    qubitPathCopy[j].insert(k+1,(qubitPathCopy[j][k]+qubitPathCopy[j][k+1])/2)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    done=False
                    while not done:
                        for k in range(len(qubitPathCopy[j])):
                            try:
                                done=False
                                while not done:
                                    for l in range(len(qubitPathCopy[j][k])):
                                        if qubitPathCopy[j][k][l]==qubitPathCopy[j][k][l+1]:
                                            qubitPathCopy[j][k].pop(l)
                                            break
                                        if l==len(qubitPathCopy[j][k])-2:
                                            done=True
                                            break
                            except TypeError:
                                if qubitPathCopy[j][k]==qubitPathCopy[j][k+1]:
                                    qubitPathCopy[j].pop(k)
                                    break
                                if k==len(qubitPathCopy[j])-2:
                                    done=True
                                    break
                    #print 'qubitPath= '+str(qubitPathCopy[j])

                #Qubit Indices
                for j in range(len(qubitIndexCopy)):
                    for k in columns[i+1]:
                        if qubitIndexCopy[j][k]==1:
                            qubitIndexCopy[j][k-qubitsPerCell]=1
                            qubitIndexCopy[j][k]=0
                    #print 'qubitIndex= '+str(t)

                #Qubit Couplers
                for j in range(len(qubitCouplersCopy)):
                    for k in range(len(qubitCouplersCopy[j])):
                        if qubitCouplersCopy[j][k] in columns[i+1]:
                            qubitCouplersCopy[j][k]-=qubitsPerCell
    done=False
    while not done:
        for i in range(len(qubitCouplersCopy)):
            if abs(qubitCouplersCopy[i][0]-qubitCouplersCopy[i][1])==2*rowDifference or abs(qubitCouplersCopy[i][0]-qubitCouplersCopy[i][1])==2*qubitsPerCell:
                high=max(qubitCouplersCopy[i][0],qubitCouplersCopy[i][1])
                low=min(qubitCouplersCopy[i][0],qubitCouplersCopy[i][1])
                mid=(qubitCouplersCopy[i][0]+qubitCouplersCopy[i][1])/2
                toAppend.append(mid)
                qubitCouplersCopy.insert(i+1,[mid,high])
                qubitCouplersCopy.insert(i+1,[low,mid])
                qubitCouplersCopy.pop(i)
                break
            if i == len(qubitCouplersCopy)-1:
                done=True
                break

    #print 'qubitCouplers= '+str(qubitCouplersCopy)
    #print 'toAppend= '+str(toAppend)
    #print '\n'

    return [qubitListCopy,qubitPathCopy,qubitIndexCopy,qubitCouplersCopy,toAppend]

#Defines the couplers in each seam
def seamCouplers():
    couplers=[[] for i in range(2*(gridSize+1))]
    for i in range(1,gridSize):
        base=halfQPC+(i-1)*qubitsPerCell
        for k in range(gridSize):
            for j in range(halfQPC):
                couplers[i].append([base+j,base+j+qubitsPerCell])
            base+=gridSize*qubitsPerCell

    for i in range(1,gridSize):
        base=gridSize*qubitsPerCell*(i-1)
        for k in range(gridSize):
            for j in range(halfQPC):
                couplers[i+gridSize+1].append([base+j,base+j+gridSize*qubitsPerCell])
            base+=qubitsPerCell

    return couplers
