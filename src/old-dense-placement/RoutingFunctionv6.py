#Miguel Aroca-Ouellette
#Last modified: 10/07/2013
#note: routes are no longer loaded from an array, they are passed from the calling function
#routes can start and end on same qubit

import sys, math, os

 ### change path to current directory
os.chdir(os.path.dirname(os.path.realpath(__file__)))

####MORE GLOBAL
global numberOfQubits

numberOfQubits=512
####MORE GLOBAL

def initVar():

    global paths,is_used,is_shared,baseCOST,hist_cost,sharingCOST,allPaths,breaker,sharingINC
    global historyCOST,forgetRATE,J,h,numberOfQubits,routes,curr_used,numQu,breakCOST

    routes= []
    J=dict()
    h=[0]*numberOfQubits

    baseCOST=1;
    is_used=[0]*numberOfQubits      #value at index represents number of paths using the node
    is_shared=[0]*numberOfQubits    #0 if not shared, 1 if shared
    curr_used=[0]*numberOfQubits    #prevents overlapping paths being created
    sharingCOST=1.0;
    sharingINC=1.0;
    historyCOST=1.0;
    forgetRATE=0.001                 #smaller number=forgets slower
    breakCOST=9999                  #value returned if cannot be routed (maybe make it so its a sum?)
    paths=[]
    allPaths=[] #has no cost attached to it
    numQu=[]

    #check if breaker has been defined
    try:
        temp=breaker
    except NameError:
        breaker=150     #program will return error after this many itirations

    try:
        temp=hist_cost
    except NameError:
        hist_cost=[0]*numberOfQubits

######################## MAIN  ########################

def Routing(inputRoutes,coupler_file,filename, tofile=True):
    global sharingCOST,allPaths,is_used,is_shared,routes,numQu

    ### load data ###
    initVar()

    try:
        pass1=couplers[1]
        #else file has already been initialized
    except NameError:
        loadArrays(coupler_file,[None])      #input here coupler file you want to load from

    routes=[]   #reset routes
    for i in xrange(len(inputRoutes)):
        routes.append(list(inputRoutes[i]))

    ### Checks for Errors ###
    for i in xrange(len(routes)):
        temp_routes=routes.pop(i)
        if temp_routes in routes:
            print "Routing ERROR: DUPLICATE ROUTE(S): ",temp_routes
            sys.exit()
        if len(temp_routes)!=2:
            print "Routing ERROR: ROUTE ",temp_routes," IS NOT OF LENGTH 2"
            sys.exit()
        if temp_routes[0]>=numberOfQubits or temp_routes[1]>numberOfQubits or temp_routes[0]<0 or temp_routes[1]<0:
            print "Routing ERROR: INVALID ROUTE ",temp_routes
            sys.exit()
        routes.insert(i,temp_routes)

    check=0
    for i in xrange(len(inputRoutes)):
        for j in xrange(2):
            for k in xrange(len(couplers_group[inputRoutes[i][j]])):
                if hist_cost[couplers_group[inputRoutes[i][j]][k]]!=10000000:
                    check=1
            if check==0 and not(inputRoutes[i][1-j] in couplers_group[inputRoutes[i][j]]):
                #print "disableQubits ERROR: All the qubits around qubit ",inputRoutes[i][j],"have been disabled, no route could be found. BreakCOST returned."
                return breakCOST
            check=0

    ### Checks for Errors ###

    ### order routes from shortest to longest###
    #for i in xrange(len(routes)):
    #    temp_rou=routes.pop(i)
    #    temp_rou.insert(0,estCoup(temp_rou[0],temp_rou[1]))
    #    routes.insert(i,temp_rou)
    #
    #routes.sort()
    #for i in xrange(len(routes)):
    #    temp_rou=routes.pop(i)
    #    temp_rou=temp_rou[1:]
    #    routes.insert(i,temp_rou)
    #print "Initial ordered routes:",routes
    ### order routes ###

    ### Negotiated Congestion and Routing ###
    is_shared[0]=1
    while sum(is_shared)!=0:
        is_shared=[0]*numberOfQubits
        allPaths=[]     #RIP UP ROUTES
        is_used=[0]*numberOfQubits  #since all routes are ripped up, nothing is used
        for i in xrange(len(routes)):
            allPaths.append(bestPath(routes[i]))
            genUsedShared(allPaths) #updates shared

        #print "SHARED:",is_shared
        #print "USED:",is_used
        #print "hist_cost:",hist_cost
        #print"ALLPATHS: ", allPaths, "\n"
        #for i in xrange (numberOfQubits):
        #    print nodeCost(i),
        #print "\n"
        genHist()
        #reorderRoutes() #reorders routes to favor shared and shorter routes
        """MIGHT WANT TO REMOVE LATER"""
        #if sharingCOST!=1:
        ###print "On iteration #",sharingCOST," we routed",round(perRoute(allPaths)*100.00,2),"% of all routes"
       # numQu=0
        #numCheck=[]
        #append2=numCheck.append
        #for k in xrange(len(allPaths)):
        #    for m in xrange(len(allPaths[k])):
        #        if allPaths[k][m] not in numCheck:
        #            append2(allPaths[k][m])
        #            numQu+=1
        #if sharingCOST%5==0:
        #     print "On iteration #",sharingCOST," we routed",round(perRoute(allPaths)*100.00,2),"% of all routes. Used ",numQu," qubits."
        if tofile:
			writeToFile(filename)
        sharingCOST+=sharingINC   #this has to be increased slowly
        if sharingCOST>breaker:
            break

    ### Negotiated Congestion and Routing ###

    ### Check to make sure that paths don't use a disabled qubit ###
    try:
        for i in xrange(len(allPaths)):
            for j in xrange(len(allPaths[i])-2):
                if allPaths[i][j+1] in qubits_off:
                    allPaths=[]
                    return breakCOST
    except NameError:
        pass
    ### Check to make sure that paths don't use a disabled qubit ###

    ### Return total cost of all routes ###
    totCost=0.0
    for i in xrange(len(allPaths)):
        for j in xrange(1,len(allPaths[i])):
            totCost+=nodeCost(allPaths[i][j])
    ### Return total cost of all routes ###

    if sharingCOST>breaker:
        print "Routing BREAK ERROR: the routing timed out and no solution was found"
        #return [0,sum(numQu)/len(numQu),breakCOST]
        return breakCOST


    """MIGHT WANT TO REMOVE LATER"""
    #if sharingCOST!=2:
    ###print "Routing SUCCESS: Took ",sharingCOST-1,"re-routes for all paths to be routed."
    #return [1,sum(numQu)/len(numQu),totCost]
    return totCost

######################## MAIN  ########################

######################## LOAD ARRAYS  ########################

def loadArrays(input,dis_list):
    #usually 'couplers.txt'
    #dist_list is a list of lists
    global couplers,couplers_group

    couplers=[]
    couplers_group=[None]*numberOfQubits

    ###Check for errors
    try:
        if dis_list!=[None]:
            check=dis_list[0][0]
    except TypeError:
        print "loadArrays DISABLE COUPLERS ERROR: Invalid input (",dis_list,")"
        return

    try:
        file_val=int(input[-6])
        file_val=file_val=int(input[-6:-4])
    except ValueError:
        try:
            file_val=int(input[-5])
        except ValueError:
            if input=='couplers.txt':
                file_val=8
            else:
                print "loadArrays ERROR: Coupler file (",input,") does not correspond with number of qubits (",numberOfQubits,")"
                sys.exit()

    if pow(numberOfQubits/8,0.5)!=file_val:
        print "loadArrays ERROR: Coupler file (",input,") does not correspond with number of qubits (",numberOfQubits,")"
        sys.exit()

    ###Check for errors

    ### Load Files ###
    filename = input
    file = open(filename,'r')
    data = file.readlines()

    ### Create Coupler Array ###
    check=0
    append1=couplers.append
    for i in xrange(len(data)):
        j=data[i].index(' ')
        temp_check=[int(data[i][0:j])]+[int(data[i][j+1:])]
        if temp_check in dis_list or temp_check[::-1] in dis_list:
            check+=1
            continue
        else:
            append1(temp_check)

    if check!=len(dis_list) and dis_list!=[None]:
        print "loadArrays DISABLE COUPLERS ERROR: Input coupler does not exit (",dis_list,")"
        couplers=[]
        couplers_group=[None]*numberOfQubits
        return

    file.close()

    for i in xrange(numberOfQubits):
        temp_couplers2 = []
        temp_append=temp_couplers2.append
        for j in xrange(len(couplers)):
            if couplers[j][0] == i:
                temp_append(couplers[j][1])
            if couplers[j][1] == i:
                temp_append(couplers[j][0])
        couplers_group[i]=temp_couplers2


def disableQubits(input):
    #input is list
    global hist_cost,qubits_off

    qubits_off=list(input)
    hist_cost=[0]*numberOfQubits
    for i in xrange(len(input)):
        hist_cost[input[i]]=10000000 #some arbitrarily high number

def resetQubits():
    global hist_cost

    hist_cost=[0]*numberOfQubits

def disableCouplers(input,coupler_file):
    global numberOfQubits

    loadArrays(coupler_file,input)

def resetCouplers(coupler_file):
    global numberOfQubits

    loadArrays(coupler_file,[None])

def setNumberQubits(num):
    """
    Number of Qubits is assumed to be 512 unless otherwise stated.
    WARNING: Make sure you load the appropriate coupler file.
    """
    global numberOfQubits,couplers,couplers_group

    try:
        del couplers,couplers_group #so that you don't accidentally use an old couplers file
    except NameError:
        pass

    if math.pow(num/8,0.5)%1!=0:
        print "setNumberQubits ERROR: Invalid input (",num,") - does not correspond to a square grid."
    numberOfQubits=num

######################## LOAD ARRAYS  ########################

######################## BEST PATH FUNCTIONS ########################

### cost of node ###
def nodeCost(index):
    global is_used, is_shared,hist_cost,sharingCOST

    if is_shared[index]>0:
        output=is_used[index]*(baseCOST+hist_cost[index])*(sharingCOST)
    else:
        output=is_used[index]*(baseCOST+hist_cost[index])
    return output

### cost of node ###

### function which returns index of cheapest path ####
def indexOfPath():
    global paths

    index=None
    min_val=0
    for i in xrange(len(paths)):
        if (paths[i][0]<min_val) or i==0:
            min_val=paths[i][0]
            index=i
    if index==None:
        print "InternalRoutingError: ERROR INDEX RETURNS NONE. PATHS ARE: ",paths
    return index

### function which returns index of cheapest path ####

### function expands specified path and increases cost ###
def expandPath(index):
    #NOTE: first value of paths is COST
    global paths, baseCOST,is_used,curr_used

    hold_paths=list(paths)  #locks path for iteration
    temp_extend=couplers_group[hold_paths[index][-1]]     #find qubits adjacent to last qubit in path
    insert1=paths.insert
    for j in xrange(len(temp_extend)):
        if temp_extend[j] in hold_paths[index][1:-1]:     #prevents node loops from forming, each node can only be visited once per path
            continue

        temp_new=list(hold_paths[index])
        temp_new.append(temp_extend[j])
        if curr_used[temp_extend[j]]>0: #if a path is already passing through a node then don't expand as other path will be shorter
            continue

        curr_used[temp_extend[j]]=1
        is_used[temp_extend[j]]+=1        #marks the node as used for the purpose of finding the cost
        temp_new[0]+=nodeCost(temp_extend[j])
        is_used[temp_extend[j]]-=1       #we do not know if the node is used until we reach the goal

        insert1(j,temp_new)     #doesn't really matter where you insert it
    paths.remove(hold_paths[index])     #removes old path

### function expands specified path and increases cost###

##### MAIN Best Path #####
def bestPath(input):
    #first value of input is start, 2nd value is goal
    global paths,curr_used

    if input[0]==input[1]:
        print "bestPath ERROR: start is same as goal!"
        return 0

    new_path=[0,input[0]]       #NOTE: first value of paths is COST
    paths.append(new_path)
    check=0      #input[1]+1   #initalize to value not equal to goal
    goal=input[1]

    curr_used[input[0]]=1
    while check!=1:
        #print paths
        expandPath(indexOfPath())
        for i in xrange(len(paths)):
            if paths[i][-1]==goal:
                check=1
                output=paths[i][1:]
                break           #stops as soon as answer is found

    #resets values
    paths=[]
    curr_used=[0]*numberOfQubits

    return output

##### MAIN Best Path ######

######################## BEST PATH FUNCTIONS ########################

######################## NEGOTIATED CONGESTION FUNCTIONS ########################

### update values of is_shared and is_used ###
def genUsedShared(allpaths):
    #handles having shared starting and end qubits
    global is_used,is_shared,hist_cost

    is_used=[0]*numberOfQubits #resets is_used
    for i in xrange (numberOfQubits):
        for j in xrange(len(allpaths)):
            if i in allpaths[j]:
                is_used[i]+=1

    end_list=[0]*numberOfQubits
    for i in xrange(len(allpaths)):
        end_list[allpaths[i][0]]+=1
        end_list[allpaths[i][-1]]+=1

    for i in xrange (numberOfQubits):
        if end_list[i]>1:
            is_used[i]-=(end_list[i]-1)

    for i in xrange (numberOfQubits):
        if is_used[i]>1:
            is_shared[i]=1
        else:
            is_shared[i]=0

### update values of is_shared and is_used ###

### update values of hist_cost ###
def genHist():
    global hist_cost,histCOST,forgetRATE

    for i in xrange (numberOfQubits):
        if hist_cost[i]>0:
            hist_cost[i]=hist_cost[i]*math.exp(-forgetRATE)    #forgets history cost over time
        if is_shared[i]>0:
            hist_cost[i]+=historyCOST      #increases hist_cost so we don't get stuck in loops

### update values of hist_cost ###

######################## NEGOTIATED CONGESTION FUNCTIONS ########################

######################## OTHER FUNCTIONS ########################

### initializes couplers ###

def getCouplers():
    J=dict()
    for i in xrange(len(allPaths)):
        for j in xrange(len(allPaths[i])-1):
            J[(allPaths[i][j],allPaths[i][j+1])]=1

    return J

### intializes couplers ###

### finds % of paths routed ###

def perRoute(input):
    #returns decimal

    #problem 1:makes sure that first item in list is longest
    #CURRENT PROBLEM: SWITCHING ORDER OF ROUTES RESULTS IN DIFFERENT ANSWERS:
            # we want to make sure that we check longest route first, becasue this will be the one which eliminates thte most nodes longest->shortest

    per_paths=list(input)   #prevents overwrite
    per_paths=orderPaths(per_paths,1)
    total=len(per_paths)
    cross=0
    temp_shared=list(is_shared)
    for i in xrange(total):
        temp_cross=0
        for j in xrange(len(per_paths[i])):
            if temp_shared[per_paths[i][j]]==1:
                temp_shared[per_paths[i][j]]=0
                temp_cross=1
        if temp_cross==1:
            cross+=1

    output=(total-cross)/(total*1.0)  #cast to float
    return output

### finds % of paths routed ###

### order the inputted paths by number of couplers used###

def orderPaths(input,order):
    #based on number of couplers used
    #input=0 for ascending order,input=1 for descending order
    #much faster than using: sorted([[1,1],[2]],key=len)

    output=[]
    output=list(input)
    length=len(output)
    for i in xrange(length): #find all lengths
        temp_output=list(output.pop(i))     #because every list index poitns to a new list, don't want referencing
        temp_output.insert(0,len(temp_output)-1)
        output.insert(i,temp_output)
    #print allPaths

    output.sort()
    if order==1:
        output.reverse()

    for i in xrange (length):
        temp_output=output.pop(i)
        temp_output=temp_output[1:]   #removes length value
        output.insert(i,temp_output)

    return output

### order the inputted paths by number of couplers used###

### estimate num of couplers used ###

def estCoup(num1,num2):
    output=abs(num1/64-num2/64) #difference in rows
    output+=abs((num1%64)/8-(num2%64)/8)  #difference in columns

    #bring everything to cell 0
    num1=(num1%64)%8
    num2=(num2%64)%8

    #difference in qubits, only two options
    if abs(num1/4-num2/4)==0:
        output+=2
    else:
        output+=1

    return output

### estimate num of couplers used ###

### reorder routes ###

def reorderRoutes():
    global routes,allPaths

    shared=[]
    lone=[]
    for i in xrange(len(routes)):
        temp_all=allPaths[i]
        check=0
        for j in xrange(len(temp_all)):
            if is_shared[temp_all[j]]==1:
                shared.append(temp_all)
                check=1
                break
        if check==0:
            lone.append(temp_all)

    shared=orderPaths(shared,0)
    lone=orderPaths(lone,0)

    #convert from paths to routes
    routes=[]
    for i in xrange(len(shared)):
        temp_rou=shared[i]
        new_rou=[]
        new_rou.append(temp_rou[0])
        new_rou.append(temp_rou[-1])
        routes.append(new_rou)
    for i in xrange(len(lone)):
        temp_rou=lone[i]
        new_rou=[]
        new_rou.append(temp_rou[0])
        new_rou.append(temp_rou[-1])
        routes.append(new_rou)

### modifyRoute ###
def modifyRoute(old_qubit,new_qubit,replace):
    #now moves all routes connected to old qubit to a new qubit
    #need to have run routing beforehand
    global allPaths,sharingCOST,is_shared,is_used,hist_cost,breakerCOST


    ###store old variables
    origPaths=[]
    missPaths=[]    #stores value but with missing routes

    map(origPaths.append,list(allPaths))
    map(missPaths.append,list(allPaths))

    origSharing=sharingCOST
    orig_is_shared=list(is_shared)
    orig_is_used=list(is_used)
    orig_hist_cost=list(hist_cost)

    if new_qubit>=numberOfQubits or new_qubit<0:
        print "modifyRoute ERROR: NEW QUBIT (",new_qubit,") IS NOT WITHIN RANGE of [0,511]"
        return 0
    if new_qubit==old_qubit:
        print "modifyRoute ERROR: NEW QUBIT (",new_qubit,") IS EQUAL TO OLD QUBIT, NO CHANGE IN ROUTING WAS PERFORMED"
        return 0

    check=0
    new_routes=[]
    append1=new_routes.append
    for i in xrange(len(allPaths)):  #checks that old qubit is actually within allPaths as well as stores them
        if (allPaths[i][0]==old_qubit):
            if allPaths[i][-1]==new_qubit:
                print "modifyRoute ERROR: NEW QUBIT CAUSES ROUTE TO START AND END ON SAME QUBIT"
                return 0
            check=1
            missPaths.remove(allPaths[i])
            allPaths[i]=[allPaths[i][-1]]
            allPaths[i].insert(0,new_qubit)
            append1(list(allPaths[i]))
        if (allPaths[i][-1]==old_qubit):
            if allPaths[i][0]==new_qubit:
                print "modifyRoute ERROR: NEW QUBIT CAUSES ROUTE TO START AND END ON SAME QUBIT"
                return 0
            check=1
            missPaths.remove(allPaths[i])
            allPaths[i]=[allPaths[i][0]]
            allPaths[i].insert(len(allPaths[i]),new_qubit)
            append1(list(allPaths[i]))

    if check==0:
        print "modifyRoute ERROR: OLD QUBIT(",old_qubit,") DOES NOT EXIST WITHIN ROUTES!"
        return 0

    is_shared[0]+=1
    while sum(is_shared)!=0:
        is_shared[0]-=1
        allPaths=list(missPaths)    #restore values
        is_used=list(orig_is_used)
        for i in xrange(len(new_routes)):
            allPaths.append(bestPath(new_routes[i]))
            genUsedShared(allPaths) #updates shared
        genHist()
        sharingCOST+=sharingINC   #this has to be increased slowly
        if sharingCOST>breaker+origSharing:
            break

    ###newCost
    newCost=0.0
    for i in xrange(len(allPaths)-len(new_routes),len(allPaths)):
        for j in xrange(1,len(allPaths[i])):
            newCost+=nodeCost(allPaths[i][j])

    if sharingCOST>breaker+origSharing:
        print "modifyRoute BREAK ERROR: the routing timed out and no solution was found"
        newCost=breakCOST

    ###restore old values
    if replace==0 or sharingCOST>breaker+origSharing:
        sharingCOST=origSharing
        allPaths=list(origPaths)
        is_shared=list(orig_is_shared)
        is_used=list(orig_is_used)
        hist_cost=list(orig_hist_cost)

    #print "SHARED:",is_shared
    #print "USED:",is_used
    #print "hist_cost:",hist_cost
    #print "after replace or restored ",allPaths

    return newCost

### modifyRoute ###


### get route ### for testing purposes only
def getRoute():
    global allPaths

    if allPaths==[] or sharingCOST>breaker:
        return breakCOST
    return allPaths

### set breaker value ###
def setBreaker(input):
    global breaker
    breaker=input

### GET DENSITY ###
def getDensity():
    #returns list of length 64 with numbering from the least dense to the most dense tiles
    #if tiles are equally dense then will return the same number in the tiles

    output=[None]*(numberOfQubits/8)    #64
    sum_tile=[0]*(numberOfQubits/8)

    for i in xrange(len(allPaths)):
        for j in xrange(len(allPaths[i])):
            sum_tile[allPaths[i][j]/8]+=1

    min_val=min(sum_tile)
    count=0
    while min_val!=numberOfQubits+1:
        for j in range(len(sum_tile)):
          if sum_tile[j] == min_val:
             sum_tile[j]=numberOfQubits+1  #no sum could ever be greater than 8
             output[j]=count
        min_val=min(sum_tile)
        count+=1

    return output

def getPopularity():
    #returns list of length (numQubits) with start and end qubits with number of outgoing connections

    output=[0]*(numberOfQubits)

    for i in xrange(len(allPaths)):
        output[allPaths[i][0]]+=1
        output[allPaths[i][-1]]+=1

    return output

def writeToFile(filename):
    #filename.write('new:\n')
    filename.write('new:\r\n')
    for u in xrange(len(allPaths)):
        for v in xrange(len(allPaths[u])):
            if v==0:
                filename.write(str(allPaths[u][v]))
            else:
                filename.write(' '+str(allPaths[u][v]))
        #filename.write('\n')
        filename.write('\r\n')
    filename.write('shared:')
    for i in xrange(len(is_shared)):
        if is_shared[i]==1:
            filename.write(' '+str(i))
    #filename.write('\n')
    filename.write('\r\n')

def getNotRouted():
    #returns list of (numQubits) which shows number of unrouted paths from each start and qubit

    output=[0]*numberOfQubits

    for i in range(len(allPaths)):
        for j in range(len(allPaths[i])):
            if is_shared[allPaths[i][j]]==1:
                output[allPaths[i][0]]+=1
                output[allPaths[i][-1]]+=1
    return output

def getInTheWay():
    #returns list of (numQubits) which shows number of qubits going through each start/end qubit

    obstruct=[0]*numberOfQubits
    for i in xrange(len(allPaths)):
        obstruct[allPaths[i][0]]=is_used[allPaths[i][0]]
        obstruct[allPaths[i][-1]]=is_used[allPaths[i][-1]]

    return obstruct

def setRoute(input_routes,coupler_file):
    """
    Input is list of lists of completed routes and couplers file.
    This function WILL NOT CHECK IF YOUR ROUTES ARE CORRECT.
    Will reset all variables if Routing has been previously called.
    """
    global allPaths
        ### load data ###

    #check for errors:
    try:
        temp=input_routes[0][0]
    except TypeError:
        print "setRoutes ERROR: INVALID INPUT (",input_routes,") - Input must be list of lists."
        return 0

    initVar()

    try:
        pass1=couplers[1]
        #else file has already been initialized
    except NameError:
        loadArrays(coupler_file,[None])      #input here coupler file you want to load from

    allPaths=[]
    for i in range(len(input_routes)):
        allPaths.append(list(input_routes[i]))

    genUsedShared(allPaths)

######################## OTHER FUNCTIONS ########################
#disableQubits([0,64])
#file3=open('routes.txt','w')
##setRoute([[4,12,20]],'couplers.txt')
##print getRoute()
##setBreaker(5)
#setNumberQubits(392)
#print Routing([[0,64]],'couplers7.txt',file3)
#print getRoute()
##modifyRoute(20,5,1)
##print getRoute()

##file3=open('routes.txt','w')
##setNumberQubits(32)
##Routing([[7,15]],'couplers1.txt',file3)
##print getRoute()
