import math

def convert(index,old,new):
    """
    Converts the index in a square qubit array containing 'old' number of qubits to
    the corresponding index in a square qubit array containing 'new' number of qubits.
    """
    if index>=old:
        print "convertERROR: Invalid input - index(",index,") is greater or equal to old # of qubits (",old,")"
        return None
    sizeOld=math.sqrt(old/8)
    sizeNew=math.sqrt(new/8)

    if sizeOld%1!=0:
        print "convertERROR: Invalid input - old number of qubits (",old,") does not correspond to a square grid"
        return None
    if sizeNew%1!=0:
        print "convertERROR: Invalid input - new number of qubits (",new,") does not correspond a square grid"
        return None

    row=int((index/8)/sizeOld)
    col=int((index/8)%sizeOld)
    wrtOld=index%8

    if row>sizeNew-1 or col>sizeNew-1:
        print "convertERROR: Invalid input - index (",index,") is outside of new grid."
        return None

    tileNew=(sizeNew*row+col)*8
    return int(tileNew+wrtOld)

def genCouplers(size):
    """
    Generates a coupler list in 'couplers'+size based on the input size.
    """
    data=open('couplers'+str(int(size))+'.txt','w')

    for i in range(size*size*8):
        #if i==16:
        #    break
        if i%8<4:
            if i/(size*8)>0:
                #data.write(str(i-(size*8))+" "+str(i)+'\n')
                data.write(str(i-(size*8))+" "+str(i)+'\r\n')
        else:
            if (i%(size*8))/8>0:
                #data.write(str(i-8)+" "+str(i)+'\n')
                data.write(str(i-8)+" "+str(i)+'\r\n')
            for j in [0,1,2,3]:
                #data.write(str((i-i%8)+j)+" "+str(i)+'\n')
                data.write(str((i-i%8)+j)+" "+str(i)+'\r\n')
    data.close()