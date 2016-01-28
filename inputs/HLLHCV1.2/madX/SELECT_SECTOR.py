import sys
assert len(sys.argv)==2
fname = sys.argv[1]


ifile=open(fname,'r')

ofile=open(fname+".RING",'w')
prevPoint = "START"
firstPoint = None
isMaxtrix = None
for line in ifile:
    if line[0]=="$" or line[0]=="@":
        continue
    l = line.split()
    if l[0]=="*":
        l = l[1:]
        isMatrix=False
    else:
        isMatrix=True
    
    
    if isMatrix:
        ofile.write("/ " + prevPoint + " -> " + l[0] + "\n")
        prevPoint=l[0]
        if firstPoint == None:
            firstPoint = prevPoint
        elif firstPoint==l[0]:
            #MadX repeats the lattice twice??
            break
        ofile.write("MATRIX")
    print l[0], l[1], "R=",
    theMatrix = []
    for i in xrange(8,44):
        print l[i],
        if isMatrix==True:
            theMatrix.append(str(l[i]))
            #ofile.write(" "+str(l[i]))
    if len(theMatrix)==6*6:
        #MadX transposes the R matrix from sectormap relative to
        # the RE from one-turn-map -> transpose it back for the file
        for i in xrange(6):
            for j in xrange(6):
                ofile.write(" "+theMatrix[i+j*6])
                
    #Prepare for the next output
    print
    if isMatrix==True:
        ofile.write("\n")
ifile.close()
ofile.close()
