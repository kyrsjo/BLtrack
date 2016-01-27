import sys
assert len(sys.argv)==2
fname = sys.argv[1]


ifile=open(fname,'r')

ofile=open(fname+".RING",'w')
prevPoint = "START"
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
    
    print l[0], l[1], "R=",
    if isMatrix:
        ofile.write("/ " + prevPoint + " -> " + l[0] + "\n")
        prevPoint=l[0]
        ofile.write("MATRIX")
    for i in xrange(8,44):
        print l[i],
        if isMatrix==True:
            ofile.write(" "+str(l[i]))
    print
    if isMatrix==True:
        ofile.write("\n")
ifile.close()
ofile.close()
