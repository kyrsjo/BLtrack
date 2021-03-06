import numpy as np
import sys
import os

import util
import Machine
import Beam


inputFile = None
if len(sys.argv) != 2:
    print "Usage: BLtrack.py inputfile"
    exit(1)
inputFile = sys.argv[1]
print "Got inputfile '"+inputFile+"'"

#Parse the input file; SixTrack-like input
myRing=None
myBeam=None
nTurns=None
ifile = open(inputFile,'r')
myBuffer = None
bufferStatus=None
for line in ifile:
    #print "LINE: '"+line.strip()+"'"
    #print "bufferStatus=",bufferStatus
    if line[0]== "/": # Comment
        continue
    
    #Tokenize
    lsp = line.split()
    
    if bufferStatus:
        if lsp[0]=="NEXT":
            #Construct ring/beam/... depending on the contents of BufferStatus
            if bufferStatus=="RING":
                myRing=Machine.Ring(myBuffer)
            elif bufferStatus=="BEAM":
                myBeam=Beam.Beam(myBuffer)
            else:
                print "Unknown bufferStatus '"+str(bufferStatus)+"'"

            myBuffer=None
            bufferStatus=None
            continue
        elif lsp[0] == "IMPORTFILE":
            assert len(lsp) == 2
            ifname2=os.path.join(os.path.dirname(inputFile), lsp[1])

            ifile2 = open(ifname2,'r')
            for line2 in ifile2:
                if line2[0]=="/":
                    continue
                assert line2[:4] != "NEXT"
                myBuffer+=line2
            ifile2.close()
        else:
            myBuffer+=line #Should include the newline...
        continue

    if lsp[0]=="RING":
        bufferStatus="RING"
        myBuffer=""
        continue
    elif lsp[0]=="BEAM":
        bufferStatus="BEAM"
        myBuffer=""
        continue
    elif lsp[0]=="TURNS":
        nTurns = int(lsp[1])
        continue
    elif lsp[0]=="END":
        break
    else:
        print "Error while parsing file; line not understood:"
        print line
        exit(1)
ifile.close()

#Sanity check input
if myRing==None:
    print "ERROR: No ring was set up"
    exit(1)
if myBeam==None:
    print "No beam was set up"
    exit(1)
print

print "Ring:"
print myRing
print "Total matrix:"
totalRE = myRing.getTotalMatrix()
print util.prettyPrint66(totalRE)
print "Eigenvalues and vectors:"
(w, v) = np.linalg.eig(totalRE)
for i in xrange(6):
    print "#%i : abs= %16.10g, angle= %16.10g [2pi] -> "%(i,np.absolute(w[i]), np.angle(w[i])/(2*np.pi)), w[i]
    #print "\t", v[i] #Complex vector, needs some better formatting...
print

print "Beam:"
print myBeam
print

#Track!
print "Tracking!"
print "# particles = ", myBeam.getNumParticles()
print "# elements  = ", myRing.getNumElements()
print "# turns     = ", nTurns
myRing.track(myBeam,nTurns)

#Postprocess? Close files etc.
