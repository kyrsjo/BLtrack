import numpy as np

class Beam:
    bunches    = None  #Bunch objects
    bunches_z0 = None #Time of z=0 for the bunch objects, relative to some particle on the design orbit
    
    def __init__(self):
        "Construct an empty Beam object"
    def __init__(self,initStr):
        "Construct a Beam object from an input file fragment"
        
        self.bunches    = []
        self.bunches_z0 = []
        
        print "BEAM"
        myBuffer = None
        bufferStatus=None
        for line in initStr.splitlines():
            print "l=",line
            if line[:5] == "BUNCH":
                l = line[5:].split()
                assert len(l)==8, \
                    "Expected format: BUNCH t Nparticles sigx sigxp sigy sigyp sigz sigzp"
                self.bunches_z0.append(float(l[0]))
                N=int(l[1])
                lf = map(float,l[2:])
                self.bunches.append(Bunch(N, lf[0],lf[1],lf[2],lf[3],lf[4],lf[5]))
                #print self.bunches[-1].particles
                
            else:
                print "Error in Beam::__init__(initStr) while parsing line '"+line+"'"
                exit(1)
            
    def sortBunches(self):
        pass
        

class Bunch:
    "A group of particles which are assumed to be independent, i.e. no short-range wake"
    
    #6xN matrix of all particles in the bunch.
    particles = None
    N = None
    
    def __init__(self,N):
        "Construct an empty Bunch object"
        self.N=N
        self.particles=np.zeros((6,N))
    def __init__(self,N,sigx,sigxp,sigy,sigyp,sigz,sigE):
        self.N = N
        #self.particles=np.empty((6,N))
        self.makeGaussian(sigx,sigxp,sigy,sigyp,sigz,sigE)
        print self.particles
    def makeGaussian(self,sigx,sigxp,sigy,sigyp,sigz,sigE):
        mean=[0.0]*6
        cov=np.diag(v=(sigx,sigxp,sigy,sigyp,sigz,sigE))
        print mean
        print cov
        self.particles=np.random.multivariate_normal(mean,cov,self.N).transpose()
    
    def getMeans(self):
        pass
    def getSigmaMatrix(self):
        pass
    
class Ring:
    elements=None

    def __init__(self):
        "Construct an empty Ring object"
    def __init__(self,initStr):
        "Construct a Ring object from an input file fragment"
        print "RING"
        self.elements=[]
        
        for line in initStr.splitlines():
            print "line=",line
            if line[:6] == "MATRIX":
                l = line[6:].split()
                print l, len(l)
                assert len(l)==6*6
                lf = map(float,l)
                self.elements.append(SectorMapMatrix(lf))
            elif line[:9]=="PRINTMEAN":
                l = line.split()
                if len(l) == 1:
                    self.elements.append(PrintMean(None))
                elif len(l) == 2:
                    self.elements.append(PrintMean(l[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTMEAN while parsing line '"+line+"'"
                    exit(1)
                    
            elif line[:10]=="PRINTBUNCH":
                self.elements.append(PrintBunch())
            else:
                print "Error in Ring::__init__(initStr) while parsing line '"+line+"'"
                exit(1)
    def track(self,beam,turns):
        for i in xrange(turns):
            for element in self.elements:
                for bunch in beam.bunches:
                    bunch.particles=element.track(bunch)
    
class Element:
    "Base class for all elements"
    
    def track(self,bunch):
        return bunch # But it should be modified in-place

class SectorMapMatrix(Element):
    def __init__(self,RE):
        matrix=np.asarray(RE)
        self.matrix=np.reshape(matrix,(6,6))
        print self.matrix
    def track(self,bunch):
        return np.dot(self.matrix,bunch.particles)

class SectorMapTensor(Element):
    def __init__(self):
        pass
    def __init__(self,initStr):
        "Construct a SectorMapTensor object from an input file fragment"

class RFCavity(Element):
    def __init__(self):
        pass
    def __init__(self,initStr):
        "Construct a RFCavity object from an input file fragment"
class CrabCavity(Element):
    def __init__(self):
        pass
    def __init__(self,initStr):
        "Construct a CrabCavity object from an input file fragment"
class DumpParticles(Element):
    def __init__(self):
        pass
class PrintMean(Element):
    fname=None
    ofile=None
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
    def track(self,bunch):
        if self.ofile:
            self.ofile.write("%g %g %g %g %g %g\n" % tuple(np.mean(bunch.particles,axis=1)))
        else:
            print np.mean(bunch.particles,axis=1)
        return bunch.particles
class PrintBunch(Element):
    fname=None
    ofile=None
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
    def track(self,bunch):
        print
        print bunch.particles
        return bunch.particles


if __name__=="__main__":
    import sys
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

        if bufferStatus:
            if line[:4]=="NEXT":
                #Construct ring/beam/... depending on the contents of BufferStatus
                if bufferStatus=="RING":
                    myRing=Ring(myBuffer)
                elif bufferStatus=="BEAM":
                    myBeam=Beam(myBuffer)
                else:
                    print "Unknown bufferStatus '"+str(bufferStatus)+"'"
                
                myBuffer=None
                bufferStatus=None
                continue
            myBuffer+=line #Should include the newline...
            continue
            
        if line[:4]=="RING":
            bufferStatus="RING"
            myBuffer=""
            continue
        elif line[:4]=="BEAM":
            bufferStatus="BEAM"
            myBuffer=""
            continue
        elif line[:5]=="TURNS":
            nTurns = int(line[5:])
            continue
        elif line[:3]=="END":
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
    print

    #Track!
    myRing.track(myBeam,nTurns)
    
    #Postprocess? Close files etc.
