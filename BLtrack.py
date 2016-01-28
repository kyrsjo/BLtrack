import numpy as np

def prettyPrint66(M):
    assert M.shape == (6,6)
    ret = ""
    for i in xrange(6):
        for j in xrange(6):
            ret += " %15.5g" %(M[i,j],)
        ret += "\n"
    return ret

class Beam:
    E0         = None #Beam energy
    bunches    = None #Bunch objects
    bunches_z0 = None #Time of z=0 for the bunch objects, relative to some particle on the design orbit
    
    def __init__(self):
        "Construct an empty Beam object"
    def __init__(self,initStr):
        "Construct a Beam object from an input file fragment"
        
        self.bunches    = []
        self.bunches_z0 = []
        
        #print "BEAM"
        
        myBuffer = None
        bufferStatus=None
        for line in initStr.splitlines():
            #print "l=",line
            if line[:5] == "BUNCH":
                l = line[5:].split()
                assert len(l)==8, \
                    "Expected format: BUNCH t Nparticles sigx sigxp sigy sigyp sigz sigzp"
                self.bunches_z0.append(float(l[0]))
                N=int(l[1])
                lf = map(float,l[2:])
                self.bunches.append(Bunch(N, lf[0],lf[1],lf[2],lf[3],lf[4],lf[5]))
                self.bunches[-1].E0 = self.E0
                #print self.bunches[-1].particles
                
            elif line[:6] == "ENERGY":
                l=line[6:].split()
                assert len(l)==1,\
                    "Expected format: 'ENERGY E0[eV]', got: "+str(len(l))+" "+str(l)
                assert self.E0==None
                self.E0 = float(l[0])
                
            else:
                print "Error in Beam::__init__(initStr) while parsing line '"+line+"'"
                exit(1)
        assert self.E0!=None
            
    def sortBunches(self):
        pass
    def getNumParticles(self):
        N = 0
        for bunch in self.bunches:
            N += bunch.N
        return N
        

class Bunch:
    "A group of particles which are assumed to be independent, i.e. no short-range wake"
    
    #6xN matrix of all particles in the bunch.
    # Variables [madX units]:
    # X [m], PX [px/p0], Y [m], PY [py/p0], T [=-ct, m], PT [=DeltaE/(p_s*c)]
    particles = None
    N         = None
    
    E0        = None
    m0        = 938.272e9 #[eV/c^2]
    
    def __init__(self,N):
        "Construct an empty Bunch object"
        self.N=N
        self.particles=np.zeros((6,N))
    def __init__(self,N,sigx,sigxp,sigy,sigyp,sigz,sigE):
        self.N = N
        #self.particles=np.empty((6,N))
        self.makeGaussian(sigx,sigxp,sigy,sigyp,sigz,sigE)
        #print self.particles
    def makeGaussian(self,sigx,sigxp,sigy,sigyp,sigz,sigE):
        mean=[0.0]*6
        cov=np.diag(v=(sigx,sigxp,sigy,sigyp,sigz,sigE))
        #print mean
        #print cov
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
        #print "RING"
        self.elements=[]
        
        for line in initStr.splitlines():
            #print "line=",line
            if line[:6] == "MATRIX":
                l = line[6:].split()
                #print l, len(l)
                assert len(l)==6*6
                lf = map(float,l)
                self.elements.append(SectorMapMatrix(lf))
            elif line[:5]=="RFCAV":
                l=line[5:].split()
                assert len(l) == 3, \
                    "ERROR in Ring::__init__(initStr)::RFCAV while parsing line '"+line+"'"
                self.elements.append(RFCavity(float(l[0]),float(l[1]),float(l[2])))
                
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
            if (i+1)%1000 == 0 or (i+1)==turns:
                print "turn %15i /%15i (%5.1f %% )" %(i+1,turns,float(i+1)/float(turns)*100)
            for element in self.elements:
                for bunch in beam.bunches:
                    bunch.particles=element.track(bunch,i)

    def getTotalMatrix(self):
        ret = np.eye(6)
        for element in self.elements:
            ret = np.dot(element.getMatrix(),ret)
        return ret

    def __str__(self):
        ret = ""
        for element in self.elements:
            ret += str(element) + "\n"
        return ret
    
    def getNumElements(self):
        return len(self.elements)
    
class Element:
    "Base class for all elements"
    
    def track(self,bunch,turn):
        return bunch # But it should be modified in-place
    def getMatrix(self):
        return np.eye(6)
    def __str__(self):
        return "ELEMENT WITHOUT __STR__"

class SectorMapMatrix(Element):
    RE = None #Matrix for tracking (sector- or one-turn-map, as it comes out of MadX)
    def __init__(self,RE):
        self.RE=np.reshape(np.asarray(RE),(6,6))
    def track(self,bunch,turn):
        return np.dot(self.RE,bunch.particles)
    def __str__(self):
        ret = "SectorMapMatrix:\n"
        ret += prettyPrint66(self.RE)
        return ret
    def getMatrix(self):
        return self.RE

# class SectorMapTensor(Element):
#     def __init__(self):
#         pass
#     def __init__(self,initStr):
#         "Construct a SectorMapTensor object from an input file fragment"

class RFCavity(Element):
    voltage=None
    wavelength=None
    phase=None
    def __init__(self,voltage,wavelength,phase):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
    def track(self,bunch,turn):
        bunch.particles[5,:] += self.voltage*np.sin(-bunch.particles[4,:]/(2*np.pi*self.wavelength) + self.phase) \
                                / np.sqrt(bunch.E0**2-bunch.m0**2)
        return bunch.particles
    def __str__(self):
        ret = "RFCavity:\n"
        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad]\n\n" % (self.voltage,self.wavelength,self.phase)
        return ret
class CrabCavity(Element):
    voltage    = None # Transverse voltage Vcc [V]
    wavelength = None # 
    phase      = None # Phase offset, radians
    HoV        = None # Horizontal('H') or Vertical('V')?
    
    Vlong      = None # longitudinal kick voltage = Vcc*omega/c [V/mm]
    
    def __init__(self,voltage,wavelength,phase,HorV):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
        assert HorV == "H" or HorV=="V"
        self.HoV = HoV
        
        self.Vlong      = self.voltage*2*np.pi/self.wavelength*1e3
        
    def track(self,bunch,turn):
        if self.HoV=='H':
            bunch.particles[1,:] -= self.voltage \
                                    * np.cos(-bunch.particles[4,:]/(2*np.pi*self.wavelength) + self.phase) \
                                    / np.sqrt(bunch.E0**2-bunch.m0**2)
        elif self.HoV=='V':
            bunch.particles[3,:] -= self.voltage \
                                    * np.cos(-bunch.particles[4,:]/(2*np.pi*self.wavelength) + self.phase) \
                                    / np.sqrt(bunch.E0**2-bunch.m0**2)
        else:
            print "WTF in CrabCavity.track()"
            exit(1)
        
        bunch.particles[5,:] += self.Vlong \
                                * np.cos(-bunch.particles[4,:]/(2*np.pi*self.wavelength) + self.phase) \
                                / np.sqrt(bunch.E0**2-bunch.m0**2)
        return bunch.particles
    def __str__(self):
        ret = "CrabCavity:\n\n"
#        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad]\n\n" % (self.voltage,self.wavelength,self.phase)
        return ret
        
class DumpParticles(Element):
    def __init__(self):
        pass
    def __str__(self):
        ret = "DumpParticles:\n\n"
        return ret

class PrintMean(Element):
    fname=None
    ofile=None
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
    def track(self,bunch,turn):
        if self.ofile:
            self.ofile.write("%i " % (turn,))
            self.ofile.write("%g %g %g %g %g %g\n" % tuple(np.mean(bunch.particles,axis=1)))
        else:
            print turn,
            print np.mean(bunch.particles,axis=1)
        return bunch.particles
    def __str__(self):
        ret = "PrintMean:\n"
        ret += "fname= '"+str(self.fname)+"'\n\n"
        return ret
class PrintBunch(Element):
    fname=None
    ofile=None
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
    def track(self,bunch,turn):
        print
        print bunch.particles
        return bunch.particles
    def __str__(self):
        ret = "PrintBunch:\n"
        ret += "fname= '"+str(self.fname)+"'\n\n"
        return ret

if __name__=="__main__":
    import sys
    import os
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
            elif line[:10] == "IMPORTFILE":
                ls = line.split()
                assert len(ls) == 2
                ifname2=os.path.join(os.path.dirname(inputFile), ls[1])
                
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

    print "Ring:"
    print myRing
    print "Total matrix:"
    totalRE = myRing.getTotalMatrix()
    print prettyPrint66(totalRE)
    print "Eigenvalues:"
    (w, v) = np.linalg.eig(totalRE)
    for i in xrange(6):
        print "#%i : abs= %16.10g, angle= %16.10g [2pi] -> "%(i,np.absolute(w[i]), np.angle(w[i])/(2*np.pi)), w[i]
    print

    #Track!
    print "Tracking!"
    print "# particles = ", myBeam.getNumParticles()
    print "# elements  = ", myRing.getNumElements()
    print "# turns     = ", nTurns
    myRing.track(myBeam,nTurns)
    
    #Postprocess? Close files etc.
