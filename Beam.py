import numpy as np

class Beam:
    E0         = None      #Beam energy [eV]
    beta0      = None      #Beam beta factor
    gamma0     = None      #Beam gamma factor
    m0         = 938.272e6 #Beam particles rest mass [eV/c^2] (default: proton)
    p0         = None      #Beam momentum [eV/c]
    
    bunches    = None      #Bunch objects
    bunches_z0 = None      #Time of z=0 for the bunch objects, relative to some particle on the design orbit
    
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
                self.bunches[-1].beam = self
                #print self.bunches[-1].particles
                
            elif line[:6] == "ENERGY":
                l=line[6:].split()
                assert len(l)==1,\
                    "Expected format: 'ENERGY E0[eV]', got: "+str(len(l))+" "+str(l)
                assert self.E0==None
                self.E0     = float(l[0])
                self.gamma0 = self.E0/self.m0
                self.beta0  = np.sqrt((1.0-1.0/np.sqrt(self.gamma0))*(1.0+1.0/np.sqrt(self.gamma0)))
                self.p0     = np.sqrt((self.E0-self.m0)*(self.E0+self.m0))
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

    def __str__(self):
        ret = ""
        ret += "E0     = " +str(self.E0/1e9)+ " [GeV]\n"
        ret += "gamma0 = " +str(self.gamma0)+ "\n"
        ret += "beta0  = " +str(self.beta0) + "\n"
        ret += "p0     = " +str(self.p0/1e9)+ " [GeV/c]\n"
        return ret
        

class Bunch:
    "A group of particles which are assumed to be independent, i.e. no short-range wake"
    
    #6xN matrix of all particles in the bunch.
    # Variables [madX units]:
    # X [m], PX [px/p0], Y [m], PY [py/p0], T [=-ct, m], PT [=DeltaE/(p_s*c)]
    particles = None
    N         = None
    
    beam      = None #Pointer back to the beam object

    
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
        #mean[4]+=2.0*sigz
        #mean[5]-=2.0*sigE
        cov=np.diag(v=(sigx,sigxp,sigy,sigyp,sigz,sigE))
        #print mean
        #print cov
        self.particles=np.random.multivariate_normal(mean,cov,self.N).transpose()
    
    def getMeans(self):
        pass
    def getSigmaMatrix(self):
        pass
