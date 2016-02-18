import numpy as np
import util

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
            lsp = line.split()
            if lsp[0] == "BUNCH":
                l = lsp[1:]
                assert len(l)==8, \
                    "Expected format: BUNCH t Nparticles sigx sigxp sigy sigyp sigz sigzp"
                self.bunches_z0.append(float(l[0]))
                N=int(l[1])
                lf = map(float,l[2:])
                self.bunches.append(Bunch(N, lf[0],lf[1],lf[2],lf[3],lf[4],lf[5]))
                self.bunches[-1].beam = self
                #print self.bunches[-1].particles
            elif lsp[0] == "BUNCH_TWISS":
                l = lsp[1:]
                assert len(l)==12, \
                    "Expected format: 'BUNCH_TWISS T0 N Np X XP BETX Y YP BETY EMIT sigT sigEnorm'"+\
                    "Got "+str(len(l))+" parameters."
                T0   = float(l[0])  # Position of T=0 for the bunch [m]
                N    = int(l[1])    # Number of tracked particles
                Np   = float(l[2])  # Number of electron charges in the bunch
                X    = float(l[3])  # X-offset of the bunch [m]
                XP   = float(l[4])  # XP-offset of the bunch [1~rad]
                BETX = float(l[5])  # X beta-function [m]
                Y    = float(l[6])  # Y-offset of the bunch [m]
                YP   = float(l[7])  # YP-offset of the bunch [1~rad]
                BETY = float(l[8])  # Y beta-function [m]
                EMIT = float(l[9])  # Emittance
                sigT = float(l[10]) # Bunch length [m]
                sigEnorm = float(l[11]) # Energy spread deltaE/E
                
                self.bunches_z0.append(float(l[0]))                
                self.bunches.append(Bunch.newBunch_twiss(N, Np, X, XP, BETX, Y, YP, BETY, EMIT, sigT, sigEnorm, self))
                
            elif lsp[0] == "ENERGY":
                l=lsp[1:]
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
        ret += "N      = " +str(self.getNumParticles())+"\n"
        for (bunch,z) in zip(self.bunches,self.bunches_z0):
            ret += "Bunch: t = "+str(z)+"\n"
            ret += "\t Means = "+str(bunch)+"\n"
        return ret
        

class Bunch:
    "A group of particles which are assumed to be independent."
    
    #6xN matrix of all particles in the bunch.
    # Variables [madX units]:
    # X [m], PX [px/p0], Y [m], PY [py/p0], T [=-ct, m], PT [=DeltaE/(p_s*c)]
    particles = None
    N         = None
    chargeN   = None #Number of electron charges in total
    
    beam      = None #Pointer back to the beam object
    
    
    def __init__(self,N):
        "Construct an empty Bunch object"
        self.N=N
        #self.particles=np.zeros((6,N))

    @staticmethod
    def newBunch_twiss(N, Np, X, XP, BETX, Y, YP, BETY, EMIT, sigT, sigEnorm, beam):

        eps_g = EMIT / (beam.beta0*beam.gamma0) # Geometrical emittance
        sigX  = np.sqrt(eps_g*BETX)
        sigXP = np.sqrt(eps_g/BETX)#*5
        sigY  = np.sqrt(eps_g*BETY)
        sigYP = np.sqrt(eps_g/BETY)#*5
        
        bunch = Bunch(N)
        bunch.beam=beam
        
        mean = (X,XP,Y,YP,0.0,0.0)
        cov = np.diag(v=(sigX**2,sigXP**2,sigY**2,sigYP**2,sigT**2,(sigEnorm*beam.E0/beam.p0)**2))
        #print util.prettyPrint66(cov)
        bunch.particles=np.random.multivariate_normal(mean,cov,N).transpose()

        #print "Created bunch:"
        #print "BETX=",BETX #m
        #print "BETY=",BETY #m
        #print "eps_g=",eps_g
        #print "gamma0=",beam.gamma0
        #print "beta0=",beam.beta0
        #print bunch.getMeans()
        
        return bunch
        
    # def __init__(self,N,sigx,sigxp,sigy,sigyp,sigt,sigPT):
    #     self.N = N
    #     #self.particles=np.empty((6,N))
    #     self.makeGaussian(sigx,sigxp,sigy,sigyp,sigt,sigPT)
    #     #print self.particles
    def makeGaussian(self,sigx,sigxp,sigy,sigyp,sigt,sigPT):
        mean=[0.0]*6
        #mean[4]+=2.0*sigz
        #mean[5]-=2.0*sigE
        cov=np.diag(v=(sigx,sigxp,sigy,sigyp,sigt,sigPT))
        #print mean
        #print cov
        self.particles=np.random.multivariate_normal(mean,cov,self.N).transpose()
    
    def getMeans(self):
        return np.mean(self.particles,axis=1)
    def getSigmaMatrix(self):
        pass
    
    def __str__(self):
        return str(self.getMeans())
