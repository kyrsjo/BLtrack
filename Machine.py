import numpy as np
from time import time,ctime

import util

class Ring:
    elements = None
    length   = None # Total length of the machine [m]
    
    def __init__(self):
        "Construct an empty Ring object"
    def __init__(self,initStr):
        "Construct a Ring object from an input file fragment"
        #print "RING"
        self.elements=[]
        
        for line in initStr.splitlines():
            #print "line=",line
            lsp = line.split()
            
            if lsp[0] == "MATRIX":
                l = lsp[1:]
                #print l, len(l)
                assert len(l)==6*6
                lf = map(float,l)
                self.elements.append(SectorMapMatrix(lf))
                self.elements[-1].ring=self
            elif lsp[0]=="RFCAV":
                l=lsp[1:]
                assert len(l) == 3, \
                    "ERROR in Ring::__init__(initStr)::RFCAV while parsing line '"+line+"'"
                self.elements.append(RFCavity(float(l[0]),float(l[1]),float(l[2])))
                self.elements[-1].ring=self
            elif lsp[0]=="RFCAV_MATRIX":
                l=lsp[1:]
                if len(l) == 4:
                    self.elements.append(RFCavity_Matrix(float(l[0]),float(l[1]),float(l[2]),float(l[3])))
                elif len(l)== 5:
                    self.elements.append(RFCavity_Matrix(float(l[0]),float(l[1]),float(l[2]),float(l[3]),float(l[4])))
                else:
                    print "ERROR in Ring::__init__(initStr)::RFCAV_MATRIX while parsing line '"+line+"'"
                self.elements[-1].ring=self
            elif lsp[0]=="RFCAV_LOADING":
                l=lsp[1:]
                assert len(l) == 7, \
                    "ERROR in Ring::__init__(initStr)::RFCAV_LOADING while parsing line '"+line+"'"
                Vg         = float(l[0])
                wavelength = float(l[1])
                phase      = float(l[2])
                RQ         = float(l[3])
                QL         = float(l[4])
                mode       = l[5]
                fname      = l[6]
                
                self.elements.append(RFCavity_loading(Vg,wavelength,phase,RQ,QL,mode,fname))
                self.elements[-1].ring=self
            elif lsp[0]=="PRINTMEAN":
                if len(lsp) == 1:
                    self.elements.append(PrintMean(None))
                elif len(lsp) == 2:
                    self.elements.append(PrintMean(lsp[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTMEAN while parsing line '"+line+"'"
                    exit(1)
                self.elements[-1].ring=self
            elif lsp[0]=="PRINTBUNCH":
                if len(lsp) == 1:
                    self.elements.append(PrintBunch(None))
                elif len(lsp) == 2:
                    self.elements.append(PrintBunch(lsp[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTBUNCH while parsing line '"+line+"'"
                    exit(1)
                self.elements[-1].ring=self
            elif lsp[0]=="RINGLENGTH":
                self.length = float(lsp[1])
            else:
                print "Error in Ring::__init__(initStr) while parsing line '"+line+"'"
                exit(1)
    def track(self,beam,turns):
        t0 = time()
        for i in xrange(turns):
            if (i+1)%1000 == 0 or (i+1)==turns:
                t1 = time()
                tT = (t1-t0)/1000.0
                t0 = t1
                ETA = tT*(turns-i)
                ETAh = int(ETA/3600.0)
                ETAm = int((ETA-ETAh)/60.0)
                ETAs = ETA-ETAh*3600-ETAm*60.0
                print "turn %10i /%10i (%5.1f %% )" %(i+1,turns,float(i+1)/float(turns)*100), ";",
                print "Time/turn=%9.5f, remaining=%4i:%02i:%02i"% (tT, ETAh,ETAm,int(ETAs)), ";",
                print "ETA=", ctime(time()+ETA)
            for element in self.elements:
                for bunch in beam.bunches:
                    element.track(bunch,i)

        print "Tracking finished; now postProcessing:"
        for element in self.elements:
            element.postProcess()

    def getTotalMatrix(self):
        ret = np.eye(6)
        for element in self.elements:
            ret = np.dot(element.getMatrix(),ret)
        return ret

    def __str__(self):
        ret = ""
        for element in self.elements:
            ret += str(element) + "\n\n"
        return ret
    
    def getNumElements(self):
        return len(self.elements)
    
class Element:
    "Base class for all elements"
    
    ring = None #Pointer back to Ring object
    
    def track(self,bunch,turn):
        pass
    def getMatrix(self):
        return np.eye(6)
    def __str__(self):
        return "ELEMENT WITHOUT __STR__"
    def postProcess(self):
        "Optional postprocessing"
        pass

class SectorMapMatrix(Element):
    RE = None #Matrix for tracking (sector- or one-turn-map, as it comes out of MadX)
    def __init__(self,RE):
        self.RE=np.reshape(np.asarray(RE),(6,6))
    def track(self,bunch,turn):
        bunch.particles = np.dot(self.RE,bunch.particles)
    def __str__(self):
        ret = "SectorMapMatrix:\n"
        ret += util.prettyPrint66(self.RE)
        return ret
    def getMatrix(self):
        return self.RE

# class SectorMapTensor(Element):
#     def __init__(self):
#         pass
#     def __init__(self,initStr):
#         "Construct a SectorMapTensor object from an input file fragment"

class RFCavity(Element):
    voltage    = None # [V]
    wavelength = None # [m]
    phase      = None # [rad]
    
    def __init__(self,voltage,wavelength,phase):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
        
    def track(self,bunch,turn):
        bunch.particles[5,:] += self.voltage*np.sin(2*np.pi * bunch.particles[4,:] / self.wavelength + self.phase) / bunch.beam.p0
    def __str__(self):
        ret = "RFCavity: "
        ret += "voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad]" % (self.voltage,self.wavelength,self.phase)
        return ret

class RFCavity_Matrix(Element):
    voltage    = None # [V]
    wavelength = None # [m]
    phase      = None # [rad]
    
    E0 = None # Beam energy [eV]
    m0 = None # Beam particle mass [eV/c^2]
    p0 = None # Reference momentum
    
    RE = None
    
    def __init__(self,voltage,wavelength,phase, E0, m0=938.272e6):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
        
        self.E0 = E0
        self.m0 = m0
        self.p0 = np.sqrt((self.E0-self.m0)*(self.E0+self.m0))
        
        self.RE = np.eye(6)
        self.RE[5,4] = self.voltage*(2*np.pi/self.wavelength) / self.p0

    def track(self,bunch,turn):        
        bunch.particles=np.dot(self.RE,bunch.particles)
    def getMatrix(self):
        return self.RE
    def __str__(self):
        ret = "RFCavity_Matrix: "
        ret += "voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad], E0 = %10g[GeV], m0 = %10g[MeV/c^2]" % (self.voltage,self.wavelength,self.phase, self.E0/1e9, self.m0/1e6)
        ret += util.prettyPrint66(self.RE)
        return ret

class RFCavity_loading(Element):
    Vg         = None # Generator voltage [V]
    phase      = None # Generator phase   [rad]
    
    wavelength = None # [m]
    RQ         = None
    QL         = None
    
    mode       = None
    
    Vb         = None #Beam loading voltage from the previous bunch (complex)
    
    fname      = None
    of         = None #Output file pointer (optional)

    def __init__(self,Vg,wavelength,phase,RQ,QL,mode,fname=None):
        self.Vg         = Vg
        self.wavelength = wavelength
        self.phase      = phase
        self.RQ         = RQ
        self.QL         = QL
        
        self.mode       = mode
        # Accepted modes:
        # - 'simple':
        #   Keep Vg steady, apply Vb on top
        # - 'noload1':
        #   The beam loading has no effect on the beam (but it is still calculated)
        #
        assert self.mode=='simple' or self.mode=='noload1', "got mode='"+str(mode)+"'"
        
        if fname=="nofile":
            self.fname = None
        if fname!=None:
            self.of = open(fname,'w')
            self.of.write("# RFCavity_loading output file\n")
            self.of.write("# " +str(self)+"\n")
            self.of.write("#\n")
            self.of.write("#     Turns       |Vb|    ang(Vb)\n")
        self.Vb         = 0j
        
    def track(self,bunch,turn):
        #Calculate the current beam loading voltage at T=0:
        L = self.ring.length # Distance from previous bunch [m]
        assert type(L) == float, "ring length should be a float, is it defined?"
        self.Vb *= np.exp(-1j*2*np.pi*L/self.wavelength)*np.exp(- np.pi * L /(self.QL*self.wavelength))
        #Beam loading voltage for a single particle
        Vb0 = self.RQ * 2*np.pi*util.c/self.wavelength * util.e*1e11/bunch.N # TODO: Bunch charge hard-coded to 10^11...
        
        if self.mode == 'simple':
            #Sort the particles and iterate (this kills the performance in Python :( ):
            tKey = np.argsort(bunch.particles[4,:])
            for pind in tKey[::-1]: #Loop from (largest T -> smallest) to (smallest T -> largest t):
                bunch.particles[5,pind] += (  self.Vg*np.sin(2*np.pi * bunch.particles[4,pind] / self.wavelength + self.phase)
                                              - Vb0/2.0 
                                              + np.real(self.Vb*np.exp(1j*2*np.pi*bunch.particles[4,pind]/self.wavelength)) 
                                           ) / bunch.beam.p0
                self.Vb -=  Vb0*np.exp(-1j*2*np.pi*bunch.particles[4,pind]/self.wavelength)
        elif self.mode == 'noload1':
            bunch.particles[5,:] += self.Vg*np.sin(2*np.pi * bunch.particles[4,:] / self.wavelength + self.phase) / bunch.beam.p0
            self.Vb -=  np.sum(Vb0*np.exp(-1j*2*np.pi*bunch.particles[4,:]/self.wavelength))
        
        if self.of:
            self.of.write(" %10i %10g %10g \n" %(turn, np.absolute(self.Vb), np.angle(self.Vb)) )
        
    def __str__(self):
        ret = "RFCavity_loading: "
        ret += "Generator voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad], R/Q = %10g, QL = %10g, mode='%s'" % (self.Vg,self.wavelength,self.phase,self.RQ,self.QL,self.mode)
        return ret
    def postProcess(self):
        if self.of:
            self.of.close()

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
        
    def __str__(self):
        ret = "CrabCavity: "
#        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad]\n\n" % (self.voltage,self.wavelength,self.phase)
        return ret
        
class PrintMean(Element):
    fname=None
    ofile=None

    totalMean = None
    nTurns    = None
    
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
        self.totalMean = np.zeros(6)
        self.nTurns = 0
    def track(self,bunch,turn):
        means = bunch.getMeans()
        if self.ofile:
            self.ofile.write("%i " % (turn,))
            self.ofile.write("%g %g %g %g %g %g\n" % tuple(means))
        else:
            print turn,
            print means
        
        self.totalMean += means
        self.nTurns += 1
        
    def __str__(self):
        ret = "PrintMean: "
        ret += "fname= '"+str(self.fname)
        return ret

    def postProcess(self):
        if self.ofile:
            self.ofile.close()
        self.totalMean /= float(self.nTurns)
        print "Orbit estimate from PrintMean: %g %g %g %g %g %g" % tuple(self.totalMean)

class PrintBunch(Element):
    fname=None
    ofile=None
    def __init__(self,fname):
        if fname:
            self.fname=fname
            self.ofile=open(fname,'w')
            self.ofile.write("# turn ID X[m] PX[px/p0] Y[m] PY[py/p0] T[=-ct, m] PT[=DeltaE/(p_s*c)]\n")
    def track(self,bunch,turn):
        for i in xrange(bunch.particles.shape[1]):
            outstring = "%i %i %g %g %g %g %g %g" % ((turn,i)+tuple(bunch.particles[:,i]))
            if self.ofile!=None:
                self.ofile.write(outstring+"\n")
            else:
                print outstring
        
    def __str__(self):
        ret = "PrintBunch: "
        ret += "fname= '"+str(self.fname)
        return ret

    def postProcess(self):
        if self.ofile:
            self.ofile.close()
