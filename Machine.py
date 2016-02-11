import numpy as np

import util

class Ring:
    elements=None
    length=None # Total length of the machine [m]
    
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
            elif lsp[0]=="RFCAV":
                l=lsp[1:]
                assert len(l) == 3, \
                    "ERROR in Ring::__init__(initStr)::RFCAV while parsing line '"+line+"'"
                self.elements.append(RFCavity(float(l[0]),float(l[1]),float(l[2])))
            elif lsp[0]=="RFCAV_MATRIX":
                l=lsp[1:]
                if len(l) == 4:
                    self.elements.append(RFCavity_Matrix(float(l[0]),float(l[1]),float(l[2]),float(l[3])))
                elif len(l)== 5:
                    self.elements.append(RFCavity_Matrix(float(l[0]),float(l[1]),float(l[2]),float(l[3]),float(l[4])))
                else:
                    print "ERROR in Ring::__init__(initStr)::RFCAV_MATRIX while parsing line '"+line+"'"
            elif lsp[0]=="RFCAV_LOADING":
                l=lsp[1:]
                assert len(l) == 4, \
                    "ERROR in Ring::__init__(initStr)::RFCAV_LOADING while parsing line '"+line+"'"
                self.elements.append(RFCavity_loading(float(l[0]),float(l[1]),float(l[2]),float(l[3])))
                
            elif lsp[0]=="PRINTMEAN":
                if len(lsp) == 1:
                    self.elements.append(PrintMean(None))
                elif len(lsp) == 2:
                    self.elements.append(PrintMean(lsp[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTMEAN while parsing line '"+line+"'"
                    exit(1)
            elif lsp[0]=="PRINTBUNCH":
                if len(lsp) == 1:
                    self.elements.append(PrintBunch(None))
                elif len(lsp) == 2:
                    self.elements.append(PrintBunch(lsp[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTBUNCH while parsing line '"+line+"'"
                    exit(1)
            elif lsp[0]=="RINGLENGTH":
                self.length = float(lsp[1])
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
        return bunch.particles
    def __str__(self):
        ret = "RFCavity:\n"
        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad]\n\n" % (self.voltage,self.wavelength,self.phase)
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
        return np.dot(self.RE,bunch.particles)
    def getMatrix(self):
        return self.RE
    def __str__(self):
        ret = "RFCavity_Matrix:\n"
        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad], E0 = %10g[GeV], m0 = %10g[MeV/c^2]\n\n" % (self.voltage,self.wavelength,self.phase, self.E0/1e9, self.m0/1e6)
        ret += util.prettyPrint66(self.RE)
        return ret

class RFCavity_loading(Element):
    voltage    = None # [V]
    wavelength = None # [m]
    phase      = None # [rad]
    RQ         = None
    
    #self.omegaC = None # [radians/second]
    
    def __init__(self,voltage,wavelength,phase,RQ):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
        self.RQ         = RQ
        
    def track(self,bunch,turn):
        voltage = - self.RQ * 2*np.pi*299792458/self.wavelength * 1.60217662e-19*1e11/bunch.N * np.exp(1j*bunch.particles[4,:]*2*np.pi/self.wavelength)
        #print np.absolute(voltage), np.angle(voltage)*180/np.pi
        
        bunch.particles[5,:] += self.voltage*np.sin(2*np.pi * bunch.particles[4,:] / self.wavelength + self.phase) / bunch.beam.p0
        return bunch.particles
    def __str__(self):
        ret = "RFCavity_loading:\n"
        ret += " voltage = %10g[V], wavelength = %10g[m], phase = %10g[rad], R/Q = %10g\n\n" % (self.voltage,self.wavelength,self.phase,self.RQ)
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
            self.ofile.write("%g %g %g %g %g %g\n" % tuple(bunch.getMeans()))
        else:
            print turn,
            print bunch.getMeans()
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
            self.ofile.write("# turn ID X[m] PX[px/p0] Y[m] PY[py/p0] T[=-ct, m] PT[=DeltaE/(p_s*c)]\n")
    def track(self,bunch,turn):
        for i in xrange(bunch.particles.shape[1]):
            outstring = "%i %i %g %g %g %g %g %g" % ((turn,i)+tuple(bunch.particles[:,i]))
            if self.ofile!=None:
                self.ofile.write(outstring+"\n")
            else:
                print outstring
        return bunch.particles
    def __str__(self):
        ret = "PrintBunch:\n"
        ret += "fname= '"+str(self.fname)+"'\n\n"
        return ret
