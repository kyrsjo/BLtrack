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
                l = line.split()
                if len(l) == 1:
                    self.elements.append(PrintBunch(None))
                elif len(l) == 2:
                    self.elements.append(PrintBunch(l[1]))
                else:
                    print "Error in Ring::__init__(initStr)::PRINTBUNCH while parsing line '"+line+"'"
                    exit(1)
            elif line[:10]=="RINGLENGTH":
                self.length = float(line[10:])
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
    voltage=None
    wavelength=None
    phase=None
    def __init__(self,voltage,wavelength,phase):
        self.voltage    = voltage
        self.wavelength = wavelength
        self.phase      = phase
    def track(self,bunch,turn):
        bunch.particles[5,:] += self.voltage*np.sin(2*np.pi * bunch.particles[4,:] / self.wavelength + self.phase) / bunch.beam.p0
        #print np.sqrt(bunch.E0**2-bunch.m0**2)
        #print self.voltage*np.sin(-bunch.particles[4,:]/(2*np.pi*self.wavelength) + self.phase)
        #print bunch.particles[5,:]
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
