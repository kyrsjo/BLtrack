import numpy as np

def prettyPrint66(M):
    assert M.shape == (6,6)
    ret = ""
    for i in xrange(6):
        for j in xrange(6):
            ret += " %15.5g" %(M[i,j],)
        ret += "\n"
    return ret
    
#Constants
c = 299792458.0    #Speed of light [m/s]
e = 1.60217662e-19 #Electron charge [C]
