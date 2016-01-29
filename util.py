import numpy as np

def prettyPrint66(M):
    assert M.shape == (6,6)
    ret = ""
    for i in xrange(6):
        for j in xrange(6):
            ret += " %15.5g" %(M[i,j],)
        ret += "\n"
    return ret
    
