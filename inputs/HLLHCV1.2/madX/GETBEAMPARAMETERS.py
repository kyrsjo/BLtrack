import sys
assert len(sys.argv)==3, "USAGE: python GETBEAMPARAMETERS.py tfs-file element"
fname = sys.argv[1]
assert fname.endswith(".tfs")
ifile=open(fname,'r')

header = {}

wasFound = False
X = None
Y = None
betaX = None
betaY = None

for line in ifile:
    #Look for header
    if line.startswith("*"): #Header!
        assert len(header) == 0
        ls = line.split()[1:]
        for i in xrange(len(ls)):
            header[ls[i]] = i
    elif line.startswith("@") or line.startswith("$"):
        continue
    else:
        assert len(header) > 0
        ls = line.split()
        if ls[0][1:-1] == sys.argv[2]:
            assert not wasFound
            wasFound=True
            print "Found element '" + sys.argv[2] + "':"
            print line
            X = ls[header["X"]]
            PX = ls[header["PX"]]
            betaX = ls[header["BETX"]]
            alfaX = ls[header["ALFX"]]
            
            Y = ls[header["Y"]]
            PY = ls[header["PY"]]
            betaY = ls[header["BETY"]]
            alfaY = ls[header["ALFY"]]
print
print "Extracted data:"
print "X     = "+X
print "PX    = "+PX
print "betaX = "+betaX
print "alfaX = "+alfaX
print
print "Y     = "+Y
print "PY    = "+PY
print "betaY = "+betaY
print "alfaY = "+alfaY

        
