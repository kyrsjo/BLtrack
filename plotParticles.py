import numpy as np
import matplotlib.pyplot as plt
import sys,os
assert len(sys.argv)==2, "USAGE: plotParticles.py particledumpfile"

particledumpfile = sys.argv[1]

# turn ID X[m] PX[px/p0] Y[m] PY[py/p0] T[=-ct, m] PT[=DeltaE/(p_s*c)]
(turn,ID,x,xp,y,yp,t,pt) = np.loadtxt(particledumpfile,unpack=True)

plotsdir = particledumpfile+"_plots"
if os.path.isdir(plotsdir):
    print "deleting", plotsdir
    for name in os.listdir(plotsdir):
        os.remove(os.path.join(plotsdir,name))
    os.rmdir(plotsdir)
print "creating", plotsdir
os.mkdir(plotsdir)

print "Getting min/max:"
xmin = x.min()
xmax = x.max()
xpmin = xp.min()
xpmax = xp.max()
ymin = y.min()
ymax = y.max()
ypmin = yp.min()
ypmax = yp.max()
tmin = t.min()
tmax = t.max()
ptmin = pt.min()
ptmax = pt.max()

for T in xrange(int(min(turn)),int(max(turn))):
#for T in xrange(0,int(max(turn)),10):
    print "plotting turn=",T
    
    plt.figure(1)
    plt.plot(x[turn==T],xp[turn==T],'.',label=str(T))
    plt.xlim(xmin,xmax)
    plt.ylim(xpmin,xpmax)
    plt.savefig(os.path.join(plotsdir,"xxp_%010i.png"%(T,)))
    plt.clf()
                
    plt.figure(3)
    plt.plot(t[turn==T],pt[turn==T],'.',label=str(T))
    plt.xlim(tmin,tmax)
    plt.ylim(ptmin,ptmax)
    plt.savefig(os.path.join(plotsdir,"tpt_%010i.png"%(T,)))
    plt.clf()

plt.figure(1)
plt.legend()
plt.savefig(os.path.join(plotsdir,"xxp_all.png"))
plt.figure(3)
plt.legend()
plt.savefig(os.path.join(plotsdir,"tpt_all.png"))
