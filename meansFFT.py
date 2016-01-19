import numpy as np
import matplotlib.pyplot as plt

(x,xp,y,yp,z,E) = np.loadtxt("means.dat",unpack=True)
T = np.arange(1,len(x)+1)

plt.figure()
plt.plot(T,x);
#plt.show()

F = np.fft.rfftfreq(len(T),1.0)

xFFT = np.fft.rfft(x)
xpFFT = np.fft.rfft(xp)
yFFT = np.fft.rfft(y)
ypFFT = np.fft.rfft(yp)

plt.figure()
plt.plot(F,np.abs(xFFT))
plt.plot(F,np.real(xFFT))
plt.plot(F,np.imag(xFFT))
plt.savefig("FFT1.png")

plt.figure()
plt.plot(F,np.abs(xFFT), label="x")
#plt.plot(F,np.abs(xpFFT), label="xp")
plt.plot(F,np.abs(yFFT), label="y")
#plt.plot(F,np.abs(ypFFT), label="yp")
plt.legend()
plt.savefig("FFT2.png")
plt.xlim(0.3,0.33)
plt.savefig("FFT2-2.png")

plt.figure()
plt.semilogy(F,np.abs(xFFT), label="x")
#plt.semilogy(F,np.abs(xpFFT), label="xp")
plt.semilogy(F,np.abs(yFFT), label="y")
#plt.semilogy(F,np.abs(ypFFT), label="yp")
plt.legend()
plt.savefig("FFT3.png")
plt.xlim(0.3,0.33)
plt.savefig("FFT3-2.png")


plt.show()
