# import the required libraries
import matplotlib.pyplot as plt
import numpy as np

# define the constants
c = 3e8  # speed of light
L = 1.  # length of the optical cavity
F=150 # Finesse

# calculate the free spectral range (fsr) of the cavity
fsr = c/(2*L)

# generate a range of frequencies for plotting
f = np.arange(-.4*fsr,.4*fsr,fsr/1e4)

# define the cavity's reflectivity and transmissivity coefficients
r1 = .99  # input mirror reflectivity
r2 = .98  # end mirror reflectivity
t1 = (1-r1**2)**.5  # input mirror transmissivity
t2 = (1-r2**2)**.5

# calculate the complex reflection coefficient of the cavity
R = (r1-(r1**2+t1**2)*r2*np.e**(2j*2*np.pi*f*L/c))/(1-r1*r2*np.e**(2j*2*np.pi*f*L/c))

# define the frequency of the sideband signal
fm = 56e6

# calculate the complex reflection coefficients with modulation sidebands
Rffm = (r1-(r1**2+t1**2)*r2*np.e**(2j*2*np.pi*(f+fm)*L/c))/(1-r1*r2*np.e**(2j*2*np.pi*(f+fm)*L/c))
Rfnfm = (r1-(r1**2+t1**2)*r2*np.e**(2j*2*np.pi*(f-fm)*L/c))/(1-r1*r2*np.e**(2j*2*np.pi*(f-fm)*L/c))
pdh = R*np.conjugate(Rffm)-np.conjugate(R)*Rfnfm
plt.subplot(311)
plt.plot(f/1e6,100*np.abs(R)**2,'#880088')
plt.axis([-.4*fsr/1e6,.4*fsr/1e6,0,102])
plt.ylabel('Reflected power (\%)')
plt.subplot(312)
plt.plot(f/1e6,180*np.angle(R)/np.pi,'#880088')
plt.axis([-.4*fsr/1e6,.4*fsr/1e6,-35,35])
plt.ylabel('Reflected phase (deg.)')
plt.subplot(313)
plt.plot(f/1e6,np.imag(pdh),'#880088')
plt.ylabel('PDH readout (arb.)')
plt.xlabel('$f-f_\mathrm{res}$ (MHz)')
plt.show()