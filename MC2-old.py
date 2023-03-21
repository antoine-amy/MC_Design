# import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import sys


def Tnm(F,g): # Transmission of TEM modes
  return 1/(1+(2*F/np.pi)**2*np.sin(i*np.arccos(np.sqrt(g)))**2)

def Airy(F,f, L): # Transmission of sidebands (F=finesse, f=frequency), as a function of frequency
  return 1/(1+(2*F/np.pi)**2*np.sin(2*np.pi*f*L/c)**2)

# Constants
c=2.99792e8  # speed of light
carrier=1064e-9 # Carrier wavelength
f=c/carrier
fm=[6e6, 56e6] #SB frequencies
sidebands_trans_min=0.6 # Minimum transmission for the higher fm
R1 = 0.99; T1=np.sqrt(1-R1) # mirror 1
R2 = 0.99; T2=np.sqrt(1-R2) # mirror 2
n=1.1 # index of cavity
TEM_limit=0.01 # Max transmission allowed for higher order TEM modes


## Limits for F
F_min=5 # Minimum finesse of the cavity for now. Will be higher when computing the required attenuation factor of higher modes
losses_max=1/100 #Maximum finesse, from losses limits
P=10/1000000
F_max=np.pi*losses_max/(4*P)


## Limits for R
g=0.9 #a calculer grace aux conditions provenants du taux dattenuation souhaite pour les modes superieurs


## Limits for L 
# Sidebands transmission
Fl_min=np.sqrt(((1-sidebands_trans_min)*c**2)/(16*sidebands_trans_min*fm[1]**2))
L_min=Fl_min/F_max; L_max=Fl_min/F_min
print("L min (mm)=", L_min*10**3)
print("L max (mm)=", L_max*10**3)

# Carrier transmission, find the possibilities between L_min and L_max
Ls=[]
q=np.ceil((2*L_min/carrier)-(1/2)-(np.arccos(np.sqrt(g))/np.pi)) # find minimum value
L_poss=(carrier/2)*(q+1/2+np.arccos(np.sqrt(g))/np.pi)
while L_poss<L_max: # loop all the possibilities until L_max
  Ls.append(L_poss)
  q+=1
  L_poss=(carrier/2)*(q+1/2+np.arccos(np.sqrt(g))/np.pi)
print("Number of posibilities: ",  len(Ls))
L=Ls[len(Ls)-1] #We keep the bigger one for now
print("L (mm)=", L*10**3)

T=[]
for i in range(len(fm)):
    T.append(1/(1+(2*F_min/np.pi)**2*np.sin(2*np.pi*fm[i]*L/c)**2))
print("Transmission of the SB:",T)


#Transmission higher order T modes (carrier and sidebands?) aller jusqua ordre m+n=10
TEM=[]
nm_max=10
rhol=200*10**(-6) # max determinable a partir des conditions sur la transmission des modes d'ordre sup
w0=(carrier/n*np.pi)*np.sqrt(rhol)
fsep=c/(2*np.pi*n)*(1/np.sqrt(rhol))
print(fsep/10**6)
for i in range(nm_max+1):
    TEM.append(1/(1+(2*F_max/np.pi)**2*np.sin(i*np.arccos(np.sqrt(g)))**2))
print("Transmission of the TEM modes:",TEM)

#Plotting transmission(frequency)
fsr = c/(2*n*L) # free spectral range
f_range = np.linspace(-1*fsr, +1*fsr, 10000) # frequency range prbleme si on centre sur f?
Tcav_carrier=1/(1+(2*F_min/np.pi)**2*np.sin(2*np.pi*f_range*L/c)**2)

k=2*np.pi*f/c
print("FSR (MHz)=",fsr*10**-6)
plt.plot(f_range/1e6, Airy(F_min,f_range,L))
plt.plot(0/1e6, 1, '.', color='black')
plt.plot((-6e6)/1e6, T [0], '.', color='black')
plt.plot((-56e6)/1e6, T [1], '.', color='black')
#print(f[np.where(np.max(I))]-c/carrier)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Transmitted Intensity')
plt.title('Transmitted Intensity of a Fabry-Perot Cavity')
plt.show()




