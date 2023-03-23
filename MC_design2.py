# import the required libraries
import matplotlib.pyplot as plt
import numpy as np
import sys


def Tnm(F,g,i): # Transmission of TEM modes
  return 1/(1+(2*F/np.pi)**2*np.sin(i*np.arccos(np.sqrt(g)))**2)

def Airy(F,f, L): # Transmission as a function of frequency
  return 1/(1+(2*F/np.pi)**2*np.sin(2*np.pi*f*L/c)**2)

def L00(lbda,q,g):
   return (lbda/2)*(q+1/2+np.arccos(np.sqrt(g))/np.pi)


# Constants
c=2.99792e8  # speed of light
carrier=1064e-9 # Carrier wavelength
f=c/carrier
k=2*np.pi*f/c
fm=[6e6, 56e6] #SB frequencies
R1 = 0.99; T1=np.sqrt(1-R1) # mirror 1
R2 = 0.99; T2=np.sqrt(1-R2) # mirror 2
n=1.45 # index of cavity
nm_max=10 # Max n+m order
SB_limit=0.90 # Minimum transmission for the higher fm
TEM_limit=0.90 # Max transmission allowed for higher order TEM modes
T=[[],[]]; Tbefore=[[1]*len(fm),[1]*nm_max]
L_lim=[5e-2, 30e-2] # Optical length limits
F_lim=[5, 50] # Finesse limits
precisions=[1, 1e-2, 1e-1] # Precisions on F, L (m), and r (m)


print("-----Parameters intervals-----")
## Losses limitations
losses_max=1/100; P=10/1000000
F_max=np.pi*losses_max/(4*P)
if F_max<F_lim[1]: F_lim[1]=F_max
print(round(F_lim[0],4),"<F<",round(F_lim[1],4))
print()

for F in np.arange(F_lim[0],F_lim[1],precisions[0]):
  FL_max=(c*np.sqrt(1-SB_limit))/(4*fm[len(fm)-1]*np.sqrt(SB_limit))
  if FL_max/F<L_lim[1]:
    L_max=FL_max/F
  else:
    L_max=L_lim[1]
  print (round(100-((((F_lim[1]-F_lim[0])/precisions[0])-F)/(F_lim[1]-F_lim[0])/precisions[0])*100,1),"%", end="\r")

  for L in np.arange(L_lim[0], L_max, precisions[1]):
    rL_max=TEM_limit*(2*F*L/(n*np.pi))**2
    r_max=rL_max/L

    for r in np.arange(50e-2,r_max, precisions[2]):
        g=1-L/r
        for i in range(len(fm)):
            T[0].append(Airy(F,L,fm[i])) # Transmissions of SB
        for i in range(1,nm_max+1):
            T[1].append(Tnm(F,g,i)) # Transmissions of TEM modes
        dT=(np.max(Tbefore[1])-np.max(T[1]))+(np.max(T[0])-np.max(Tbefore[0]))
        if dT>0:
           Parameters=[F,L,r]
           Tparameters=T
        Tbefore=T
        T=[[],[]]
print()
print("-----Parameters choosed-----")
print("Parameters: ", Parameters)
print("Transmissions: ", Tparameters)
print()
g=1-Parameters[1]/Parameters[2]

# Carrier transmission
print("-----Length determination-----")
q_min=np.floor((2*Parameters[1]/carrier)-(1/2)-(np.arccos(np.sqrt(g))/np.pi)) # find minimum value
print("q_min=",q_min)
L_min=L00(carrier,q_min,g); L_max=L00(carrier,q_min+1,g)
if np.abs(Parameters[1]-L_min)<np.abs(Parameters[1]-L_max):
  Parameters[1]=L_min
else:
  Parameters[1]=L_max
print("L (mm)=", L*10**3)
print()



# Frequency of the TEM modes
print("-----Other stuff-----")
df=np.arccos(np.sqrt(g))*c/(2*np.pi*L)
print("df (MHz)=", df/1e6)

f=0 #a commenter pour centrer sur f

#Plotting transmission(frequency)
fsr = c/(2*n*Parameters[1]) # free spectral range
limits=[df*nm_max,np.max(fm),0.5*fsr]
f_range = np.linspace(f-np.max(limits), f+np.max(limits), 10000) # frequency range prbleme si on centre sur f?
print("FSR (MHz)=",fsr*10**-6)
plt.grid()
plt.plot(f_range/1e6, Airy(Parameters[0],f_range,Parameters[1]), label="Cavity transmission")
plt.ylim(0,1.2)
plt.plot(f/1e6, 1, '.', color='black', label="Carrier")

for i in range(len(fm)):
  plt.arrow(f/1e6-fm[i]/1e6, 0, 0, Tparameters[0][i],head_width=10, head_length=0.03, color='red', alpha=1-i/len(fm), length_includes_head=True, label="SB ("+str(fm[i]/1e6)+"MHz)")
  plt.arrow(f/1e6+fm[i]/1e6, 0, 0, Tparameters[0][i],head_width=10, head_length=0.03,color='red', alpha=1-i/len(fm), length_includes_head=True)

plt.arrow(f/1e6+df/1e6, 0, 0, Tparameters[1][1],head_width=10, head_length=0.03, color='green',length_includes_head=True, label="TEM modes")
for i in range(1,len(Tparameters[1])):
  plt.arrow(f/1e6+i*df/1e6, 0, 0, Tparameters[1][i],head_width=10, head_length=0.03, color='green',length_includes_head=True)

#print(f[np.where(np.max(I))]-c/carrier)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Transmission')
plt.title('Transmission of a Mode Cleaner Cavity')
plt.legend()
plt.show()



# Rajouter calcul du waist


