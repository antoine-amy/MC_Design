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
fm=[6e6, 56e6] #SB frequencies
R1 = 0.99; T1=np.sqrt(1-R1) # mirror 1
R2 = 0.99; T2=np.sqrt(1-R2) # mirror 2
n=1.1 # index of cavity
sidebands_trans_min=0.95 # Minimum transmission for the higher fm
TEM_limit=0.05 # Max transmission allowed for higher order TEM modes
Ls=[]; T=[]; TEM=[]
k=2*np.pi*f/c
nm_max=10
L_min=40e-3
F_min=5 # Minimum finesse of the cavity

print("-----Parameters intervals-----")
## Losses limitations
losses_max=1/100 #Maximum finesse, from losses limits
P=10/1000000
F_max=np.pi*losses_max/(4*P)
print(round(F_min,4),"<F<",round(F_max,4))

## Sidebands transmission
Fl_max=np.sqrt(((1-sidebands_trans_min)*c**2)/(16*sidebands_trans_min*fm[1]**2))
print(Fl_max)
Fl_min=L_min*F_min
L_max=Fl_max/F_min
print(round(L_min,4),"<L<",round(L_max,4))

## Filtering of higher modes
rl_max=TEM_limit*(2*Fl_max/(n*np.pi))**2
r_min=rl_max/L_max; r_max=rl_max/L_min
g_min=np.cos(1/n)**2 # From the conditions that the max n+m TEM should not go above the FSR
print(round(g_min,4),"<g<",round(1-L_min/r_max,4))


F=10
# Determine the best g to attenuate the higher order TEM
maxbefore=1
for i in np.arange(g_min,0.99,0.01):
  TEM=[]
  for j in range(1,nm_max+1):
    TEM.append(Tnm(F,i,j))
  if np.max(TEM)<maxbefore:
    maxbefore=np.max(TEM)
    g_final=i
g=g_final
print("Optimal g=",g)
TEM=[]
for i in range(nm_max+1):
    TEM.append(Tnm(F,g,i))
print()


# Carrier transmission, find the possibilities between L_min and L_max
q=np.ceil((2*L_min/carrier)-(1/2)-(np.arccos(np.sqrt(g))/np.pi)) # find minimum value
qmin=q
L_poss=L00(carrier,q,g)
while L_poss<L_max: # loop all the possibilities until L_max
  Ls.append(L_poss)
  q+=1
  L_poss=L00(carrier,q,g)

L=Ls[len(Ls)-1] #We keep the bigger one for now
print("-----Length determination-----")
print("Number of posibilities: ",  len(Ls))
print("L (mm)=", L*10**3)
print()

#Transmissions
print("-----Transmissions-----")
for i in range(len(fm)):
    T.append(Airy(F,L,fm[i]))
print("Transmission of the SB:",T)
print("Transmission of the TEM modes:",TEM)
print()

# Frequency of the TEM modes
print("-----Other stuff-----")
df=np.arccos(np.sqrt(g))*c/(2*np.pi*L)
print("df (MHz)=", df/1e6)

f=0 #a commenter pour center sur f

#Plotting transmission(frequency)
fsr = c/(2*n*L) # free spectral range
limits=[df*nm_max,np.max(fm),0.5*fsr]
f_range = np.linspace(f-np.max(limits), f+np.max(limits), 10000) # frequency range prbleme si on centre sur f?
print("FSR (MHz)=",fsr*10**-6)
plt.grid()
plt.plot(f_range/1e6, Airy(F,f_range,L), label="Cavity transmission")
plt.ylim(0,1.2)
plt.plot(f/1e6, 1, '.', color='black', label="Carrier")

for i in range(len(fm)):
  plt.arrow(f/1e6-fm[i]/1e6, 0, 0, T [i],head_width=10, head_length=0.03, color='red', alpha=1-i/len(fm), length_includes_head=True, label="SB ("+str(fm[i]/1e6)+"MHz)")
  plt.arrow(f/1e6+fm[i]/1e6, 0, 0, T [i],head_width=10, head_length=0.03,color='red', alpha=1-i/len(fm), length_includes_head=True)

plt.arrow(f/1e6+df/1e6, 0, 0, TEM[1],head_width=10, head_length=0.03, color='green',length_includes_head=True, label="TEM modes")
for i in range(1,nm_max+1):
  plt.arrow(f/1e6+i*df/1e6, 0, 0, TEM[i],head_width=10, head_length=0.03, color='green',length_includes_head=True)

#print(f[np.where(np.max(I))]-c/carrier)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Transmission')
plt.title('Transmission of a Mode Cleaner Cavity')
plt.legend()
plt.show()



# Rajouter calcul du waist


