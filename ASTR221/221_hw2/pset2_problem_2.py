'''
garcia, gil
ASTR221 HW
2/7/2019




In this program, we fit a blackbody curve to given data
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


#plotting our data
mag = np.array([9.31, 8.94, 8.11, 7.93, 7.84]) #b, v, j, h, k mags
f_v0 = np.array([4000,3600,1570,1020,636]) #flux of star in jansky
wavlen = np.array([445,551,1220,1630,2190]) #wavelength of each band respectively
wavlen = wavlen * 10**(-9) #convert nm to m
f_v0 = f_v0 * (10**(mag/-2.5)) #converting our maggies into flux values


plt.plot(wavlen,f_v0,'k.',label='Data Points') #plotting our data
plt.xlabel('Wavelength [m]') #labeling x axis
plt.ylabel(r'Flux $[Wm^{-2}Hz^{-1}]$') #labeling y axis

#defining the black body curve function
def black_body(lam,t):
    #defining constants
    c=2.9979*10**8 #speed of light [m/s]
    k=1.3806*10**(-23) #Boltzmann constant [m^2kgs^-2K^-1]
    h=6.6261*10**(-34) #Planck constant [m^2kgs^-1]
    #the equation
    term1 = (2.*h*c**2.)/(lam**5.)
    exp = ((h*c)/(lam*k*t))
    bttm = np.exp(exp) - 1.
    term2 = (1./bttm)
    return term1 * term2

lams = np.linspace(10**(-20),0.0000025,1000) #creating an array for wavlengths
curve = [] #creating an empty list to append fluxes into
scalar = 1/(1.9e12) #scalar to line up our curve with our data
temp = 3700 # guess and check until we line up the peaks of our data
for i in lams:
    curve.append(scalar*black_body(i,temp)) #finding the intensity for every wavelength in the array

plt.plot(lams,curve,color='k',label='Blackbody Curve') #plotting our curve
plt.legend()
plt.title('Given Data vs Blackbody Curve')
plt.show()
