#garcia, gil
#astr 231 - hw1 prob 3
#2/6/2020


#in this script, we take the given data and fit a blackbody curve to it in order
# to determine the temperature of the source

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

df = pd.read_csv('pset1_q3_data.csv')

#functions

'''
Planck function for creating blackbody curves.
Input temperature in kelvins and wavelengths in meters unless you change
the default parameters to fit your units.
h = planck constant, default units in J*s
c = speed of light, default units in micron/s
k = boltzmann constant, default units in J/K.
Returns intensity in W/m^2 unless default parameters changed.
'''
def planck(wavelength,temperature,scale):
    #default parameters / constants
    h=6.626e-34
    c=2.998e14
    k=1.3806e-23
    #function
    term1 = (2 * h * c**2) / wavelength**5
    exponential_term = (h*c) / (wavelength * k * temperature)
    term2_bttm = np.exp(exponential_term) - 1
    term2 = 1 / term2_bttm
    return term1 * term2 * scale


#constants
c = 2.99e14 #speed of light in meters per seconds

#data manipulation
    #using zero mag fluxes & mags to solve for flux in frequency
df['f_nu(Jy)'] = df['zero_mag_flux_density(Jy)']*10**(df['magnitude'] / (-2.5))
    # solving for flux in W/m^2
df['f_lambda(W/m2/micron)'] = (df['f_nu(Jy)'] * (10**-26) * c) / (df['eff_wav(micron)']**2)

'''
fitting our data w the plank function
'''
pp = 5000,8.5e-9 #an initial guess for the temperature and scale factor
popt,pcov = curve_fit(planck,df['eff_wav(micron)'],df['f_lambda(W/m2/micron)'],pp) #scipy curve fitting

print('the temperature[k] of star is: \t',popt[0])
print('the scale factor is: \t',popt[1])

#plotting our fit line and data points
wavelengths = np.linspace(0.5,2.3,1000)
intensities = planck(wavelengths,popt[0],popt[1])
plt.plot(np.log10(wavelengths),intensities,color='black',label=r'%0.2E $\cdot$ best fit blackbody curve'%popt[1],ls='-')
plt.scatter(np.log10(df['eff_wav(micron)']),df['f_lambda(W/m2/micron)'],s=8,color='black',label='data')
plt.axvline(x=np.log10(df['eff_wav(micron)'].iloc[1]),color='k',ls='--',label=r'$\lambda_{peak}$= %0.2f $\mu$m'%df['eff_wav(micron)'].iloc[1])
plt.xlabel(r'$\log_{10}(\lambda)[\mu$m]')
plt.ylabel(r'$B_{\lambda}(T)$ [W m$^{-2}\mu$m$^{-1}$]')
plt.legend()
plt.title('Blackbody Curve')
plt.ylim(0, 2.2e-14)
plt.show()
print('temperature using peak wavelength and wein displacement law:', 2898 / df['eff_wav(micron)'].iloc[1] )
