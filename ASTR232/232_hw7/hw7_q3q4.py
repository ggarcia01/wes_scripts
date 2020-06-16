'''
ASTR232
HW7 - questions 3 and 4
date: 12-11-19
garcia, gil
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#creating our  simpson integrator fxn
def simp_intergration(function,a,b,steps):
    h = (b-a)/steps
    x = np.linspace(a,b,steps+1)
    f_x = function(x)
    integral_evaluation = h/3 * np.sum(f_x[0:-1:2] + 4*f_x[1::2] + f_x[2::2])
    return integral_evaluation


######################
######################
######## #3a #########
######################
######################

#defining some constants
Ho = 73 #hubble const in km/s/Mpc
c = 3e5 # speed of light in km/s
dh = c/Ho
omega_r,omega_m,omega_lam,omega_0 = 0,0.27,0.73,1

#proper distance integral
def func(z):
    dh = c/Ho
    a = (1+z)**(-1)
    denom = omega_r *a**-4 + omega_m*a**-3 + omega_lam + (1-omega_0)*a**-2
    return dh/np.sqrt(denom)

# using simpson on our proper distance integral to get proper distance
def dp(z,steps=2000,func=func):
    a=0
    return  simp_intergration(func,a,z,steps)

#using proper distance to compute luminosity distance
def lum_dist(z,steps=2000,func=func):
    return (1+z) * dp(z,steps,func)


#testing integrator on x^2
def x_squared(x):
    return x*x

value = simp_intergration(x_squared,0,7,700)
#print('value',value)
# it works!!



######################
######################
######## #3b #########
######################
######################

print('#3b-------------------------------------------')

#creating empty lists to append our values, to be plotted
lum_dist_lst =[]
lum_dist_lst_hubble =[]
z_lst = np.linspace(0,7,5000) #range of z values to plug in

counter = 0
for z in z_lst: #calculating luminosity distance for 5,000 values between 0 and 7 redshifts
    lum_dists = lum_dist(z) #lum dist calculation
    lum_dist_hubble = (c*z)/Ho #hubble approximation calculation
    lum_dist_lst += [lum_dists] #append lum dist to lum dist lst
    lum_dist_lst_hubble += [lum_dist_hubble] #append hubble approx to hubble approx lst
    ratio = (lum_dists - lum_dist_hubble) / lum_dist_hubble #finding ratio of distances
    if ratio > 0.1 and counter == 0: #prints out when it crosses 10% difference
        print('z at which they differ by more than 10%: ',z,'\t ratio: ',ratio)
        counter +=1

#converting lists to array - dividing by dh
lum_dist_array = np.array(lum_dist_lst)
lum_dist_array = lum_dist_array / dh


lum_dist_array_hubble = np.array(lum_dist_lst_hubble)
lum_dist_array_hubble = lum_dist_array_hubble / dh



######################
######################
######## #3c #########
######################
######################


#plotting the benchmark + hubble approx models
plt.title('Benchmark Model')
plt.plot(z_lst,lum_dist_array,c='k',linestyle=':',label='Benchmark')
plt.plot(z_lst,lum_dist_array_hubble,c='k',linestyle='--',label='Hubble Approx')
plt.axvline(0.13,color='red',label='z=0.13')
plt.xlabel('z')
plt.ylabel('dL / dH')
plt.legend()
plt.close()


######################
######################
######## #4a #########
######################
######################


print('#4a-------------------------------------------')

#constants
energy = [0.16e3,3.5e3] #eV
h = 4.135e-15 #plancks constant in eV * s
c = 2.998e8 #speed of light in m/s

#fxn for converting energy to wavelength
def energy_to_wavelength(energy): #enrgy must be in eV. returns wavelength in nm
    h = 4.135e-15 #plancks constant in eV * s
    c = 2.998e8 #speed of light in m/s
    return (10**9*c*h) / energy

#fxn for converting observed wavelength to intrinsic (emmited) wavelength
def intrinsic_lambda(z_observed,lambda_observed):
    return lambda_observed / (1+z_observed)



#computing the intrinsic wavelength for our 2 energy values
for val in energy:
    wav = energy_to_wavelength(val)
    print('intrinsitc wavelength for energy: ',val,'is \t',intrinsic_lambda(4.3,wav))
# our values give us hard x-ray band


######################
######################
######## #4b #########
######################
######################
print('#4b-------------------------------------------')

#luminisisty fxn
def luminosity(flux,distance): #distance in cm, flux in erg * cm^-2 * s^-1
    return 4*np.pi*distance**2*flux

#constants
Ho = 73 #hubble const in km/s/Mpc
c = 3e5 # speed of light in km/s
dh = c/Ho
flux = 6.4e-13 # erg cm^-2 s^-1
omega_r,omega_m,omega_lam,omega_0 = 0,0.27,0.73,1
z = 4.3

#calculating lum dist at a z of 4.3
lum_dist_4b = lum_dist(z,steps=2000,func=func)
print('z:',z,'\t lum dist(Mpc): \t',format(lum_dist_4b,'E'))
#converting to lum dist to cm
lum_dist_4b_cm = lum_dist_4b * 3.086e24
print('z:',z,'\t lum dist(cm): \t',lum_dist_4b_cm)

#computing luminosity of our object
print('luminosity(erg/s):\t',luminosity(flux,lum_dist_4b_cm))

######################
######################
######## #4c #########
######################
######################

print('#4c-------------------------------------------')

#integral fxn for finding time since light was emmited
def func_4c(z):
    Ho = 73/3.086e19
    denom = (1+z)*np.sqrt(omega_m*(1+z)**3 + omega_lam)
    return (1/Ho)/(denom)

#integrating fxn to find time for our object
#in theory, we would use an upper bound of infinity but we cant do that
#we use an upper bound of z = 20,000
time_emit = simp_intergration(func_4c,4.3,20000,200000)
time_emit_yr = time_emit * 3.171e-8
print('time emit(s):\t',format(time_emit,'E'))
print('time emit(yr):\t',format(time_emit_yr,'E'))
#comparing our answer to hubble time
hubble_time = 1/(2.37e-18)
print('hubble time: \t',hubble_time)
print('fraction of the universe current age: \t',time_emit/hubble_time)
