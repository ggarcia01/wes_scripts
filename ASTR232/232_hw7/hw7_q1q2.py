'''
ASTR232
HW7 - questions 1 and 2
date: 12-11-19
garcia, gil
'''



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


######################
######################
######## #1a #########
######################
######################


file = '40lss.txt' #name of file that contains the data. data from SDSS.

#importing our file and adding headers to each row
columns = ['ra','z','g_mag','r_mag','i_mag'] #the column names for our file
df = pd.read_csv(file,delimiter='   ',index_col=False, header=None,names=columns)
print('file: ',file)
print('#1-----------------------------------------------')

    #calculating distance: d = z*c/h_o
c = 299792 #speed of light in km/s
h_o = 73 #hubble constant in km/s/Mpc
df['dist(Mpc)'] = (df['z'] * c) / h_o #calculating the distance in pc

#plotting the distribution of the galaxies in a polar plot
fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
c = ax.scatter(df['ra']*(np.pi/180),df['dist(Mpc)'], s=10,color='k',alpha=0.75,label='Galaxy')
ax.set_theta_offset(np.pi/4) # we are only given a quarter slice so we can cut the rest of the circle off
ax.set_thetamin(0)
ax.set_thetamax(90)
plt.title('Distribution of Galaxies our Sample')
plt.xlabel('dist[Mpc]')
plt.legend(markerscale=0.4)
plt.close()

#some squares more filled in than others. distribution of galaxies is not isotropic & homogeneous.
#this goes against the cosmological principle.


######################
######################
######## #1b #########
######################
######################

#setting up all the components required for the Bell equation

df['g-r'] = df['g_mag'] - df['r_mag'] #color
df['dist(pc)'] = df['dist(Mpc)'] * 10**6 #converting Mpc to ps
M_o = 4.75 #abs mag of sun
a_i = -0.222 #coefficient from table 7, bell et al, 2013
b_i = 0.864 #coefficient from table 7, bell et al, 2013
# i tot is the apparent magnitude of the galaxy

#computing the stellar mass using the Bell equation
df['stellar_mass(Msun)']  = 10**(0.4*( M_o - df['i_mag'] + (5*np.log10(df['dist(pc)'])) -5 + 2.5*(a_i+(b_i*df['g-r']))))
#computing total mass:
df['total_mass(Msun)'] = df['stellar_mass(Msun)'] / 0.15

#finding combined stellar mass of our wedge:
combined_stellar_mass = df['stellar_mass(Msun)'].sum()
print('combined stellar mass: \t', format(combined_stellar_mass,'E'))
#finding the combined total mass of our wedge:
total_combined_mass = df['total_mass(Msun)'].sum()
print('total combined (stellar+dark) mass: \t', format(total_combined_mass,'E'))
#calculating the volume of our wedge:
volume = 2.1*np.pi*(90/360)*(80**2 - 18**2)
print('volume(Mpc^3): \t',format(volume,'E')    )
#finding the mass density of our wedge
mass_density = total_combined_mass / volume
print('mass density of local universe now (solar masses/Mpc^3): \t ',format(mass_density,'E'))
#converting the mass density to units of kg/m^3
mass_of_sun = 1.988e30 #kg
cubic_mpc_to_cubic_m = 2.938e67
mass_density_kg_per_m_cubed = mass_density * mass_of_sun / cubic_mpc_to_cubic_m
print('mass density of local universe (kg/m^3) \t', format(mass_density_kg_per_m_cubed,'E'))

######################
######################
######## #1c #########
######################
######################

#the critical density fxn
def critical_density(H):
    G = 6.67e-11 #m^3 kg^-1 s^-2
    #g = 4.3e-3 #pc solar mass (km/s)^2
    return (3. * H**2) / (8.*np.pi*G)

#converting H0 to units of 1 / secs
Ho_in_one_over_secs = 73 / 3.086e19
#plugging in Ho to find critical density of the universe
critical_density = critical_density(Ho_in_one_over_secs)
print('critical density (kg/m^3): \t ',format(critical_density,'E'))
#fraction of calculated mass density over critical_density
print('mass density / critical density:\t ', format(mass_density_kg_per_m_cubed / critical_density,'E') )


##########
#dif plot for 1a - size corresponds to relative size
##########

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
c = ax.scatter(df['ra']*(np.pi/180), df['dist(Mpc)'], color='k', s=100*df['total_mass(Msun)']/df['total_mass(Msun)'].max(), alpha=0.4,label='Galaxy')
ax.set_thetamin(0)
ax.set_thetamax(90)
ax.set_theta_offset(np.pi/4)
plt.title('Distribution of Galaxies for our Sample')
plt.xlabel('dist[Mpc]')
plt.legend(markerscale=0.4)
plt.show()

######################
######################
######## #2a #########
######################
######################

print('#2a-----------------------------------------------')
#one over Ho is hubble time
hubble_time_sec = (1/(73/3.086e19))

'''
we are in a flat, single component (phantom dark energy) universe.so,
'''

# age of universe, eqn 1
def age_of_universe_flat_phantom_dark_energy(w,hubble_parameter=73):
    term1 = 1/(1+w)
    term2 = 1/(hubble_parameter / 3.086e19) #gives hubble time in secs
    return (2/3)*term1*term2 #returns to in seconds

w_phantom = -0.5

#calculating age of this universe
t_o_flat_phantom = age_of_universe_flat_phantom_dark_energy(w_phantom)
t_o_flat_phantom_yr = t_o_flat_phantom * 3.171e-8


#expansion law 1
def scale_factor_flat_phantom(t,w=w_phantom, t_o = t_o_flat_phantom):
    expo = 2/ (3*(1+w))
    return (t/t_o)**(expo)

#reporting age of this unvierse
print('t_o(secs) where universe is flat, phantom dark energy dominates:\t',format(t_o_flat_phantom,'E'))
print('t_o(yrs) where universe is flat, phantom dark energy dominates:\t',format(t_o_flat_phantom_yr,'E'))
print()

'''
in a flat, matter dominated universe:
'''

# age of universe, eqn 2
def age_of_universe_flat_matter(w,hubble_parameter=73):
    return (2/3) * 1/(hubble_parameter/3.086e19)

w_matter = 0

#calculating age of this universe:
t_o_flat_matter = age_of_universe_flat_matter(w_matter)
t_o_flat_matter_yr = t_o_flat_matter * 3.171e-8

#expansion law 2
def scale_factor_flat_matter(t,w=w_matter, t_o = t_o_flat_matter):
    return (t/t_o)**(2/3)

#reporting age of this universe:
print('t_o(secs) where universe is flat, matter dominates:\t', format(t_o_flat_matter,'E'))
print('t_o(yrs) where universe is flat, matter dominates:\t', format(t_o_flat_matter_yr,'E'))
print()

'''
in an empty, k = -1 (negative curvature)
'''

#age of this universe, eqn 3
def age_of_universe_neg_curve_empty(hubble_parameter=73):
    return 1/(hubble_parameter/3.086e19)

#calculating age of this universe:
t_o_neg_curve_empty = age_of_universe_neg_curve_empty()
t_o_neg_curve_empty_yr = t_o_neg_curve_empty * 3.171e-8

#expansion law 3
def scale_factor_neg_curve_empty(t,t_o= t_o_neg_curve_empty):
    return t/t_o

#reporting age of this universe:
print('t_o(secs) where universe is neg. curvature, empty :\t',format(t_o_neg_curve_empty,'E'))
print('t_o(yrs) where universe is neg. curvature, empty :\t',format(t_o_neg_curve_empty_yr,'E'))


#plotting a(t) vs Ho(t-to) for each universe model
plt.title('Expansion Law for Different Universe Models')
t = np.linspace(-1,4*t_o_flat_phantom,100000) #the time values we want to plug in
plt.ylim(0,7)
#phantom universe, using expansion law 1
plt.plot((1/hubble_time_sec)*(t-t_o_flat_phantom), scale_factor_flat_phantom(t), color='k',linestyle='--', label='flat, phantom dark energy dominated' )
#flat,  matter dominanted universe, using expansion law 2
plt.plot((1/hubble_time_sec)*(t-t_o_flat_matter), scale_factor_flat_matter(t), color='k',linestyle=':', label='flat,matter dominated' )
#negative curvature, empty universe, using expansion law 3
plt.plot((1/hubble_time_sec)*(t-t_o_neg_curve_empty), scale_factor_neg_curve_empty(t), color='k', label='negative curvature,empty' )
plt.legend()
plt.ylabel('a(t)')
plt.xlabel(r'$H_0 (t-t_0)$')
plt.show()

######################
######################
######## #2b #########
######################
######################
print('#2b-----------------------------------------------')
#we calculated this in 2a, we report them again here:
print('t_o(secs) where universe is flat, phantom dark energy dominates:\t',format(t_o_flat_phantom,'E'))
print('t_o(yrs) where universe is flat, phantom dark energy dominates:\t',format(t_o_flat_phantom_yr,'E'))
print()


######################
######################
######## #2c #########
######################
######################
print('#2c-----------------------------------------------')

#eqn for time emmited:
def time_emmited(z,w,t_o=t_o_flat_phantom):
    expon = 3/2 * (1+w)
    denom = (1+z)**expon
    return t_o / (denom)

#t_lookback = t_0 - t_e (using all values for phantom universe)
lookback = t_o_flat_phantom - time_emmited(4.3,-0.5)
#reporting our values:
print('lookback time (secs): \t', lookback)
print('lookback time (yrs):\t', format(lookback * 3.171e-8,'E'))
