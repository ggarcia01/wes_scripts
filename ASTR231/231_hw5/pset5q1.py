#garcia, gil
#astr 231 hw 5 prob 1
#4/9/2020

'''
in this script we integrate the eqns of stellar structure to determine the mass
and radius of a white dwarf of central density 10^5 gm /cm^3 under non-relativistic effects
'''



import numpy as np
import matplotlib.pyplot as plt



DR = 100 # our step size in the integrations


####################
####################             #1
####################


#some constants that will be needed
g = 6.67e-11 #m^3 kg^-1 s^-2
solar_mass = 1.98e30 #kg
solar_radius = 6.963e8 #meters

h = 6.626e-34 #J s (planck constant)
m_e = 9.109e-31 #kg (mass of electron)
m_H = 1.67e-27 #kg (mass of hydrogen)
c = 2.99e8 # m/s (speed of light)



####################
####################             #1a
####################

#we use Newton's method to solve the coupled diff eqns

#first, we define the routines that will be needed in our integrations:

k_non_rel = (h**2 / (5*m_e)) * (3/(8*np.pi))**(2/3) * (1/(2*m_H))**(5/3) #K constant in non-relativistic case

def density(pressure): #density for non-rel case
    return (pressure / k_non_rel)**(3/5)

def mass_step(density,dr): #using mass conservation diff eq to calculate one step in mass
    return 4*np.pi*r**2*density*dr

def pressure_step(mass,radius,density,dr): #using hydrostatic diff eq to calculate one step in pressure
    return ((g*mass) / radius**2 ) * density * dr


# now we initialize our variables and set them equal to the initial conditions

r = 0 #start at the core, radius is 0
m = 0 #at the core, w a radius of 0, there is no mass enclosed
d = 10**8 #initial density in kg / m^3
p = k_non_rel * (d)**(5/3) #initial pressure at the core
#p = k_ultra_rel * (d)**(5/3)

#intializing lists to store all values at each intigration step

r_lst = []
m_lst =[]
d_lst = []
p_lst=[]


#now we write the integrator using the fact that at the outer boundar, pressure will be 0:
iter=0
while p > 0:
    #keep track of num of steps we take
    iter+=1

    #updating our lsts at each step
    r_lst += [r]
    m_lst += [m]
    d_lst += [d]
    p_lst += [p]

    #now we update the values:
    r = r + DR #updating our radius value by our radius step size
    d = density(p) #calculating the new pressure value
    dm = mass_step(d,DR) #calculating the change in mass
    m = dm + m #updating our mass value
    dp = pressure_step(m,r,d,DR) #calculating the change in pressure
    p = p - dp #updating our pressure value. pressure decreases as we increase radius.



#    print(r,m,d,p)
#    print()



#results of the integration:
print('# of steps taken:',iter)

print('final radius in meters:\t', format(r,"E") )
print('final mass in kg:\t', format(m,'E') )


print('final radius in r_sun:\t', r / solar_radius )
print('final mass in m_sun\t', m / solar_mass)



#plotting our results:
r_arr = np.array(r_lst) / solar_radius
m_arr = np.array(m_lst) / solar_mass
d_arr = np.array(d_lst)
p_arr = np.array(p_lst)

print('average density in kg per m^3:\t' ,  format(np.average(d_arr),'E')  )

plt.subplot(311)
plt.title(r'properties of a star with $\gamma = 4 / 3$')
plt.plot(r_arr,m_arr,color='k',label=r'M$_r$')
plt.ylabel(r'mass [M$_{\odot}$]')
#plt.axhline(1.44,ls='--',color='k',label=r'Chandrasekhar Limit: 1.44 M$_{\odot}$')
#plt.ylim(-0.01,1.6)
left,right = plt.xlim()
plt.xlim(0,right)
plt.grid()
#plt.legend()


plt.subplot(312)
plt.plot(r_arr,np.log10(p_arr),color='k',label='pressure')
plt.ylabel(r'log(pressure [kg m$^{-1}$ s$^{-2}$] )')
left,right = plt.xlim()
plt.xlim(0,right)
plt.grid()



plt.subplot(313)
plt.plot(r_arr,np.log10(d_arr),color='k',label='density')
plt.ylabel(r'log(density [kg / m$^{3}$] )')
left,right = plt.xlim()
plt.xlim(0,right)
plt.grid()



plt.xlabel(r'radius[R$_{\odot}$]')
#plt.tight_layout()
plt.show()
print(k_non_rel)
