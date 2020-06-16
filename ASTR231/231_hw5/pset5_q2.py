#garcia, gil
#astr 231 hw 5 prob 1
#4/9/2020

'''
in this script we integrate the eqns of stellar structure to determine the maximum mass
and radius of a white dwarf experiencing ultra-relativistic effects
'''



import numpy as np
import matplotlib.pyplot as plt



DR = 10 # our step size in the integrations


####################
####################             #2
####################


#some constants that will be needed
g = 6.67e-11 #m^3 kg^-1 s^-2
solar_mass = 1.98e30 #kg
solar_radius = 6.963e8 #meters

h = 6.626e-34 #J s (planck constant)
m_e = 9.109e31 #kg (mass of electron)
m_H = 1.67e-27 #kg (mass of hydrogen)
c = 2.99e8 # m/s (speed of light)



####################
####################             #2a
####################

#we use Newton's method to solve the coupled diff eqns

#first, we define the routines that will be needed in our integrations:

k_ultra_rel = ((h*c)/4) * (3/ (8*np.pi))**(1/3) * (1/ (2*m_H))**(4/3) #K constant in ultra-relativistic case
#print(format(k_ultra_rel,'E'))


def density(pressure): #density for ultra-rel case
    return (pressure / k_ultra_rel)**(3/4)

def mass_step(density,r,dr): #using mass conservation diff eq to calculate one step in mass
    return 4*np.pi*r**2*density*dr

def pressure_step(mass,radius,density,dr): #using hydrostatic diff eq to calculate one step in pressure
    return ((g*mass) / radius**2 ) * density * dr



#we define a fxn that will be doing the intergrating

def newton_intergration(initial_density):
    # now we initialize our variables and set them equal to the initial conditions

    r = 0 #start at the core, radius is 0
    m = 0 #at the core, w a radius of 0, there is no mass enclosed
    d = initial_density #initial density in kg / m^3. this is the variable we will pass the newton_intergration fcn.
    p = k_ultra_rel * (d)**(4/3) #initial pressure at the core

        #now we write the integrator using the fact that at the outer boundar, pressure will be 0:
    iter=0
    while p > 0:
        #keep track of num of steps we take
        iter+=1

        #now we update the values:
        r = r + DR #updating our radius value by our radius step size
        d = density(p) #calculating the new pressure value
        dm = mass_step(d,r,DR) #calculating the change in mass
        m = dm + m #updating our mass value
        dp = pressure_step(m,r,d,DR) #calculating the change in pressure
        p = p - dp #updating our pressure value. pressure decreases as we increase radius.
    return m,r


#the range of density values that we will use in the intergration
d_range = np.linspace(10**12,10**18,10000)


#initialize empty lists where we will add the radius and mass for each density from each iteration
mass_lst=[]
radius_lst=[]

#running through our newton_intergration for each value in our d_range
for d_value in d_range:
    mass,radius =newton_intergration(d_value)
    mass_lst+=[mass]
    radius_lst +=[radius]
    print(d_value,mass,radius)






#plotting our results:
r_arr = np.array(radius_lst) / solar_radius
m_arr = np.array(mass_lst) / solar_mass
d_arr = np.log10(np.array(d_range))



plt.subplot(211)
plt.title('mass and radius dependance on central density')
plt.plot(d_arr,m_arr,color='k')
plt.ylabel(r'mass [M$_{\odot}$]')
plt.axhline(1.4416,ls='--',color='k',label=r'Chandrasekhar Limit: 1.44 M$_{\odot}$')
#plt.axhline(1.44,ls='--',color='k',label=r'Chandrasekhar Limit: 1.44 M$_{\odot}$')
plt.ylim(0,1.6)
#left,right = plt.xlim()
#plt.xlim(0,right)
plt.grid()
plt.legend()



plt.subplot(212)
plt.plot(d_arr,r_arr,color='k')
plt.ylabel(r'radius [R$_{\odot}$]')
#plt.ylabel('log(pressure)')
#plt.yscale('log')
plt.grid()

#left,right = plt.xlim()
#plt.xlim(0,right)

plt.xlabel(r'log( central density [kg / m$^{3}$] )')
plt.tight_layout()
plt.show()
