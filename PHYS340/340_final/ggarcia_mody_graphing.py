import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

t,x,y,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,tot_energy,potential,kinetic = np.genfromtxt('molecular_dynamics_2d_test.out').T
i,velocities_x,velocities_y = np.genfromtxt('mo_dy_velocities.out').T


plt.figure(1)
plt.plot(t,kinetic,color='red',label='kinetic',linestyle='-')
plt.plot(t,potential,color='brown',label='potential',linestyle='-')
plt.plot(t,tot_energy,color='magenta',label='total',linestyle='-')
plt.title('Total Energy over Time')
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.legend()
plt.show()

'''
#plotting time vs the y component of position
plt.figure(2)
plt.plot(t,y,color='red',linestyle='-',label='1')
plt.plot(t,y1,color='brown',linestyle='-',label='2')
#plt.plot(t,y2,color = 'magenta',linestyle='-',label='3')
plt.legend(title='Particle')
plt.title('Y Projection of Particle Motion')
plt.xlabel('Time')
plt.ylabel('Y motion')
plt.show()


#plotting time vs the x Component of position
plt.figure(3)
plt.plot(t,x,color='red',linestyle='-',label='1')
plt.plot(t,x1,color='blue',linestyle='-',label='2')
#plt.plot(t,x2,color = 'green',linestyle='-',label='3')
plt.legend(title='Particle')
plt.title('X Projection of Particle Motion')
plt.xlabel('Time')
plt.ylabel('X motion')
plt.show()
'''

#x and y trajectory of a particle
plt.figure(4)
plt.plot(x,y,color='red',linestyle='-',label='1')
plt.plot(x1,y1,color='blue',linestyle='-',label='2')
plt.plot(x2,y2,color = 'green',linestyle='-',label='3')
plt.plot(x3,y3,color='orange',linestyle='-',label='4')
plt.plot(x4,y4,color='black',linestyle='-',label='5')
plt.plot(x5,y5,color = 'magenta',linestyle='-',label='6')
plt.plot(x6,y6,color='brown',linestyle='-',label='7')
plt.plot(x7,y7,color='maroon',linestyle='-',label='8')
plt.plot(x8,y8,color = 'yellow',linestyle='-',label='9')
plt.plot(x9,y9,color='pink',linestyle='-',label='10')
plt.legend(title='Particle')
plt.title('X and Y Motion of Particle')
plt.xlabel('X Motion')
plt.ylabel('Y motion')
plt.show()


#doing our fitting here, learned to do this with the help of Ben
N = 100

param_bounds=([0,0,-10], [5,10,10])
def f_x(v, m, E, shift):
    return((1 / np.sqrt(2 * np.pi * m * ((2 * E)/(3*N))) * np.e ** (-0.5 * m * (v-shift)**2 * (3 * N)/(2*E))))

p0 = 0.2, 3.5, 0.5
fit = curve_fit(f_x,i,velocities_x, p0,bounds=param_bounds)

m = fit[0][0]
E = fit[0][1]
shift = fit[0][2]

xs = []
therm = []
for v in np.arange(-10,20,0.1):
    therm += [1 / np.sqrt(2 * np.pi * m * ((2 * E)/(3*N))) * np.e ** (-0.5 * m * (v-shift)**2 * (3 * N)/(2*E))]
    xs += [v]

plt.plot(xs,therm,'-',color = 'k',label='gaussian fit')
plt.plot(i,velocities_x,'.',color = 'red')
plt.ylabel('Density')
plt.xlabel(r'$V_x$ Bin')
plt.title('Velocities Distrbution for x Component')
plt.xlim(-11,11)
plt.show()

param_bounds=([0,0,-10], [5,10,10])
def f_x(v, m, E, shift):
    return((1 / np.sqrt(2 * np.pi * m * ((2 * E)/(3*N))) * np.e ** (-0.5 * m * (v-shift)**2 * (3 * N)/(2*E))))
p0 = 0.2, 3.5, 0.5
fit = curve_fit(f_x,i,velocities_y, p0,bounds=param_bounds)

m = fit[0][0]
E = fit[0][1]
shift = fit[0][2]

xs = []
therm = []
for v in np.arange(-10,20,0.1):
    therm += [1 / np.sqrt(2 * np.pi * m * ((2 * E)/(3*N))) * np.e ** (-0.5 * m * (v-shift)**2 * (3 * N)/(2*E))]
    xs += [v]

plt.plot(xs,therm,color = 'k',linestyle='-',label='gaussian fit')
plt.plot(i,velocities_y,'.',color='red')
plt.ylabel('Density')
plt.xlabel(r'$V_y$ Bin')
plt.title('Velocities Distrbution for y Component')
plt.xlim(-11,11)
plt.show()
