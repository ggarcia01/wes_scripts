'''
Gil Garcia
ASTR232 - Cosmology
Sept 11, 2019
HW1
'''

### routines needed for the assignment ###


#To find the line of best fit's slope and y-int, I created a method of least squares fxn
def least_squares(x,y):
	N=float(len(x))
	#We know will use the equations given in class to find A and B
	#A= y-int, B= slope
	delta=(N*(np.sum(x**2)))-((np.sum(x))**2)
	A_top=((np.sum(y))*(np.sum(x**2)))-((np.sum(x))*(np.sum(x*y)))
	A=(A_top)/(delta)
	B_top=((N*(np.sum(x*y)))-((np.sum(x))*(np.sum(y))))
	B=(B_top)/(delta)
	return A,B

### importing the required libraries ###
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#####################################
#####################################
#######          #2          ########
#####################################
#####################################

    #manually entering the data for table 1
period = [35.551341, 10.150730, 9.842425, 7.594904, 7.012877, 5.773380, 5.366270, 4.470916, 4.435462, 3.728190]
parallax = [2.01, 2.78, 3.14, 2.28, 3.00, 2.13, 3.66, 2.81, 1.90, 2.40]
v_mag = [3.652, 3.881, 3.731, 4.607, 4.526, 5.593, 3.950, 5.342, 5.632, 5.414]
a_v = [0.52,0.06,0.25,0.37,0.58,0.67,0.23,0.64,0.34,0.20]

	#converting the lists into a python dictionary
dict = {
'period(d)':period,
'parallax(mas)':parallax,
'v_mag':v_mag,
'a_v':a_v
}

	#converting the dictionary into a pandas dataframe for easier handling
df = pd.DataFrame.from_dict(dict)

#####################################
#####################################
#######          #2b         ########
#####################################
#####################################

    # creating a new column for the corrected magnitudes using  v_intrinsic = v_mag - a_v
df['v_intrinsic'] = df['v_mag'] - df['a_v']
    # using parallax (mas) to find distance in pc and making a new column for it
df['distance(pc)'] = 1000/ (df['parallax(mas)'])
    # calculating the abs mags using distance modulus and making a new column for it
df['abs_mag'] = -5*np.log10(df['distance(pc)']/10) + df['v_intrinsic']

#####################################
#####################################
#######          #2c         ########
#####################################
#####################################

	#plotting period against abs mag
plt.scatter( np.log10(df['period(d)'])-1, df['abs_mag'] ,color='k',s=12)
plt.xlabel(r'$  \log(P) -1 $')
plt.ylabel(r'$M_v$')
plt.close()

	#using the least_squares routine on our period and abs mag data to find the line of best fit
yint, slope = least_squares( np.log10(df['period(d)'])-1, df['abs_mag'] )

	#creating the line of best fit in order to plot it over our data
x = np.linspace( -0.5,0.7,100  )
plt.plot(x,(slope*x)+yint,'black',label=r'$M_v$ = '+str(yint)[:5]+' + '+str(slope)[:5]+'*log(P)-1'   )
plt.plot( np.log10(df['period(d)'])-1, df['abs_mag'] ,'.',label='Galactic Cepheids',c='#F5313F',marker = "*")
plt.xlabel(r'$  \log(P) -1 $')
plt.ylabel(r'$M_v$')
plt.legend()
plt.title('PLR Relation for Galactic Cepheids')
	#inverting the y axis so that brighter stars are higher in the plt
plt.gca().invert_yaxis()
	#saving plt as a jpg file
plt.savefig('2c.jpg')
plt.close()

#####################################
#####################################
#######          #3          ########
#####################################
#####################################

	#reading in table 2 into a pandas dataframe
table2 = pd.read_csv('hw1_table2.txt',sep='\s+',header=0,index_col=None)

#####################################
#####################################
#######          #3a         ########
#####################################
#####################################

	#creating a new column of a_v values found using a_v = 3.1 * E(B-V)
table2['A_v'] = 3.1 * table2['E(B-V)']
	#creating a new column of corrected v magnitudes
table2['v_intrinsic'] = table2['V_obs'] - table2['A_v']

#####################################
#####################################
#######          #3b         ########
#####################################
#####################################

	#creating a column of abs mag values using the PLR best fit line found in question 2
table2['abs_mag'] = yint + slope*(np.log10(table2['P(d)']) -1)
	#creating a column of distances using the distance modulus
table2['distance(pc)'] = 10** (( table2['v_intrinsic'] - table2['abs_mag'] +5 )/5)
	#converting the distances from pc to Mpc
table2['distance(Mpc)'] = table2['distance(pc)'] / 1e6
	#using pandas inherit fxns to find the mean and st.dev of our distances
print('avg distance to stars in Mpc: ',table2['distance(Mpc)'].mean())
print('std distance to stars in Mpc: ',table2['distance(Mpc)'].std())
