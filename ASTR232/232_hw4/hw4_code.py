'''
Gil Garcia
ASTR232 - Cosmology
HW4
10 - 11 - 19
'''

#importing required packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#importing csv with spectral line data for 6 galaxies
spectra_lines = pd.read_csv('spectra_lines.csv',index_col=False)
#importing line ratio data for HII galaxies
h2 = pd.read_csv('C:/Users/Gilberto/Documents/Wesleyan/senior_year/fall_term/ASTR232/hw4/hw4_zip/hw4_moodle/h2.txt',sep='\s+',index_col=False)

spectra_lines_latex = spectra_lines.to_latex(index=False)
#print(spectra_lines)

#computing line ratios for spectra_lines file:
spectra_lines['log([OIII]5007/h-beta)'] = np.log10(spectra_lines['[OIII]5007'] / spectra_lines['h-beta'])
spectra_lines['log([NII]6583/h-alpha)'] = np.log10(spectra_lines['[NII]6583']  / spectra_lines['h-alpha'])
spectra_lines['log([SII]6716 + [SII]6731/h-alpha)'] = np.log10( (spectra_lines['[SII]6716'] + spectra_lines['[SII]6731'] )  / spectra_lines['h-alpha'])
spectra_lines['log([OI]6300/h-alpha)'] = np.log10(spectra_lines['[OI]6300']  / spectra_lines['h-alpha'])

#convert the names column into a lst
spectra_lst = list(spectra_lines['galaxy_spectra'])
spectra_lst = [string.upper() for string in spectra_lst]

#plot 1

for i,spectra in enumerate(spectra_lst):
    x_coord = spectra_lines['log([NII]6583/h-alpha)'].iloc[i]
    y_coord = spectra_lines['log([OIII]5007/h-beta)'].iloc[i]
    plt.scatter(x_coord,y_coord,color='red',s=10)
    plt.text(x_coord+0.01,y_coord+0.01,spectra,fontsize=9)

plt.scatter(spectra_lines['log([NII]6583/h-alpha)'],spectra_lines['log([OIII]5007/h-beta)'],s=10,color='red',label='unknown galaxies')
plt.scatter(h2['N2/Ha'],h2['#O3/Hb'],color='k',s=2,alpha=0.5,label='H2 Galaxies')
plt.ylabel(r'$\log(\frac{[\mathrm{OIII}]5007}{H\beta})$')
plt.xlabel(r'$\log(\frac{[\mathrm{NII}]5007}{H\alpha})$')

x = np.linspace(-2.5,-0.1,1000)
plt.plot(x, [(((0.61 / (b-0.05))) +1.3) for b in x], linestyle = '--',color = 'blue',label='maximal starbust curve 1'  )

plt.legend()
plt.title('BPT Diagrams')
plt.savefig('bpt1.jpg')
plt.close()


#plot 2

for i,spectra in enumerate(spectra_lst):
    x_coord = spectra_lines['log([SII]6716 + [SII]6731/h-alpha)'].iloc[i]
    y_coord = spectra_lines['log([OIII]5007/h-beta)'].iloc[i]
    plt.scatter(x_coord,y_coord,color='red',s=10)
    plt.text(x_coord+0.01,y_coord+0.01,spectra,fontsize=9)

plt.scatter(spectra_lines['log([SII]6716 + [SII]6731/h-alpha)'],spectra_lines['log([OIII]5007/h-beta)'],s=10,color='red',label='unknown galaxies')
plt.scatter(h2['S2/Ha'],h2['#O3/Hb'],color='k',s=2,alpha=0.5,label='H2 Galaxies')
plt.ylabel(r'$\log(\frac{[\mathrm{OIII}]5007}{H\beta})$')
plt.xlabel(r'$\log(\frac{ [\mathrm{SII}]6716 + [\mathrm{SII}]6731 }{H\alpha})$')

x1 = np.linspace(-2.5,0.05,1000)
plt.plot(x1, [(((0.72 / (x-0.32))) +1.30) for x in x1], linestyle = '--',color = 'green',label='maximal starbust curve 2'  )

plt.legend()
#plt.title('BPT Diagram')
plt.savefig('bpt2.jpg')
plt.close()

#plot 3

for i,spectra in enumerate(spectra_lst):
    x_coord = spectra_lines['log([OI]6300/h-alpha)'].iloc[i]
    y_coord = spectra_lines['log([OIII]5007/h-beta)'].iloc[i]
    plt.scatter(x_coord,y_coord,color='red',s=10)
    plt.text(x_coord+0.01,y_coord+0.01,spectra,fontsize=9)

plt.scatter(spectra_lines['log([OI]6300/h-alpha)'],spectra_lines['log([OIII]5007/h-beta)'],s=10,color='red',label='unknown galaxies')
plt.scatter(h2['O1/Ha'],h2['#O3/Hb'],color='k',s=2,alpha=0.5,label='H2 Galaxies')
plt.xlabel(r'$\log(\frac{[\mathrm{OI}]6300}{H\alpha})$')
plt.ylabel(r'$\log(\frac{[\mathrm{OIII}]5007}{H\beta})$')


x2 = np.linspace(-3,-0.8,1000)
plt.plot(x2, [(((0.73 / (a+0.59))) +1.33) for a in x2], linestyle = '--',color = 'orange',label='maximal starbust curve 3'  )

plt.legend()
#plt.title('BPT Diagram')
plt.savefig('bpt3.jpg')
plt.close()


##############################
##############################
##############################
###########  1c ##############
##############################
##############################
##############################
##############################

### temperature ###

# to measure the 4959 + 5007 line, assume = 4/3 * 5007. Thus,
spectra_lines['(4959+5007)/4363'] = ((4/3)*spectra_lines['[OIII]5007'])/spectra_lines['[OIII]4363']
spectra_lines['T_electron(K)'] = [15200,14000,14000,14000,9000,14000]

### density ###
spectra_lines['S2(6716/6731)'] = spectra_lines['[SII]6716'] / spectra_lines['[SII]6731']


spectra_lines['electron_density(cm-3)_uncorrected'] = [180,180,100,700,700,350]
spectra_lines['electron_density(cm-3)_corrected'] = spectra_lines['electron_density(cm-3)_uncorrected'] * (10**4/spectra_lines['T_electron(K)'])**(0.5)

#print(spectra_lines)

spectra_lines_latex = spectra_lines[['galaxy_spectra','S2(6716/6731)','electron_density(cm-3)_uncorrected','electron_density(cm-3)_corrected']].to_latex(index=False)


##############################
##############################
##############################
###########  2  ##############
##############################
##############################
##############################
##############################


###luminosity calculation ###

names = ['fits_file','num1','num2','num3','distance','dist_unit']
distances = pd.read_csv('C:/Users/Gilberto/Documents/Wesleyan/senior_year/fall_term/ASTR232/hw4/hw4_zip/hw4_moodle/distances.txt',sep='\s+',index_col=False,skiprows=1,names=names)

distances['dist(cm)'] = distances['distance'] * 3.086e24 #convert Mpc to cm
spectra_lines['lum_h-beta(erg/s)'] = spectra_lines['h-beta']*10**(-17) * distances['dist(cm)']**2 * 4 * np.pi

### radius calculation ###

#some constants:
alpha_Hb = (2.59e-13)/8.5 # for temp = 10^4
h = 6.626e-27 #erg/sec
nu_Hb = (2.997e10) / (4.86e-5) #cm/sec
epsilon = 10**(-2)

#solving for radius:

spectra_lines['radius_nlr(cm)'] = ((3*spectra_lines['lum_h-beta(erg/s)']) / (4*np.pi*alpha_Hb*h*nu_Hb*epsilon*spectra_lines['electron_density(cm-3)_corrected']**2))**(1/3)
spectra_lines['radius_nlr(cm)'] = ( (3 * spectra_lines['lum_h-beta(erg/s)']) /(4*np.pi*alpha_Hb*h*nu_Hb*spectra_lines['electron_density(cm-3)_corrected']**2 *epsilon) )**(1/3)
spectra_lines['radius_nlr(pc)'] = spectra_lines['radius_nlr(cm)'] * 3.24e-19


### mass calculation ###

#more constants
m_proton = 1.67e-24 #grams

#solving for mass:

spectra_lines['mass_nlr(grams)'] = (m_proton * spectra_lines['lum_h-beta(erg/s)'] ) / (alpha_Hb*h*nu_Hb*spectra_lines['electron_density(cm-3)_corrected'])
spectra_lines['mass_nlr(solar_mass)'] = spectra_lines['mass_nlr(grams)'] / 1.99e33
print(spectra_lines)
