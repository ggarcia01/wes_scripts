# Gil Garcia
# ASTR221 - hw3 prob 2
#4/4/2019


'''
In this script, we query the Gaia DR2 catalog to find the HR diagram for M4, a globular cluster
and fit it with ZAMS
'''

# we import the required libraries
import numpy as np
import matplotlib.pyplot as plt
# we use vizier to query the Hiparcos catalog
from astroquery.vizier import Vizier
# we use angle to define our radius for our cone searches
from astropy.coordinates import Angle
# we will query gaia dr2
from astroquery.gaia import Gaia
#data management
import pandas as pd

'''
defini
'''

#The distance modulus to calculate the absolute magnitude
def distance_modulus(apparant_m,distance):
    term2 = 5*(np.log10(distance)-1)
    return (apparant_m - term2)

#read in file of isochrones data
df = pd.read_csv('minus1point07age12.csv')

'''
to find tha abs value for the G,BP,RP mags, we need to find the distance to
according to Gaia. Thus, for m = 1.2847 and M=0, we have
'''
dist_vega = 10**((0.2*5)+(0.2*1.3847))


#querying the gaia dr2 catalog
job = Gaia.launch_job_async("SELECT \
pmra, pmdec, parallax, parallax_error, parallax_over_error, phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_mag, \
bp_rp, bp_g,g_rp \
FROM gaiadr2.gaia_source \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec), CIRCLE('ICRS',245.8958,-26.5256,1))=1 \
AND parallax_error < 2 \
AND parallax_over_error > 5 \
AND pmra BETWEEN -14.5 and -10.5 \
AND pmdec BETWEEN -21 and -17 \
AND parallax BETWEEN 0.35 and 0.48")
gaia_results = job.get_results()


# creating proper motion plt
plt.figure(1)
plt.plot(gaia_results['pmra'],gaia_results['pmdec'],'.',color = 'maroon',alpha = 0.3)
plt.xlabel('Proper Motion in RA')
plt.ylabel('Proper Motion in DEC')
plt.title('Proper Motion of M4 Using Gaia Catalog')
plt.grid()
plt.savefig('gaia_m4_proper_motion.png')
#plt.show()
plt.close()
g_distance = 1/ (gaia_results['parallax']/1000.)
abs_g_mags = distance_modulus(gaia_results['phot_g_mean_mag'],g_distance)

# creating an HR diagram
plt.figure(2)
plt.plot(gaia_results['bp_rp'],abs_g_mags,'.',color = 'black',alpha=0.6,markersize = 1.5)
plt.gca().invert_yaxis()
plt.xlabel('BP - RP mag')
plt.ylabel('G mag')
plt.title('HR Diagram Using M4 in Gaia Catalog')
plt.savefig('gaia_m4_hr_bp_rp.png')
#plt.show()
plt.close()

g_avg = np.average(g_distance)
print('M4 avg distance: ', g_avg, 'number of stars: ',len(g_distance))

'''
now we want to take the HR we just created and fit isochrones to it
'''

#converting the apparent mags to abs values
tail_val =190
abs_G = distance_modulus(df['Gaia_G'].tail(tail_val),dist_vega)
abs_BP = distance_modulus(df['Gaia_BP'].tail(tail_val),dist_vega)
abs_RP = distance_modulus(df['Gaia_RP '].tail(tail_val),dist_vega)


#we then move the isochrones so that they better fir the HR
#bp_rp = (abs_BP - abs_RP) +0.55
#abs_G = abs_G +2

bp_rp = df['Gaia_BP'].tail(250) -df['Gaia_RP '].tail(250)
app_g = df['Gaia_G'].tail(250)

bp_rp = (bp_rp) +0.55
app_g = app_g +12.6

#plot

plt.figure(3)
plt.plot(bp_rp,app_g,'-',markersize=5,label='ZAMS')
plt.plot(gaia_results['bp_rp'],gaia_results['phot_g_mean_mag'],'.',color = 'black',alpha=0.6,markersize = 1.5,label='M4 HR-Diagram')
plt.gca().invert_yaxis()
plt.ylabel('G mag')
plt.xlabel('BP - RP mag')
plt.title('HR Diagram of M4 with Fitted Isochrones')
plt.legend()
plt.savefig('gaia_m4_isoschrones_fitted.png')
plt.close()
