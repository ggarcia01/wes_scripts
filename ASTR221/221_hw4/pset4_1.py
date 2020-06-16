# Gil Garcia
# ASTR221 - hw3 prob 1
#4/4/2019

'''
In this script, we query the Gaia DR2 catalog to create an HR Diagram for the open cluster NGC188
with fitted isochrones
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
import pandas as pd

def distance_modulus(apparant_m,distance):
    term2 = 5*(np.log10(distance)-1)
    return (apparant_m - term2)


'''importing isochrones data'''
dist_vega = 10**((0.2*5)+(0.2*1.3847))

df = pd.read_csv('0point0_age6point5.csv')

# sql querying the gaia catalog
job = Gaia.launch_job_async("SELECT \
pmra, pmdec, parallax, parallax_error, parallax_over_error, phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_mag, \
bp_rp, bp_g,g_rp \
FROM gaiadr2.gaia_source \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec), CIRCLE('ICRS',12.11,85.255,2))=1 \
AND parallax_error < 2 \
AND parallax_over_error > 5 \
AND pmra BETWEEN -3 and -1.5 \
AND pmdec BETWEEN -2 and 0 \
AND parallax BETWEEN 0.50 and 1.5")
gaia_results = job.get_results()

# creating proper motion plt
plt.figure(1)
plt.plot(gaia_results['pmra'],gaia_results['pmdec'],'.',color = 'maroon',alpha = 0.3)
plt.xlabel('Proper Motion in RA')
plt.ylabel('Proper Motion in DEC')
plt.title('Proper Motion of NGC188 Using Gaia Catalog')
plt.grid()
plt.savefig('gaia_ngc188_proper_motion.png')
plt.close()

#finding abs G mag using parallax and apparent G mag
g_distance = 1/ (gaia_results['parallax']/1000.)
abs_g_mags = distance_modulus(gaia_results['phot_g_mean_mag'],g_distance)

# creating an HR diagram
plt.figure(2)
plt.plot(gaia_results['bp_rp'],abs_g_mags,'.',color = 'black',alpha=0.5,markersize = 2.3)
plt.gca().invert_yaxis()
plt.xlabel('BP - RP mag')
plt.ylabel('G mag')
plt.title('HR Diagram Using NGC188 in Gaia Catalog')
plt.savefig('gaia_ngc188_hr_bp_g.png')
plt.close()

# average distance as calculated from gaia:
g_avg = np.average(g_distance)
#printing avg distance and len of list (number of stars)
print('NGC188 avg distance: ', g_avg, 'number of stars: ',len(g_distance))




'''plotting 3rd plot'''


#tail_val =361
#abs_G = distance_modulus(df['Gaia_G'],dist_vega)
#abs_BP = distance_modulus(df['Gaia_BP'],dist_vega)
#abs_RP = distance_modulus(df['Gaia_RP'],dist_vega)

bp_rp = df['Gaia_BP']-df['Gaia_RP']
app_g = df['Gaia_G']

bp_rp = (bp_rp) +0.8
app_g = app_g +15

# creating an HR diagram using Pleiades data from Gaia
plt.plot(gaia_results['bp_rp'],gaia_results['phot_g_mean_mag'],'.',color = 'black',alpha=0.6,markersize = 1.5,label='NGC HR-Diagram')
plt.plot(bp_rp,app_g,'-',markersize=5,label='ZAMS')
plt.gca().invert_yaxis()
plt.xlabel('BP - RP mag')
plt.ylabel('G mag')
plt.title('HR Diagram of NGC188 with Fitted Isochrones')
plt.legend()
plt.savefig('gaia_ngc188_isoschrones_fitted.png')
plt.close()
