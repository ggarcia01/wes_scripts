'''
garcia, gil
ASTR231 - hw3 prob 2
2/28/2019


In this program, we create an HR diagram of the
300 brightest stars and of the 300 closest stars
using Gaia data for both
'''

#importing libraries to run math operations and to plot
import numpy as np
import matplotlib.pyplot as plt
#we import Gaia so we can query Gaia DR2
from astroquery.gaia import Gaia
#we import coordinates to handle our celestial coordinates
from astropy.coordinates import Angle
# prints out all the columns available in gaia dr2 so we know what we have to work with
gaia_dr2 = Gaia.load_table('gaiadr2.gaia_source')
#for column in gaia_dr2.get_columns():
#    print(column.get_name())
#function1 - given apparent mag and distance, using the distance modulus to find abs mag
#we will use this to find abs mag of the stars that we query
def distane_modulus(apparant_m,distance):
    term2 = 5*(np.log10(distance)-1)
    return (apparant_m - term2)
#we query gaia dr2 using sql
#we first query to get the 100 brightest stars in terms of apparent mag
job = Gaia.launch_job_async("SELECT \
parallax, phot_rp_mean_mag, phot_g_mean_mag, \
bp_rp, bp_g,g_rp \
FROM gaiadr2.gaia_source \
WHERE phot_g_mean_mag < 3.5 \
ORDER BY phot_g_mean_mag\
")
#print(job)
gaia_results = job.get_results()
#only the 100 brightest stars:
gaia_results = gaia_results[0:300]
#print(gaia_results)
#finding abs G mag using parallax and apparent G mag
g_distance = 1/ (gaia_results['parallax']/1000.)

abs_g_mags = distane_modulus(gaia_results['phot_g_mean_mag'],g_distance)
#creating an HR diagram of the 100 brightest stars as seen from Earth
plt.plot(gaia_results['bp_rp'],abs_g_mags,'.',color = 'black',markersize = 4,alpha=0.8,label='brightest')
#plt.gca().invert_yaxis()
#plt.xlabel('BP - RP mag')
#plt.ylabel('G mag')
#plt.show()

#now we query the 100 closest stars
job1 = Gaia.launch_job_async("SELECT \
parallax, phot_rp_mean_mag, phot_g_mean_mag, \
bp_rp, bp_g,g_rp \
FROM gaiadr2.gaia_source \
WHERE parallax >100 \
AND parallax < 770 \
AND NOT (phot_g_mean_mag > 18 AND parallax > 100) \
ORDER BY parallax")
gaia_results1 = job1.get_results()
gaia_results1.sort('parallax')
#gaia_results1.reverse()
#only the 100 brightest stars:
gaia_results1 = gaia_results1[0:300]
#print(gaia_results1)
#print(gaia_results1['phot_g_mean_mag'])
#finding abs G mag using parallax and apparent G mag
g_distance1 = 1/ (gaia_results1['parallax']/1000.)
abs_g_mags1 = distane_modulus(gaia_results1['phot_g_mean_mag'],g_distance1)
#creating an HR diagram of the 100 brightest stars as seen from Earth
plt.plot(gaia_results1['bp_rp'],abs_g_mags1,'.',color = 'red',markersize = 4,alpha=0.8,label='closest')
plt.gca().invert_yaxis()
plt.xlabel('BP - RP mag')
plt.ylabel('G mag')
plt.legend()
plt.title("HR Diagram of Select Stars")
plt.show()
