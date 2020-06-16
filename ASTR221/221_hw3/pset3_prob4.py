'''

garcia,gil
astr231 - HW3 prob 4
2/28/2019


In this program, we find the mass - luminosity relationship
to find number density and mass density of stars in the local neighborhood
'''

#we import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
#we import the gaia catalog
from astroquery.gaia import Gaia
from astropy.coordinates import Angle

gaia_dr2 = Gaia.load_table('gaiadr2.gaia_source')

#we query the selected products from Gaia so
job = Gaia.launch_job_async("SELECT \
parallax, phot_rp_mean_mag, phot_g_mean_mag, \
bp_rp, bp_g,g_rp \
FROM gaiadr2.gaia_source \
WHERE parallax BETWEEN 100 AND 780 \
AND NOT (phot_g_mean_mag > 18 AND parallax > 100)\
AND phot_g_mean_mag < 1000\
AND bp_rp < 1000\
ORDER BY parallax DESC\
")
#printing our job results so that we can analyze
gaia_results = job.get_results()
gaia_results = gaia_results[0:100]
#finding the distance to our stars using the parallax given
g_distance = 1 / (gaia_results['parallax']/1000.)
#we find the bolometric correction using the abs mag and temperature of the stars
#finding the temp
temp = 10**((gaia_results['bp_rp'] - 14.551) / -3.684)
#calculating bol_correction
bol_correct = -8.499*(np.log10(temp) - 4)**4 + 13.421*(np.log10(temp)-4)**3 - 8.131*(np.log10(temp)- 4)**2 - 3.901*(np.log10(temp)- 4) - 0.438
#generating our bol mags
abs_mag = gaia_results['phot_g_mean_mag'] + bol_correct
#calculating luminosity
lum = 10**((abs_mag  - 4.75) / 2.5)
#calculating the mass from M-L reln
mass = (lum)**0.2875
#calculating mass in grams
mass_in_grams = mass * 1.989e33
#finding total amount of mass from our query
tot_grams = np.sum(mass_in_grams)
#printing our the total amount of mass
print("Total amount of mass in grams: ",tot_grams)
#finding the furthest star from our query
far_star = g_distance[99]
print("The furthest star in our catalog in pc: ",far_star)
#finding the volume enclosed in our query
vol = (4/3) * np.pi * (far_star ** 3)
strs_per_c_pc = 100 /vol
print("stars per cubic pc: ",strs_per_c_pc)
#counting the number of white dwarfs and main sequence stars in our query
white_dwarf = 0
main_seq = 0
for item in gaia_results:
    if item[2] > 17:
        white_dwarf += 1
    else:
         main_seq += 1
print("There is %d white dwarfs in our query and %d main sequence stars." %(white_dwarf, main_seq))
