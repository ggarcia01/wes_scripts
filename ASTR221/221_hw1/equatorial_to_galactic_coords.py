'''
Gil Garcia
ASTR221 - Galactic Astronomy
01/31/2019

In this script, we convert equatorial coordinates into galactic coordinates.
'''
#import the necessary libraries
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
#ask user to enter the equatorial coordinates
#asking for RA
print('Enter RA in decimal form:')
ra = float(input())
#asking for DEC
print('Enter DEC in decimal form:')
dec = float(input())
#we define some constants, all in degrees
ag,dg = 192.85948, 27.12825 #equatorial coordinates of north galactic pole
l_ncp = 122.93192 #galactic longtitude of the north celestial pole
#we now convert the given equatorial coordinates (ra,dec) into galactic
#   coordinates (l,b) using the following equations:
def eq_to_gal(ra,dec): # must give ra,dec in decimal degrees for fcn to work, returns coords in degrees
    rad_ag, rad_dg = np.deg2rad(192.85948),np.deg2rad(27.12825) #equatorial coordinates of north galactic pole
    rad_l_ncp = np.deg2rad(122.93192) #galactic longtitude of the north celestial pole
    rad_ra,rad_dec = np.deg2rad(ra),np.deg2rad(dec) # converting the given coords to radians
    #we use the following formula to calculate b - galactic latitude
    b = np.arcsin(np.sin(rad_dec)*np.sin(rad_dg) + np.cos(rad_dec)*np.cos(rad_dg)*np.cos(rad_ra-rad_ag))
    #the function to find l - galactic longtitude - is the arctan of (l_ncp minus the fraction(x1/x2))
    x1 = np.cos(rad_dec)*np.sin(rad_ra - rad_ag)
    x2 = np.sin(rad_dec)*np.cos(rad_dg) - np.cos(rad_dec)*np.sin(rad_dg)*np.cos(rad_ra - rad_ag)
    l = rad_l_ncp - np.arctan2(x1,x2) #use arctan2 so that we get the correct sign by calculating the angle between x1 and x2
    l,b = np.rad2deg(l),np.rad2deg(b)
    if l < 0: #we want l to be between 0 and 360 degrees
        l = l + 360
    return l,b
print()
print('Equatorial Coordinates (ra,dec):',ra,',',dec)
print()
print('Our Calculation:')
l,b = eq_to_gal(ra,dec)
print("Galactic Coordinates (l,b):",l,',',b)
print()

'''
#checking our answer with the answer given by the astropy library
# setting up the coordinates for conversion
#coord = SkyCoord(ra,dec,unit=(u.hourangle, u.deg))
coord = SkyCoord(ra,dec,frame='icrs',unit='deg')
#converting to galactic coordinates
g_coords = coord.galactic
print('Astropy Calculation:')
print("Galactic Coordinates (l,b):",g_coords.l.degree,',',g_coords.b.degree)

#we now want to make sure that our function works for all possible equatorial coordinates
ras = np.arange(0,360,1) # create an array of all possible ra's
decs = np.arange(-90,90,0.5) # create an array of all possible dec's
counter = 0
#create a for loop that checks all combinations of coordinates
for ra in ras:
    for dec in decs:
        lo,bo = eq_to_gal(ra,dec) #convert using our function
        coord = SkyCoord(ra,dec,frame='icrs',unit='deg') #convert unsing astropy
        g_coords = coord.galactic
        if abs(lo - g_coords.l.degree) > 1: #if the answer is wrong, tell me what coordinate doesn't work
            print(ra,dec)
    print('Done',counter)
    counter += 1

#we verify that it works for all coordinates so our function works!
'''
