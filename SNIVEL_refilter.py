#!/usr/bin/env python
import numpy
import math
import georinex as gr
import SNIVEL_orbits
import datetime
from SNIVEL_filedownloader import getbcorbit, getrinexhr, getrinex
from SNIVEL_tools import azi_elev, klobuchar, ecef2lla, gpsweekdow, getklobucharvalues, niell, dxyz2dneu, printProgressBar, niell_wet
from SNIVEL_tools import writesac, filter_int_VTEC
import os
from scipy.optimize import lsq_linear
#####################################################################################
#REDO THIS SECTION WITH IONOSPHERE INFORMATION
#SNIVEL.py
#This code will process the input files on the dates specified and produce velocities
#This uses the ionospheric free range combination to determine velocity
#Written by Brendan Crowell, University of Washington
#Last modified February 18, 2020
#Input files of testing sites and testing dates is required, the format is rather simple
#Follow the examples provided
#####################################################################################
homedir = '/home/jghent/SNIVEL/'
homedir = '/home/crowellb/SNIVEL_ION_CURRENT/'
event='tohoku'
c = 299792458.0 #speed of light
fL1 = 1575.42e6 #L1 frequency
fL2 = 1227.60e6 #L2 frequency
TECcons = 40.308e16*(fL1**2-fL2**2)/(fL1**2)/(fL2**2)
wL1 = c/fL1 #L1 wavelength
wL2 = c/fL2 #L2 wavelength
sitefile=homedir+'sites_'+event+'.txt'
datefile=homedir+'dates_'+event+'.txt'
#sampersec=1 #samples per second...this is now pulled from header
elevmask = 18 #elevation mask
#####################################################################################
##Preprocessing
#This section reads in the RINEX, apriori locations and orbit files and forms an
#observable file in the output folder
#I suggest only running up to 30 minutes of data at a time. You can process an entire day
#but it will take awhile. To modify this, the last two columns in dates_process.txt
#have the start time and the number of minutes you wish to process. Remember to include enough
#pre-event time and to account for leap seconds. RINEX files are in GPS time and it is 18 s
#ahead of UTC as of Feb, 2020. I perform no corrections for leap seconds.
with open(sitefile, 'rt') as g:
    rows = (line.split() for line in g)
    for grow in rows:
        site = grow[0] #looping over sites to process
        with open(datefile, 'rt') as f:
                rows2 = (line.split() for line in f)
                for grow2 in rows2:
                        year = grow2[0] #looping over years and days of year to process
                        doy = grow2[1]
                        st = grow2[2]
                        nummin = grow2[3]
                        #print('Processing station ', site, ' on year and  day ', year,  doy)
                        if not os.path.exists('output'): #if output folder doesn't exist, make it
                            os.makedirs('output')
                        if not os.path.exists('output/fig'): #if output folder doesn't exist, make it
                            os.makedirs('output/fig')
                        if not os.path.exists('output/fig/'+event): #if output folder doesn't exist, make it
                            os.makedirs('output/fig/'+event)
                        #obtain gps week and day of week for sp3 file download...not needed here
                        [gpsweek,gpsdow]=gpsweekdow(int(year),int(doy))
                        week = str(int(gpsweek)) #convert to strings
                        dow = str(int(gpsdow))

                        try:
                            obsfile = homedir+'rinex/' + site + doy + '0.' +  year[-2:] + 'o' #rinex file name
                            os.system('./teqc -R -S -E -C -J -phc -st ' + st + ' +dm ' + nummin + ' ' + obsfile + ' > example.o')
                            header=gr.rinexheader('example.o')#read the RINEX header
                            (x0,y0,z0)=header['position'] #use the a priori location of the site from the RINEX header. If its really bad, you might want to change it
                            samfreq = header['interval'] #interval between observations
                            sampersec = 1/float(15)
                            [latsta,lonsta,altsta]=ecef2lla(float(x0),float(y0),float(z0)) #station lat and lon are needed for klobuchar correction
                            print (site, latsta*180/math.pi, lonsta*180/math.pi, altsta)
                            veloutfile = homedir+'output/ion_' + site + '_' + doy + '_' + year + '.txt'
                            filter_int_VTEC(veloutfile,site,doy,year,event,sampersec,homedir)
                            print ('Station ', site, ' complete')
                        except Exception as error:
                            print ('Station ', site, ' not available on date ', str(year), ' ', str(doy))
                            print(error)
