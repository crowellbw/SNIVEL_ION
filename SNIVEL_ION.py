#!/usr/bin/env python
import numpy
import math
import georinex as gr
import SNIVEL_orbits
import datetime
from SNIVEL_filedownloader import getbcorbit, getrinexhr, getrinex, getrinexNGS
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
homedir='/home/crowellb/SNIVEL_ION_CURRENT/'
event='greenland'
c = 299792458.0 #speed of light
fL1 = 1575.42e6 #L1 frequency
fL2 = 1227.60e6 #L2 frequency
TECcons = 40.308e16*(fL1**2-fL2**2)/(fL1**2)/(fL2**2)
wL1 = c/fL1 #L1 wavelength
wL2 = c/fL2 #L2 wavelength
sitefile='sites_'+event+'.txt'
datefile='dates_'+event+'.txt'
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
                        outfile = 'output/observables_' + site + '_' + doy + '_' + year + '.txt'
                        ffo = open(outfile,'w')
                        #obtain gps week and day of week for sp3 file download...not needed here
                        [gpsweek,gpsdow]=gpsweekdow(int(year),int(doy))
                        week = str(int(gpsweek)) #convert to strings
                        dow = str(int(gpsdow))

                        getbcorbit(year, doy) #Download broadcast orbit
                        try:
                            getrinexNGS(site, year, doy) #Download RINEX

                            obsfile = 'rinex/' + site + doy + '0.' +  year[-2:] + 'o' #rinex file name
                            navfile = 'nav/brdc' + doy + '0.' +  year[-2:] + 'n' #broadcast navigation file name

                            #use teqc to convert the rinex to remove comments within and remove non gps satellites
                            #Also, the -st command sets the start time and +dm is the number of minutes to consider
                            print('here')
                            os.system('./teqc -R -S -E -C -J -phc -st ' + st + ' +dm ' + nummin + ' ' + obsfile + ' > example.o')
                            print('here')
                            #os.system('tac example.o | sed -e "/post/{N;d;}" | tac > example2.o') #remove the comment lines that teqc didn't remove
                            header=gr.rinexheader('example.o')#read the RINEX header
                            (x0,y0,z0)=header['position'] #use the a priori location of the site from the RINEX header. If its really bad, you might want to change it
                            samfreq = 30 #header['interval'] #interval between observations
                            sampersec = 1/float(samfreq)
                            [latsta,lonsta,altsta]=ecef2lla(float(x0),float(y0),float(z0)) #station lat and lon are needed for klobuchar correction
                            print (site, latsta*180/math.pi, lonsta*180/math.pi, altsta)
                            nav = gr.load(navfile) #Load the broadcast navigation file
                            [alpha,beta]=getklobucharvalues(navfile) #get klobuchar corrections from navigation file

                            obs = gr.load('example.o') #Load the RINEX file
                            
                            L1 = obs['L1'].values#phase values
                            L2 = obs['L2'].values
                            S1 = obs['S1'].values
                            S2 = obs['S2'].values

                            obs_time = obs.time.values #observation times
                            #samfreq = float(obs_time[1])-float(obs_time[0])
                            #sampersec = 1/float(samfreq)
                            print(sampersec)
                            nt = len(obs_time) #number of observation times
                            svs = obs.sv #satellites observed
                            ns = len(svs) #number of satellites observed

                            for i in range (0, nt):
                                printProgressBar(i,nt, prefix = 'Writing Obs File:',suffix = 'Complete', length = 25)
                                gps_time = (numpy.datetime64(obs_time[i]) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's') #compute gps time, seconds from Jan 6, 1980
                                gps_sow = float(gps_time) - gpsweek*604800 #determine gps second of week
                                for j in range (0, ns):
                                    total_val = float(float(L2[i,j])+float(L1[i,j])) #if any of these 2 observables don't exist, don't record
                                    if (math.isfinite(total_val)):
                                        sv1 = svs.sv[j]
                                        if "G" in str(sv1.values):
                                            [x,y,z,tr,rho,tsv]=SNIVEL_orbits.bcorbit(nav,sv1.values,gps_sow,x0,y0,z0) #compute satellite location, clock, and relativistic term with broadcast orbits
                                            [azi,elev]=azi_elev(x0,y0,z0,x,y,z) #compute azimuth and elevation angles. Broadcast orbits are good enough here
                                            Mdry=niell(elev,latsta*180/math.pi,altsta,int(doy))#this is the niell dry troposphere delay
                                            Mwet=niell_wet(elev,latsta*180/math.pi)#Niell wet tropospheric delay
                                            [dIon1,dIon2,latIPP,lonIPP]=klobuchar(float(latsta)*180/math.pi,float(lonsta)*180/math.pi,float(elev),float(azi),gps_sow,alpha,beta) #compute the klobuchar correction for broadcast
                                            rclock=(tr+tsv)*c #relativistic clock correction with broadcast clock correction
                                            #Reformat to strings to output to observation file
                                            l1 = "{0:.5f}".format(float(L1[i,j])*wL1+rclock) #L1 with clock corrections in meters
                                            l2 = "{0:.5f}".format(float(L2[i,j])*wL2+rclock) #L2 with clock corrections in meters
                                            rx0 = "{0:.9f}".format(float((x0-x)/rho)) #Direction cosine, x
                                            ry0 = "{0:.9f}".format(float((y0-y)/rho)) #Direction cosine, y
                                            rz0 = "{0:.9f}".format(float((z0-z)/rho)) #Direction cosine, z
                                            az = "{0:.2f}".format(float(azi)) #azimuth between receiver and satellite
                                            el = "{0:.2f}".format(float(elev)) #elevation angle between receiver and satellite
                                            rhotrue = "{0:.5f}".format(float(rho)) #true range between satellite and a priori position of station
                                            di1 = "{0:.5f}".format(dIon1) #Klobuchar ionospheric delay on L1
                                            di2 = "{0:.5f}".format(dIon2) #Klobuchar ionospheric delat on L2
                                            shd = "{0:.5f}".format(Mdry) #Niell slant hydrostatic delay
                                            swd = "{0:.5f}".format(Mwet) #Niell slant wet delay
                                            gpst = "{0:.2f}".format(float(gps_time)) #gps time, continuous seconds since Jan 6, 1980
                                            gpsw = "{0:.0f}".format(float(gpsweek)) #gps week
                                            gpss = "{0:.2f}".format(float(gps_sow)) #gps second of week
                                            svstrip = str(sv1.values)[1:]
                                            dxsat = "{0:.6f}".format(float((x0-x))) #distance between receiver and satellite in x
                                            dysat = "{0:.6f}".format(float((y0-y))) #distance between receiver and satellite in y
                                            dzsat = "{0:.6f}".format(float((z0-z))) #distance between receiver and satellite in z
                                            s1 = "{0:.2f}".format(float(S1[i,j]))
                                            s2 = "{0:.2f}".format(float(S2[i,j]))
                                            #s1 = "{0:.2f}".format(0)
                                            #s2 = "{0:.2f}".format(0)


                                            latipp = "{0:.4f}".format(float(latIPP*180))
                                            lonipp = "{0:.4f}".format(float(lonIPP*180))
                                            ffo.write(str(i)+' '+gpst+' '+gpsw+' '+gpss+' '+svstrip+' '+rx0+' '+ry0+' '+rz0+' '+l1+' '+l2+' '+az+' '+el+' '+di1+' '+di2+' '+shd+' '+dxsat+' '+dysat+' '+dzsat+' '+swd+' '+s1+' '+s2+' '+latipp+' '+lonipp+'\n')

                            ffo.close()
                            #####################################################################################
                            #Load the observables file, compute velocities
                            #Refer to section above about variable names
                            #####################################################################################
                            veloutfile = 'output/ion_' + site + '_' + doy + '_' + year + '.txt'
                            ffo = open(veloutfile,'w')
                            a = numpy.loadtxt(outfile)
                            tind = a[:,0]
                            gtime = a[:,1]
                            svind = a[:,4]
                            i1corr = a[:,12]
                            i2corr = a[:,13]
                            shdcorr = a[:,14]

                            rx = a[:,5]
                            ry = a[:,6]
                            rz = a[:,7]
                            l1 = a[:,8]
                            l2 = a[:,9]
                            azi = a[:,10]
                            el = a[:,11]
                            dxsat = a[:,15]
                            dysat = a[:,16]
                            dzsat = a[:,17]
                            latipp = a[:,21]
                            lonipp = a[:,22]
                            tstart = numpy.amin(tind)+1
                            tstop = numpy.amax(tind)
                            for i in range (int(tstart),int(tstop)+1):
                                a0 = numpy.where(tind == i-1)[0]
                                a1 = numpy.where(tind == i)[0]
                                l10 = l1[a0]
                                l11 = l1[a1]
                                l20 = l2[a0]
                                l21 = l2[a1]
                                sv0 = svind[a0]
                                sv1 = svind[a1]
                                dxsat0 = dxsat[a0]
                                dxsat1 = dxsat[a1]
                                dysat0 = dysat[a0]
                                dysat1 = dysat[a1]
                                dzsat0 = dzsat[a0]
                                dzsat1 = dzsat[a1]
                                rx0 = rx[a0]
                                rx1 = rx[a1]
                                ry0 = ry[a0]
                                ry1 = ry[a1]
                                rz0 = rz[a0]
                                rz1 = rz[a1]

                                el1 = el[a1]
                                az1 = azi[a1]
                                latipp1 = latipp[a1]
                                lonipp1 = lonipp[a1]
                                Vdat = list() #variometric data
                                if ((gtime[a1[0]]-gtime[a0[0]])  < 1/sampersec+0.005):
                                    for j in range (0, len(sv1)):
                                        asv = numpy.where(sv0 == sv1[j])[0]
                                        if (len(asv)>0 and el1[j] > elevmask):
                                            dran0 = math.sqrt(math.pow(dxsat0[int(asv)],2)+math.pow(dysat0[int(asv)],2)+math.pow(dzsat0[int(asv)],2))
                                            dran1 = math.sqrt(math.pow(dxsat1[j],2)+math.pow(dysat1[j],2)+math.pow(dzsat1[j],2))
                                            dran = dran1-dran0
                                            l1diff = l11[j]-l10[int(asv)]-dran#time difference of l1
                                            l2diff = l21[j]-l20[int(asv)]-dran#time difference of l2
                                            LG = l11[j]-l21[j] #geometry free combination of current epoch...for absolute TEC
                                            varvalL = [l1diff-l2diff] #variometric value of geometry free combination
                                            gpst = "{0:.2f}".format(float(gtime[a1[0]]))#reformat for output
                                            satnum = "{0:.0f}".format(float(sv1[j]))
                                            vTEC = "{0:.5e}".format(1/TECcons*float(varvalL[0])*sampersec)
                                            aTEC = "{0:.5e}".format(1/TECcons*float(LG))
                                            elev = "{0:.5f}".format(float(el1[j]))
                                            azim = "{0:.5f}".format(float(az1[j]))
                                            latpp = "{0:.4f}".format(float(latipp1[j]))#ionospheric piercing points
                                            lonpp = "{0:.4f}".format(float(lonipp1[j]))
                                            ffo.write(str(i)+' '+satnum+' '+gpst+' '+aTEC+' '+vTEC+' '+elev+' '+azim+' '+latpp+' '+lonpp+'\n')
                            ffo.close()
                            filter_int_VTEC(veloutfile,site,doy,year,event,sampersec,homedir)
                            print ('Station ', site, ' complete')
                        except Exception:
                            print ('Station ', site, ' not available on date ', str(year), ' ', str(doy))
