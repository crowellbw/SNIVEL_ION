#!/usr/bin/env python
import numpy
import datetime
import calendar
import math
import georinex as gr
import obspy
import scipy
from scipy.signal import butter, lfilter, filtfilt
import matplotlib.pyplot as plt
from obspy.io.sac import SACTrace
import matplotlib.pyplot as plt
import geopy.distance
#####################################################################################
#SNIVEL_tools.py
#Written by Brendan Crowell, University of Washington
#Last edited January 11, 2019
#####################################################################################
c = 299792458.0 #speed of light
fL1 = 1575.42e6 #L1 frequency
fL2 = 1227.60e6 #L2 frequency
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = 'â–ˆ', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def month_converter(month):
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    return months.index(month) + 1

def doy2month(doy,year):
    isleap = calendar.isleap(year)
    if str(isleap) == 'True':
        dom = [31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        dom = [31,28,31,30,31,30,31,31,30,31,30,31]

    dayind = int(dom[0])
    monind=1
    for i in range (1, 12):
        if (int(doy) < dayind):
            month = monind
        else:
            dayind = dayind+dom[i]
            monind=monind+1
            
    return month

def doy_calc(year,month,day):
    isleap = calendar.isleap(year)
    if str(isleap) == 'True':
        dom = [31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        dom = [31,28,31,30,31,30,31,31,30,31,30,31]
    doy = int(numpy.sum(dom[:(month-1)])+day)
    return(doy)

def gpsweekdow(year,doy):
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
    gpstime = (numpy.datetime64(date) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's')
    gpsweek = int(gpstime/604800)
    gpsdow = math.floor((gpstime-gpsweek*604800)/86400)                   
    return(gpsweek, gpsdow)

def gpsleapsec(gpssec):
    leaptimes = numpy.array([46828800, 78364801, 109900802, 173059203, 252028804, 315187205, 346723206, 393984007, 425520008, 457056009, 504489610, 551750411, 599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017])
    leapseconds = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])
    a1 = numpy.where(gpssec > leaptimes)[0]
    leapsec = len(a1)
    return(leapsec)

def lla2ecef(lat,lon,alt):
    lat = lat*math.pi/180
    lon = lon*math.pi/180
    a = 6378137
    e = 8.1819190842622e-2

    N = a/numpy.sqrt(1-numpy.power(e,2)*numpy.power(numpy.sin(lat),2))

    x = (N+alt)*numpy.cos(lat)*numpy.cos(lon)
    y = (N+alt)*numpy.cos(lat)*numpy.sin(lon)
    z = ((1-numpy.power(e,2))*N+alt)*numpy.sin(lat)

    return (x, y, z)


def ecef2lla(x,y,z):
    a = 6378137
    e = 8.1819190842622e-2

    b = math.sqrt(math.pow(a,2)*(1-math.pow(e,2)))
    ep = math.sqrt((math.pow(a,2)-math.pow(b,2))/math.pow(b,2))
    p = math.sqrt(math.pow(x,2)+math.pow(y,2))
    th = math.atan2(a*z,b*p)
    lon = math.atan2(y,x)
    lat = math.atan2((z+math.pow(ep,2)*b*math.pow(math.sin(th),3)),(p-math.pow(e,2)*a*math.pow(math.cos(th),3)))
    N = a/math.sqrt(1-math.pow(e,2)*math.pow(math.sin(lat),2))
    alt = p/math.cos(lat)-N

    return (lat,lon,alt)


def azi_elev(xsta,ysta,zsta,xsat,ysat,zsat):
    [latsta,lonsta,altsta]=ecef2lla(xsta,ysta,zsta)
    [latsat,lonsat,altsat]=ecef2lla(xsat,ysat,zsat)
    re = math.sqrt(math.pow(xsta,2)+math.pow(ysta,2)+math.pow(zsta,2))
    rs = math.sqrt(math.pow(xsat,2)+math.pow(ysat,2)+math.pow(zsat,2))
    gamma = math.acos(math.cos(latsta)*math.cos(latsat)*math.cos(lonsat-lonsta) + math.sin(latsta)*math.sin(latsat))
    elev = math.acos(math.sin(gamma)/math.sqrt(1 + math.pow(re/rs,2) - 2*re/rs*math.cos(gamma)))

    deltalon = lonsat-lonsta

    azi = math.atan2(math.sin(deltalon)*math.cos(latsat),math.cos(latsta)*math.sin(latsat)-math.sin(latsta)*math.cos(latsat)*math.cos(deltalon))

    azi = azi*180/math.pi

    if (azi < 0):
        azi = azi+360
    elev = elev*180/math.pi
    return(azi,elev)


def getklobucharvalues(navfile):
    navhead = gr.rinexheader(navfile) #load navigation header to obtain klobuchar
    alpha = navhead['ION ALPHA'] #klobuchar alphas
    beta = navhead['ION BETA'] #klobuchar betas
    alpha2 = alpha.replace("D", "E") #nav headers sometimes use D instead of E for power
    beta2 = beta.replace("D", "E")
    alp1= numpy.asarray(alpha2.split())
    bet1= numpy.asarray(beta2.split())
    alp = alp1.astype(numpy.float) #numpy array of klobuchar alphas
    bet = bet1.astype(numpy.float) #numpy array of klobuchar betas
    return(alp,bet)

#Klobuchar ionospheric correction

def klobuchar(latsta,lonsta,elev,azimuth,tow,alpha,beta):
    a1 = float(alpha[0])
    a2 = float(alpha[1])
    a3 = float(alpha[2])
    a4 = float(alpha[3])
    b1 = float(beta[0])
    b2 = float(beta[1])
    b3 = float(beta[2])
    b4 = float(beta[3])
    

    elev = elev/180 #elevation angle in semicircles
    azimuth = azimuth*math.pi/180

    psi = 0.0137/(elev+0.11) - 0.022 #earth centered angle

    lat_i = latsta/180 + psi*math.cos(azimuth) #subionospheric latitude

    if (lat_i > 0.416):
        lat_i = 0.416
    if (lat_i < -0.416):
        lat_i = -0.416

    lon_i = lonsta/180 + psi*math.sin(azimuth)/math.cos(lat_i*math.pi) #subionospheric longitude

    lat_m = lat_i + 0.064*math.cos((lon_i - 1.617)*math.pi) #geomagnetic lat

    t = 4.32e4*lon_i + tow
    t = t % 86400
    if (t > 86400):
        t = t - 86400
    if (t < 0):
        t = t + 86400

    sF = 1 + 16*math.pow(0.53-elev,3) #slant factor

    PER = b1 + b2*lat_m + b3*math.pow(lat_m,2) + b4*math.pow(lat_m,3)

    if (PER < 72000):
        PER = 72000

    x = 2*math.pi*(t-50400)/PER

    AMP = a1 + a2*lat_m + a3*math.pow(lat_m,2) + a4*math.pow(lat_m,3)

    if (AMP < 0):
        AMP = 0

    if (abs(x) > 1.57):
        dIon1 = sF*(5e-9)
    else:
        dIon1 = sF*(5e-9 + AMP*(1 - math.pow(x,2)/2 + math.pow(x,4)/24))

    dIon1 = c*dIon1

    dIon2 = math.pow(fL1/fL2,2)*dIon1

    return(dIon1,dIon2,lat_i,lon_i)


def niell(elev,lat,alt,doy):
    aht = 2.53e-5
    bht = 5.49e-3
    cht = 1.14e-3
    
    aavg15 = 1.2769934e-3
    bavg15 = 2.9153695e-3
    cavg15 = 62.610505e-3
    aamp15 = 0.0
    bamp15 = 0.0
    camp15 = 0.0
    
    aavg30 = 1.2683230e-3
    bavg30 = 2.9152299e-3
    cavg30 = 62.837393e-3
    aamp30 = 1.2709626e-5
    bamp30 = 2.1414979e-5
    camp30 = 9.0128400e-5

    aavg45 = 1.2465397e-3
    bavg45 = 2.9288445e-3
    cavg45 = 63.721774e-3
    aamp45 = 2.6523662e-5
    bamp45 = 3.0160779e-5
    camp45 = 4.3497037e-5

    aavg60 = 1.2196049e-3
    bavg60 = 2.9022565e-3
    cavg60 = 63.824265e-3
    aamp60 = 3.4000452e-5
    bamp60 = 7.2562722e-5
    camp60 = 84.795348e-5

    aavg75 = 1.2045996e-3
    bavg75 = 2.9024912e-3
    cavg75 = 64.258455e-3
    aamp75 = 4.1202191e-5
    bamp75 = 11.723375e-5
    camp75 = 170.37206e-5
    

    if (abs(lat) <= 15):
        aavg = aavg15
        bavg = bavg15
        cavg = cavg15
        aamp = aamp15
        bamp = bamp15
        camp = camp15


    if (abs(lat) > 15 and abs(lat) <= 30):       
        amavg = 15.0/(aavg30-aavg15)
        aavg = (abs(lat) - 15)/amavg + aavg15
        amamp = 15.0/(aamp30-aamp15)
        aamp = (abs(lat) - 15)/amamp + aamp15

        bmavg = 15.0/(bavg30-bavg15)
        bavg = (abs(lat) - 15)/bmavg + bavg15
        bmamp = 15.0/(bamp30-bamp15)
        bamp = (abs(lat) - 15)/bmamp + bamp15

        cmavg = 15.0/(cavg30-cavg15)
        cavg = (abs(lat) - 15)/cmavg + cavg15
        cmamp = 15.0/(camp30-camp15)
        camp = (abs(lat) - 15)/cmamp + camp15

    if (abs(lat) > 30 and abs(lat) <= 45):       
        amavg = 15.0/(aavg45-aavg30)
        aavg = (abs(lat) - 30)/amavg + aavg30
        amamp = 15.0/(aamp45-aamp30)
        aamp = (abs(lat) - 30)/amamp + aamp30

        bmavg = 15.0/(bavg45-bavg30)
        bavg = (abs(lat) - 30)/bmavg + bavg30
        bmamp = 15.0/(bamp45-bamp30)
        bamp = (abs(lat) - 30)/bmamp + bamp30

        cmavg = 15.0/(cavg45-cavg30)
        cavg = (abs(lat) - 30)/cmavg + cavg30
        cmamp = 15.0/(camp45-camp30)
        camp = (abs(lat) - 30)/cmamp + camp30

    if (abs(lat) > 45 and abs(lat) <= 60):       
        amavg = 15.0/(aavg60-aavg45)
        aavg = (abs(lat) - 45)/amavg + aavg45
        amamp = 15.0/(aamp60-aamp45)
        aamp = (abs(lat) - 45)/amamp + aamp45

        bmavg = 15.0/(bavg60-bavg45)
        bavg = (abs(lat) - 45)/bmavg + bavg45
        bmamp = 15.0/(bamp60-bamp45)
        bamp = (abs(lat) - 45)/bmamp + bamp45

        cmavg = 15.0/(cavg60-cavg45)
        cavg = (abs(lat) - 45)/cmavg + cavg45
        cmamp = 15.0/(camp60-camp45)
        camp = (abs(lat) - 45)/cmamp + camp45

    if (abs(lat) > 60 and abs(lat) <= 75):       
        amavg = 15.0/(aavg75-aavg60)
        aavg = (abs(lat) - 60)/amavg + aavg60
        amamp = 15.0/(aamp75-aamp60)
        aamp = (abs(lat) - 60)/amamp + aamp60

        bmavg = 15.0/(bavg75-bavg60)
        bavg = (abs(lat) - 60)/bmavg + bavg60
        bmamp = 15.0/(bamp75-bamp60)
        bamp = (abs(lat) - 60)/bmamp + bamp60

        cmavg = 15.0/(cavg75-cavg60)
        cavg = (abs(lat) - 60)/cmavg + cavg60
        cmamp = 15.0/(camp75-camp60)
        camp = (abs(lat) - 60)/cmamp + camp60

    if (abs(lat) > 75):
        aavg = aavg75
        bavg = bavg75
        cavg = cavg75
        aamp = aamp75
        bamp = bamp75
        camp = camp75

    a = aavg - aamp*math.cos(2*math.pi*(doy-28)/365.25)
    b = bavg - bamp*math.cos(2*math.pi*(doy-28)/365.25)
    c = cavg - camp*math.cos(2*math.pi*(doy-28)/365.25)

    el = math.sin(elev*math.pi/180)

    m = (1 + a/(1+b/(1+c)))/(el + a/(el+b/(el+c)))
    mh = (1 + aht/(1+bht/(1+cht)))/(el + aht/(el+bht/(el+cht)))
    dm = (1/el - mh)*alt/1000

    Mdry = m + dm

    Tropdelay = 2.3*math.exp(-0.116e-3*alt)*Mdry


    return(Tropdelay)

#This takes displacements in x, y, z and converts them to north, east up
def dxyz2dneu(dx,dy,dz,lat,lon):
	lat = lat*math.pi/180
	lon = lon*math.pi/180
	dn = -numpy.sin(lat)*numpy.cos(lon)*dx-numpy.sin(lat)*numpy.sin(lon)*dy+numpy.cos(lat)*dz
	de = -numpy.sin(lon)*dx+numpy.cos(lon)*dy
	du = numpy.cos(lat)*numpy.cos(lon)*dx+numpy.cos(lat)*numpy.sin(lon)*dy+numpy.sin(lat)*dz
	return (dn, de, du)


def niell_wet(elev, lat):

    aavg15 = 5.8021897e-4
    bavg15 = 1.4275268e-3
    cavg15 = 4.3472961e-2

    aavg30 = 5.6794847e-4
    bavg30 = 1.5138625e-3
    cavg30 = 4.6729510e-2

    aavg45 = 5.8118019e-4
    bavg45 = 1.4572752e-3
    cavg45 = 4.3908931e-2

    aavg60 = 5.9727542e-4
    bavg60 = 1.5007428e-3
    cavg60 = 4.4626982e-2

    aavg75 = 6.1641693e-4
    bavg75 = 1.7599082e-3
    cavg75 = 5.4736038e-2

    if (abs(lat) <= 15):
        aavg = aavg15
        bavg = bavg15
        cavg = cavg15

    if (abs(lat) > 15 and abs(lat) <= 30):
        amavg = 15.0/(aavg30-aavg15)
        aavg = (abs(lat) - 15)/amavg + aavg15

        bmavg = 15.0/(bavg30-bavg15)
        bavg = (abs(lat) - 15)/bmavg + bavg15

        cmavg = 15.0/(cavg30-cavg15)
        cavg = (abs(lat) - 15)/cmavg + cavg15

    if (abs(lat) > 30 and abs(lat) <= 45):
        amavg = 15.0/(aavg45-aavg30)
        aavg = (abs(lat) - 30)/amavg + aavg30

        bmavg = 15.0/(bavg45-bavg30)
        bavg = (abs(lat) - 30)/bmavg + bavg30

        cmavg = 15.0/(cavg45-cavg30)
        cavg = (abs(lat) - 30)/cmavg + cavg30

    if (abs(lat) > 45 and abs(lat) <= 60):
        amavg = 15.0/(aavg60-aavg45)
        aavg = (abs(lat) - 45)/amavg + aavg45

        bmavg = 15.0/(bavg60-bavg45)
        bavg = (abs(lat) - 45)/bmavg + bavg45

        cmavg = 15.0/(cavg60-cavg45)
        cavg = (abs(lat) - 45)/cmavg + cavg45

    if (abs(lat) > 60 and abs(lat) <= 75):
        amavg = 15.0/(aavg75-aavg60)
        aavg = (abs(lat) - 60)/amavg + aavg60

        bmavg = 15.0/(bavg75-bavg60)
        bavg = (abs(lat) - 60)/bmavg + bavg60

        cmavg = 15.0/(cavg75-cavg60)
        cavg = (abs(lat) - 60)/cmavg + cavg60

    if (abs(lat) > 75):
        aavg = aavg75
        bavg = bavg75
        cavg = cavg75

    a = aavg
    b = bavg
    c = cavg

    el = math.sin(elev*math.pi/180)

    m = (1 + a/(1+b/(1+c)))/(el + a/(el+b/(el+c)))

    Mwet = m

    return(Mwet)

def writesac(velfile,site,stalat,stalon,doy,year,samprate,event):
    a = numpy.loadtxt(velfile)
    tind = a[:,0]
    gtime = a[:,1]
    leapsec = gpsleapsec(gtime[0])
    
    #Get the start time of the file in UTC
    date = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1)
    gpstime = (numpy.datetime64(date) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's')
    stime = (gtime[0]-leapsec)*numpy.timedelta64(1, 's')+ numpy.datetime64('1980-01-06T00:00:00')
    sitem = stime.item()
    print(sitem)
    styr = sitem.year
    stdy = sitem.day
    stmon = sitem.month
    sthr = sitem.hour
    stmin = sitem.minute
    stsec = sitem.second

    
    nv = a[:,2]
    ev = a[:,3]
    uv = a[:,4]
    sr = "{0:.2f}".format(float(samprate))
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXN.sac')
    headN = {'kstnm': site, 'kcmpnm': 'LXN', 'stla': float(stalat),'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}

    sacn = SACTrace(data=nv, **headN)
    sacn.write('output/' + event + '.' + site.upper() + '.' + sr + '.LXN.sac')
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXE.sac')

    headE = {'kstnm': site, 'kcmpnm': 'LXE', 'stla': float(stalat),'stlo': float(stalon),
         'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
         'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sace = SACTrace(data=ev, **headE)
    sace.write('output/' + event + '.' + site.upper() + '.' + sr + '.LXE.sac')
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXZ.sac')

    headZ = {'kstnm': site, 'kcmpnm': 'LXZ', 'stla': float(stalat),'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sacu = SACTrace(data=uv, **headZ)
    sacu.write('output/' + event + '.' + site.upper() + '.' + sr + '.LXZ.sac')

def lowpass(data,fcorner,fsample,order,zerophase=True):
    from scipy.signal import butter, filtfilt, lfilter


def filter_int_VTEC(ionfile,site,doy,year,event,samplerate,homedir):
    b = numpy.loadtxt(homedir+'loc_' + event + '.txt')
    [xs,ys,zs]=lla2ecef(float(b[1]),float(b[0]),0)
    a = numpy.loadtxt(ionfile)
    indx = a[:,0]
    snum = a[:,1]
    gpst = a[:,2]
    aTEC = a[:,3]
    vTEC = a[:,4]
    elev = a[:,5]
    azim = a[:,6]
    latpp = a[:,7]
    lonpp = a[:,8]
    dt = gpst - numpy.amin(gpst)
    uniquesats = numpy.unique(snum)
    outfile = 'output/ioncorr_' + site + '_' + doy + '_' + year + '_' + event + '.txt'
    ffo = open(outfile,'w')
    for i in range (0,len(uniquesats)):
        a1 = numpy.where(uniquesats[i] == snum)[0]
        difft = numpy.diff(gpst[a1])
        adiff = numpy.where(difft > 1/samplerate)[0]
        if (len(adiff) == 0 and len(a1) >= 240):
            pvTEC = numpy.polyfit(indx[a1],vTEC[a1],8)
            vTECmodel = numpy.polyval(pvTEC, indx[a1])
            vTECcorr = vTEC[a1] - vTECmodel
            vTEC_int = scipy.integrate.cumtrapz(vTECcorr,gpst[a1])
            bf, af = butter(4, [0.5/1000, 10/1000],btype='band',fs=samplerate)
            vTEC_int_filt = filtfilt(bf,af,vTEC_int-vTEC_int[0])
            
            paTEC = numpy.polyfit(indx[a1],aTEC[a1],8)
            aTECmodel = numpy.polyval(paTEC, indx[a1])
            aTECcorr = aTEC[a1] - aTECmodel
            aTECcorr_filt = filtfilt(bf,af,aTECcorr-aTECcorr[0])
            plotvTEC(event,site,"{0:.0f}".format(float(uniquesats[i])),gpst[a1],vTEC_int_filt,0) 
            for j in range (1,len(a1)):
                [xpp,ypp,zpp] = lla2ecef(float(latpp[a1[j]]),float(lonpp[a1[j]]),0)
                coord1 = (float(latpp[a1[j]]),float(lonpp[a1[j]]))
                coord2 = (float(b[1]),float(b[0]))

                distPP = geopy.distance.distance(coord1,coord2).km
                #distPP = math.sqrt(math.pow(xs-xpp,2)+math.pow(ys-ypp,2)+math.pow(zs-zpp,2))/1000
                dPP = "{0:.2f}".format(float(distPP))
                loPP = "{0:.4f}".format(float(lonpp[a1[j]]))
                laPP = "{0:.4f}".format(float(latpp[a1[j]]))
                vTECfilt = "{0:.5e}".format(float(vTEC_int_filt[j-1]))
                aTECfilt = "{0:.5e}".format(float(aTECcorr_filt[j]))
                vTEC2 = "{0:.5e}".format(float(vTEC_int[j-1]))
                aTEC2 = "{0:.5e}".format(float(aTECcorr[j]))
                el = "{0:.2f}".format(float(elev[a1[j]]))
                az = "{0:.2f}".format(float(azim[a1[j]]))
                t = "{0:.2f}".format(float(gpst[a1[j]]))
                sat = "{0:.0f}".format(float(snum[a1[j]]))            
                ffo.write(t+' '+sat+' '+vTEC2+' '+vTECfilt+' '+aTEC2+' '+aTECfilt+' '+laPP+' '+loPP+' '+az+' '+el+' '+dPP+'\n')
        elif (len(adiff) > 0):
            ADIFF = list()
            ADIFF.append(0)
            for j in range(0, len(adiff)):
                ADIFF.append(adiff[j])
            ADIFF.append(len(a1))
            print(ADIFF)
            for j in range(1, len(ADIFF)):
                print(j)
                VTEC = vTEC[a1[ADIFF[j-1]:ADIFF[j]]]
                ATEC = aTEC[a1[ADIFF[j-1]:ADIFF[j]]]
                INDX = indx[a1[ADIFF[j-1]:ADIFF[j]]]
                GPST = gpst[a1[ADIFF[j-1]:ADIFF[j]]]
                LATPP = latpp[a1[ADIFF[j-1]:ADIFF[j]]]
                LONPP = lonpp[a1[ADIFF[j-1]:ADIFF[j]]]
                ELEV = elev[a1[ADIFF[j-1]:ADIFF[j]]]
                AZIM = azim[a1[ADIFF[j-1]:ADIFF[j]]]
                if (len(INDX) >= 240):
                    pvTEC = numpy.polyfit(INDX,VTEC,8)
                    vTECmodel = numpy.polyval(pvTEC, INDX)
                    vTECcorr = VTEC - vTECmodel
                    vTEC_int = scipy.integrate.cumtrapz(vTECcorr,GPST)
                    bf, af = butter(4, [0.5/1000, 10/1000],btype='band',fs=samplerate)
                    vTEC_int_filt = filtfilt(bf,af,vTEC_int-vTEC_int[0])

                    paTEC = numpy.polyfit(INDX,ATEC,8)
                    aTECmodel = numpy.polyval(paTEC, INDX)
                    aTECcorr = ATEC - aTECmodel
                    aTECcorr_filt = filtfilt(bf,af,aTECcorr-aTECcorr[0])
                    plotvTEC(event,site,"{0:.0f}".format(float(uniquesats[i])),GPST,vTEC_int,j-1) 
                    for k in range (1,len(GPST)):
                        [xpp,ypp,zpp] = lla2ecef(float(LATPP[k]),float(LONPP[k]),0)
                        distPP = math.sqrt(math.pow(xs-xpp,2)+math.pow(ys-ypp,2)+math.pow(zs-zpp,2))/1000
                        dPP = "{0:.2f}".format(float(distPP))
                        loPP = "{0:.4f}".format(float(LONPP[k]))
                        laPP = "{0:.4f}".format(float(LATPP[k]))
                        vTECfilt = "{0:.5e}".format(float(vTEC_int_filt[k-1]))
                        aTECfilt = "{0:.5e}".format(float(aTECcorr_filt[k]))
                        vTEC2 = "{0:.5e}".format(float(vTEC_int[k-1]))
                        aTEC2 = "{0:.5e}".format(float(aTECcorr[k]))
                        el = "{0:.2f}".format(float(ELEV[k]))
                        az = "{0:.2f}".format(float(AZIM[k]))
                        t = "{0:.2f}".format(float(GPST[k]))
                        sat = "{0:.0f}".format(float(uniquesats[i]))            
                        ffo.write(t+' '+sat+' '+vTEC2+' '+vTECfilt+' '+aTEC2+' '+aTECfilt+' '+laPP+' '+loPP+' '+az+' '+el+' '+dPP+'\n')


    ffo.close()


def plotvTEC(event,site,sat,t,vTECfilt,ind):
    fname = 'output/fig/'+event+'/'+event+'_'+site+'_'+sat.zfill(2)+'_'+str(ind).zfill(2)+'.png'
    tplot = t[1:len(t)]
    print (fname)
    #A=plt.magnitude_spectrum(vTECfilt, Fs=1/30)
    #print(A[1])
    plt.plot(tplot,vTECfilt)
    plt.title(event+' '+site+' G'+sat.zfill(2))
    plt.ylabel('Filtered TEC')
    plt.xlabel('GPS Time')
    plt.savefig(fname)
    plt.clf()
    plt.close()
    return
