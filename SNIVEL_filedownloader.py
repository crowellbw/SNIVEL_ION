#!/usr/bin/env python
import wget
import os
import urllib
import SNIVEL_tools
#####################################################################################
#SNIVEL_filedownloader.py
#This code will download any RINEX, nav or UNR
#time series file requested
#Written by Brendan Crowell, University of Washington
#Last edited January 8, 2021
#Broadcast navigation messages are only downloaded from CDDIS
#RINEX files will try to download from UNAVCO, then CWU, then CDDIS, then SOPAC
#Time Series files will only download from UNR, cartesian positions
#VARIABLES
#year - 4 digit string of year
#doy - 3 digit string of day of year
#site - 4 digit string of site id
#####################################################################################
#This subroutine downloads the broadcast navigation message for a given day from CDDIS
def getbcorbit(year, doy):
    if not os.path.exists('nav'): #if nav folder doesn't exist, make it
        os.makedirs('nav')
    fnameZ = 'nav/brdc' + doy + '0.' +  year[-2:] + 'n.Z'
    fnamegz = 'nav/brdc' + doy + '0.' +  year[-2:] + 'n.gz'
    fname2 = 'nav/brdc' + doy + '0.' +  year[-2:] + 'n'
    if (os.path.isfile(fname2) == True):
        print ('Navigation file ' + fname2 + ' already exists')
    else:
        if (int(year) > 2020):
            ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
            ftps.login(user='anonymous', passwd='snivel@uw.edu')
            ftps.prot_p()
            ftps.cwd('gnss/data/daily/' + year + '/' + doy + '/' + year[-2:] + 'n/')
            ftps.retrbinary("RETR " + 'brdc' + doy + '0.' +  year[-2:] + 'n.gz', open('brdc' + doy + '0.' +  year[-2:] + 'n.gz', 'wb').write)
            os.system('mv '+ 'brdc' + doy + '0.' +  year[-2:] + 'n.gz  nav')
            os.system('gunzip' + ' ' + fnamegz)
        elif (int(year) == 2020 and int(doy) >334):
            ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
            ftps.login(user='anonymous', passwd='snivel@uw.edu')
            ftps.prot_p()
            ftps.cwd('gnss/data/daily/' + year + '/' + doy + '/' + year[-2:] + 'n/')
            ftps.retrbinary("RETR " + 'brdc' + doy + '0.' +  year[-2:] + 'n.gz', open('brdc' + doy + '0.' +  year[-2:] + 'n.gz', 'wb').write)
            os.system('mv '+ 'brdc' + doy + '0.' +  year[-2:] + 'n.gz  nav')
            os.system('gunzip' + ' ' + fnamegz)
        else:
            ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
            ftps.login(user='anonymous', passwd='snivel@uw.edu')
            ftps.prot_p()
            ftps.cwd('gnss/data/daily/' + year + '/' + doy + '/' + year[-2:] + 'n/')
            ftps.retrbinary("RETR " + 'brdc' + doy + '0.' +  year[-2:] + 'n.Z', open('brdc' + doy + '0.' +  year[-2:] + 'n.Z', 'wb').write)
            os.system('mv '+ 'brdc' + doy + '0.' +  year[-2:] + 'n.Z nav')
            os.system('gunzip' + ' ' + fnameZ)

def getrinexNGS(site,year,doy):
    if not os.path.exists('rinex'): #if rinex folder doesn't exist, make it
        os.makedirs('rinex')
    fnamegz = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.gz'
    fnameogz = 'rinex/' + site + doy + '0.' +  year[-2:] + 'o.gz'
    fnameZ = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.Z'
    fnamebz2 = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.bz2'
    fnamed = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd'
    fnameo = 'rinex/' + site + doy + '0.' +  year[-2:] + 'o'
    if (os.path.isfile(fnameo) == True):
        print ('Rinex file ' + fnameo + ' already exists')
    else:
        try:
            url = 'https://geodesy.noaa.gov/corsdata/rinex/' + year + '/' + doy + '/' + site + '/' + site + doy + '0.' +  year[-2:] + 'd.gz'
            print ('Attempting to download ' + fnamed + ' from NGS')
            wget.download(url, out='rinex/')
            os.system('gunzip' + ' ' + fnamegz)
            os.system('./crx2rnx' + ' ' + fnamed)
            os.remove(fnamed)
        except Exception:
            pass

#This subroutine will download RINEX files given the station, year and day of year. 
def getrinex(site, year, doy):
    if not os.path.exists('rinex'): #if rinex folder doesn't exist, make it
        os.makedirs('rinex')
    fnamegz = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.gz'
    fnameZ = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.Z'
    fnamebz2 = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd.bz2'
    fnamed = 'rinex/' + site + doy + '0.' +  year[-2:] + 'd'
    fnameo = 'rinex/' + site + doy + '0.' +  year[-2:] + 'o'
    if (os.path.isfile(fnameo) == True):       
        print ('Rinex file ' + fnameo + ' already exists')
    else:
        try:
            url = 'ftp://data-out.unavco.org/pub/rinex/obs/' + year + '/' + doy + '/' + site + doy + '0.' +  year[-2:] + 'd.Z'
            print ('Attempting to download ' + fnamed + ' from UNAVCO')
            wget.download(url, out='rinex/')
            os.system('gunzip' + ' ' + fnameZ)
            os.system('./crx2rnx' + ' ' + fnamed)
            os.remove(fnamed)
        except Exception:
            print ('File not at UNAVCO, checking CWU')
            try:
                url = 'https://www.geodesy.cwu.edu/data_ftp_pub/data/' + year+ '/' + doy + '/30sec/' + site + doy + '0.' +  year[-2:] + 'd.bz2'
                print ('Attempting to download ' + fnamed + ' from CWU')
                wget.download(url, out='rinex/')
                os.system('bzip2 -d' + ' ' + fnamebz2)
                os.system('./crx2rnx' + ' ' + fnamed)
                os.remove(fnamed)
            except Exception:
                print ('File not at CWU, checking CSN')
                try:
                    url = 'http://gps.csn.uchile.cl/data/' + year+ '/' + doy + '/' + site + doy + '0.' +  year[-2:] + 'd.Z'
                    print ('Attempting to download ' + fnamed + ' from CSN')
                    wget.download(url, out='rinex/')
                    os.system('gunzip' + ' ' + fnameZ)
                    os.system('./crx2rnx' + ' ' + fnamed)
                    os.remove(fnamed)
                except Exception:
                    print ('File not at CDDIS, checking SOPAC')
                    try:
                        url = 'ftp://garner.ucsd.edu/pub/rinex/' + year+ '/' + doy + '/' + site + doy + '0.' +  year[-2:] + 'd.Z'
                        print ('Attempting to download ' + fnamed + ' from SOPAC')
                        wget.download(url, out='rinex/')
                        os.system('gunzip' + ' ' + fnameZ)
                        os.system('./crx2rnx' + ' ' + fnamed)
                        os.remove(fnamed)
                    except Exception:
                        try:
                            url = 'ftp://ftp.earthobservatory.sg/SugarData/' + year+ '/' + site + '/RINEX2-HATANAKA/'  + site + doy + '0.' +  year[-2:] + 'd.Z'
                            print ('Attempting to download ' + fnamed + ' from EOS')
                            wget.download(url, out='rinex/')
                            os.system('gunzip' + ' ' + fnameZ)
                            os.system('./crx2rnx' + ' ' + fnamed)
                            os.remove(fnamed)
                        except Exception:
                            try:
                                url = 'https://geodesy.noaa.gov/corsdata/rinex/' + year + '/' + doy + '/' + site + '/' + site + doy + '0.' +  year[-2:] + 'd.gz'
                                print ('Attempting to download ' + fnamed + ' from NGS')
                                wget.download(url, out='rinex/')
                                os.system('gunzip' + ' ' + fnamegz)
                                os.system('./crx2rnx' + ' ' + fnamed)
                                os.remove(fnamed)
                            except Exception:
                                pass
                            print ('File not found, moving onto next station')


#This subroutine will download highrate (1-Hz) RINEX files
def getrinexhr(site, year, doy):
    if not os.path.exists('rinex_hr'): #if rinex highrate folder doesn't exist, make it
        os.makedirs('rinex_hr')
    fnameZ = 'rinex_hr/' + site + doy + '0.' +  year[-2:] + 'd.Z'
    fnamebz2 = 'rinex_hr/' + site + doy + 'i.' +  year[-2:] + 'd.bz2'
    fnamecwuo = 'rinex_hr/' + site + doy + 'i.' +  year[-2:] + 'o'
    fnamecwud = 'rinex_hr/' + site + doy + 'i.' +  year[-2:] + 'd'
    fnamed = 'rinex_hr/' + site + doy + '0.' +  year[-2:] + 'd'
    fnameo = 'rinex_hr/' + site + doy + '0.' +  year[-2:] + 'o'
    if (os.path.isfile(fnameo) == True):       
        print ('Rinex file ' + fnameo + ' already exists')
    else:
        try:
            url = 'ftp://data-out.unavco.org/pub/highrate/1-Hz/rinex/' + year + '/' + doy + '/' + site + '/' + site + doy + '0.' +  year[-2:] + 'd.Z'
            print ('Attempting to download ' + fnamed + ' from UNAVCO')
            wget.download(url, out='rinex_hr/')
            os.system('gunzip' + ' ' + fnameZ)
            os.system('./crx2rnx' + ' ' + fnamed)
            os.remove(fnamed)
        except Exception:
            print ('File not at UNAVCO checking CWU')
            try:
                url = 'https://www.geodesy.cwu.edu/data_ftp_pub/data/' + year+ '/' + doy + '/01sec/' + site + doy + 'i.' +  year[-2:] + 'd.bz2'
                print ('Attempting to download ' + fnamed + ' from CWU')
                wget.download(url, out='rinex_hr/')
                os.system('bzip2 -d' + ' ' + fnamebz2)
                os.system('./crx2rnx' + ' ' + fnamecwud)
                os.rename(fnamecwuo, fnameo)
                os.remove(fnamecwud)
            except Exception:
                print ('File not at CWU, moving on')



#Examples
##getrinexhr('lwck','2018','002')
##getrinex('p494','2018','002')
##getbcorbit('2018','002')
##gettseries('p494')
##
#getsp3file('2018','002')


