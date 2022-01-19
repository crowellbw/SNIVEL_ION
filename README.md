# SNIVEL_ION

Follow the directions within SNIVEL for running the code. You will need an additional file, loc_event.txt, that defines the location of the source so distance to ionospheric piercing points can be determined. The format for this file is a single line with longitude and latitude. See the example I uploaded for Calbuco.

There will be two output files named the following:

ion_site_doy_year.txt

ioncorr_site_doy_year_event.txt


The ion file has the following columns:

index time, satellite number, GPS time, absolute TEC, variometric TEC, elevation angle, azimuth, latitude piercing point, longitude piercing point

The ioncorr file has the following columns:

GPS time, satellite, variometric TEC, filtered variometric TEC, absolute TEC, filtered absolute TEC, latitude piercing point, longitude piercing point, azimuth, elevation angle, distance to piercing point from source

The filtered time series use a 2 to 10 mHz bandpass filter.

