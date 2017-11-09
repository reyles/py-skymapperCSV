#!/usr/bin/python
# Photometry calibration script v1.0
# Scrapes skymapper CSV files
# Created by Adam Higgins & modified by Rob Eyles
# Email: abh13@le.ac.uk & raje1@leicester.ac.uk

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from xml.etree import ElementTree
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob
import csv

sys.path.append('/home/r/raje1/Documents/')
import requests
from bs4 import BeautifulSoup  

""" Calculates the magnitude zero-point offset for an image on non-photometric
nights. It currently does this comparing field star instrumental magnitudes
against measured catalogue values and determines the offset! The script
currently uses csvs from Skymapper. It can calibrate images in
normal filters (B,V,R) and SDSS filters (u,g,r,i,z) """

""" Ask for user inputs... """
directory = str(input('Please input the target directory = '))
targetra = float(input('Target Source RA = '))
targetdec = float(input('Target Source Dec = '))
arcmin = float(input('Search Radius (Arcmin) = '))
mag_zero = float(input('Your initial mag zero point = '))
waveband = str(input('Please choose your filter [B/V/R/u/g/r/i/z] = '))
if (waveband != 'V' and waveband != 'B' and waveband != 'R' and waveband != 'u' and waveband != 'g' and waveband != 'r' and waveband != 'i' and waveband != 'z'):
    print 'Script does not calibrate for this waveband!'
    sys.exit()

""" Run commands in bash for SExtractor """
os.chdir(directory)
os.system('sex new-image.fits -c photpars.txt -checkimage_type APERTURES,BACKGROUND -checkimage_name aper.fits,back.fits')

""" Takes results from SExtractor """
file = 'photresults.cat'
data = np.genfromtxt(file,skip_header=6,dtype=float)
rmag = data[:,0]
rmagerr = data[:,1]
ra = data[:,2]
dec = data[:,3]
xpixel = data[:,4]
ypixel = data[:,5]
  
""" Set of masks to remove bad detections/noisy edge pixels """
rmagm = np.ma.masked_where(rmag < 90, rmag)
rmag = rmag[rmagm.mask]
rmagerr = rmagerr[rmagm.mask]
ra = ra[rmagm.mask]
dec = dec[rmagm.mask]
xpixel = xpixel[rmagm.mask]
ypixel = ypixel[rmagm.mask]

"""Create variable lists"""
star_ra = []
star_dec = []
star_mag = []
star_magerr = []

"""Scrapes a skymapper csv file and adds them to the lists to compare"""
with open('skymapperCatalogue.csv','rb') as cat:
	next(cat)
	
    	reader = csv.reader(cat, delimiter=',')

    	for row in reader:
		if waveband == 'r':
			if row[66] != '':
				star_mag.append(float(row[66]))
				star_magerr.append(float(row[67]))
				star_ra.append(float(row[2]))
				star_dec.append(float(row[3]))
		if waveband == 'i':
			if row[70] != '':
				star_mag.append(float(row[70]))
				star_magerr.append(float(row[71]))
				star_ra.append(float(row[2]))
				star_dec.append(float(row[3]))
		if waveband == 'z':
			if row[74] != '':
				star_mag.append(float(row[74]))
				star_magerr.append(float(row[75]))
				star_ra.append(float(row[2]))
				star_dec.append(float(row[3]))

""" For catalog star mags with given error as 0, change this to 0.1 mag """
star_magerr = [0.1 if (x == 0.0 or x == -0.0) else x for x in star_magerr]

""" Crossmatch SExtractor detections output with nearby known stars
for calibration and print to file """
cm_error_region = 1
detect_ra = []
detect_dec = []
real_mag = []
real_magerr = []
fake_mag = []
fake_magerr = []

orig_stdout = sys.stdout
resultf = open('cm_results.txt', "w")
sys.stdout = resultf
print('DetectionRA   DetectionDec  CatalogRA   CatalogDec  MagReal  MagRealErr  MagMeas  MagMeasErr')
for i in range(len(rmag)):
    targetra_diff = (ra[i] - targetra)/0.000277778
    targetdec_diff = (dec[i] - targetdec)/0.000277778
    targetdiff = (targetra_diff**2 + targetdec_diff**2)**0.5
    if targetdiff > 5:
        for j in range(len(star_ra)):
            ra_diff = (ra[i] - star_ra[j])/0.000277778
            dec_diff = (dec[i] - star_dec[j])/0.000277778
            dist = (ra_diff**2 + dec_diff**2)**0.5
            if dist <= cm_error_region:
		if (ra[i] in detect_ra and dec[i] in detect_dec):
                    continue
                else:
                    detect_ra.append(ra[i])
                    detect_dec.append(dec[i])
                    real_mag.append(star_mag[j])
                    real_magerr.append(star_magerr[j])
                    fake_mag.append(rmag[i])
                    fake_magerr.append(rmagerr[i])

"""Removes very wrong points based on the difference in their matched magnitudes compared to the mean"""
"""for the population. There's probably a better way but it works."""
total_mag_diff = 0.0
mean_mag_diff = 0.0

for i in range(len(detect_ra)):
    total_mag_diff = total_mag_diff + abs(real_mag[i] - fake_mag[i])

mean_mag_diff = total_mag_diff / len(detect_ra)

k = 0
l = len(detect_ra)

while k < l:
    if abs(mean_mag_diff - abs(real_mag[k] - fake_mag[k])) > 1:
        detect_ra.pop(k)
	detect_dec.pop(k)
	real_mag.pop(k)
        real_magerr.pop(k)
        fake_mag.pop(k)
        fake_magerr.pop(k)
	l = l - 1
    else:
	print(detect_ra[k],' ',detect_dec[k],' ',ra[k],' ',dec[k],' ',real_mag[k],' ',real_magerr[k],' ',fake_mag[k],' ',fake_magerr[k])
	k = k + 1

sys.stdout = orig_stdout
resultf.close()

orig_stdout = sys.stdout
resultf = open('ds9_regions.reg', "w")
sys.stdout = resultf

for i in range(len(detect_ra)):
    print "WCS; circle ",detect_ra[i]," ",detect_dec[i]," ",0.0008

sys.stdout = orig_stdout
resultf.close()

""" Exit Programme if too few field stars to calibrate offset """
if len(real_mag) < 3:
    print('Not enough field stars to calibrate offset')
    sys.exit()

""" Plotting data and sigma clipping """
sigma = float(input('Significance level for data clipping = '))

def calib_plotting(xdata,ydata,xerrors,yerrors,sigma,figname):    
    """ Find calibration offset """
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Real (V) Mag')
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    ax1.errorbar(x=xdata,y=ydata,xerr=xerrors,yerr=yerrors,fmt='.')

    """ Just fitting a straight line for offset """
    def line1(x,c):
        return x + c

    def line2(x,c):
        return 0*x + c

    param, pcov = curve_fit(line1, xdata, ydata, sigma=yerrors)
    paramerr = np.sqrt(np.diag(pcov))

    """ Plot real mag vs instrumental mag """
    xx = np.linspace(min(xdata),max(xdata),1000)
    ymod1 = line1(xx,param)
    ymodhigh1 = line1(xx,param+sigma*paramerr)
    ymodlow1 = line1(xx,param-sigma*paramerr)
    ax1.plot(xx,ymod1,color='red')
    ax1.plot(xx,ymodhigh1,color='red',linestyle='--')
    ax1.plot(xx,ymodlow1,color='red',linestyle='--')

    """ Plot delta mag (offset) vs instrumental mag """
    delta = []
    deltaerrors = []
    for i in range(len(xdata)):
        delta.append(ydata[i] - xdata[i])
        deltaerrors.append(np.sqrt(xerrors[i]**2 + yerrors[i]**2))
    ymod2 = line2(xx,param)
    ymodhigh2 = line2(xx,param+sigma*paramerr)
    ymodlow2 = line2(xx,param-sigma*paramerr)
    ax2 = fig.add_subplot(212,sharex=ax1)
    ax2.errorbar(x=xdata,y=delta,xerr=xerrors,yerr=deltaerrors,fmt='.')
    ax2.plot(xx,ymod2,color='red')
    ax2.plot(xx,ymodhigh2,color='red',linestyle='--')
    ax2.plot(xx,ymodlow2,color='red',linestyle='--')
    ax2.set_ylabel('$\Delta$ Mag')
    ax2.set_xlabel('Instrumental (V) Mag')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.savefig(directory+figname)
    
    return param, paramerr

""" Plot raw data """
param, paramerr = calib_plotting(fake_mag,real_mag,fake_magerr,real_magerr,sigma,'calib.png')
print 'Magnitude Offset (no clipping) =',param,'+-',paramerr

def data_clipping(xdata,ydata,xerrors,yerrors,param,paramerr,sigma):
    """ Clips any data that is X-sigma or further away from the best fitting
    line """
    xdatac = []
    xerrorc = []
    ydatac = []
    yerrorc = []
    for i in range(len(xdata)):
        yhigh = xdata[i]+float(param)+sigma*float(paramerr)
        ylow = xdata[i]+float(param)-sigma*float(paramerr)
        if (ylow < ydata[i] < yhigh): 
            ydatac.append(ydata[i])
            yerrorc.append(yerrors[i])
            xdatac.append(xdata[i])
            xerrorc.append(xerrors[i])
        else:
            continue
    return xdatac, xerrorc, ydatac, yerrorc

""" Plot clipped data """
fake_magc, fake_magerrc, real_magc, real_magerrc = data_clipping(fake_mag,real_mag,fake_magerr,real_magerr,param,paramerr,sigma) 
if (len(fake_magc) == len(fake_mag) or len(fake_magc) < 3):
    print 'Fit cannot be further improved by clipping!'
    sys.exit()
paramc, paramerrc = calib_plotting(fake_magc,real_magc,fake_magerrc,real_magerrc,sigma,'calib_c.png')
print 'Magnitude Offset (one set of clipping) =',paramc,'+-',paramerrc
    
""" Plot doubly clipped data """
fake_magcc, fake_magerrcc, real_magcc, real_magerrcc = data_clipping(fake_magc,real_magc,fake_magerrc,real_magerrc,paramc,paramerrc,sigma) 
if (len(fake_magcc) == len(fake_magc) or len(fake_magcc) < 3):
    print 'Fit cannot be further improved by clipping twice!'
    sys.exit()
paramcc, paramerrcc = calib_plotting(fake_magcc,real_magcc,fake_magerrcc,real_magerrcc,sigma,'calib_dc.png')
print 'Magnitude Offset (two sets of clipping) =',paramcc,'+-',paramerrcc

""" Print the corrected zero point """
try:
    mag_zero_corr = mag_zero + paramc
except NameError:
    mag_zero_corr = mag_zero + param
print '\nCorrect zero point =',mag_zero_corr
