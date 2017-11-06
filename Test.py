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

#Test user input
waveband = str(input('Please choose your filter [B/V/R/u/g/r/i/z] = '))

star_ra = []
star_dec = []
star_mag = []
star_magerr = []

with open('skymapperCatalogue.csv','rb') as cat:
	next(cat)
	
    	reader = csv.reader(cat, delimiter=',')

    	for row in reader:
		if waveband == 'r':
			if row[66] != '':
				star_mag.append(row[66])
				star_magerr.append(row[67])
				star_ra.append(row[2])
				star_dec.append(row[3])
		if waveband == 'i':
			if row[70] != '':
				star_mag.append(row[70])
				star_magerr.append(row[71])
				star_ra.append(row[2])
				star_dec.append(row[3])
		if waveband == 'z':
			if row[74] != '':
				star_mag.append(row[74])
				star_magerr.append(row[75])
				star_ra.append(row[2])
				star_dec.append(row[3])

	print star_mag
	print star_magerr
	print star_ra
	print star_dec
		
