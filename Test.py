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

with open('skymapperCatalogue.csv','rb') as cat:
	reader = csv.reader(cat)
	for row in reader:
		print row
