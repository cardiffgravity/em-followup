# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 11:08:12 2018

@author: Tilly
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.close('all')

def fit2png(filename):    
    fits_data = fits.getdata(filename)
    plt.imshow(fits_data, cmap='Reds', norm=LogNorm())
    imagename = filename.replace('.fits', '.png')
    plt.savefig(imagename)
    plt.close('all')
    pass

