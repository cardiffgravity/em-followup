#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 09:56:42 2018

@author: lewis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as cm
from scipy import ndimage
from scipy import signal
import scipy.signal as sg
from astropy.io import fits
from PIL import Image
import astropy.stats
import numpy.ma
import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D
from skimage import img_as_bool, io, color, morphology
from skimage.morphology import skeletonize
from skimage import measure
from numpy import random
from scipy import signal
from skimage import measure
import copy
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import ICRS, Galactic
from astropy import units as u
import matplotlib.ticker as mticker
import math

'''
--------------------------------------------------------------------
'''

#SECTION 1 - LOAD IN IMAGES AND EXTRACT COORDINATES

filename = get_pkg_data_filename('12_ogg2m001-fs02-20160321-0059-e90.fits')
hdu_list=fits.open(filename)
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
image_data1=hdu.data
hdu_list.close()

filename = get_pkg_data_filename('5_ogg2m001-fs02-20160309-0087-e90.fits')
hdu_list=fits.open(filename)
hdu = fits.open(filename)[0]
wcs2 = WCS(hdu.header)
image_data2=hdu.data
hdu_list.close()


#takes 2 images 
#convert pixels of i1 into coordinates 

def coords(data1,data2,system,system2): 
    x=data1.shape[1]
    y=data1.shape[0]
    grid=np.ones((y,x))
    points=np.where(grid==1)
    pointsy,pointsx=points[0],points[1]
    coords1=system.wcs_pix2world(pointsx,pointsy,1)
    coordsy,coordsx=coords1[1],coords1[0]
    c_join=(coordsy,coordsx)
    return c_join,points

coordinates1,points1 = coords(image_data1,image_data2,wcs,wcs2)
coordinates2,points2= coords(image_data2,image_data1,wcs2,wcs)

'''
--------------------------------------------------------------------
'''

#SECTION 2 - TRIM BOTH IMAGES SO SHOW SAME REGION OF SKY


#first recall the pixel coordinates - split x and y
points1y,points1x=points1[0],points1[1]
points2y,points2x=points2[0],points2[1]
#split coords into x and y arrays
ycord1,xcord1=coordinates1[0],coordinates1[1]
ycord2,xcord2=coordinates2[0],coordinates2[1]
#min and max values of image 1 
My1,my1=max(ycord1),min(ycord1)
Mx1,mx1=max(xcord1),min(xcord1)
#min and max values of image 2
My2,my2=max(ycord2),min(ycord2)
Mx2,mx2=max(xcord2),min(xcord2)
#set absolute coord values that niether image can overflow
MAXy=min(My1,My2)
MINy=max(my1,my2)
MAXx=min(Mx1,Mx2)
MINx=max(mx1,mx2)
#crop both images to these values 
#back into pixel coordinates
points1x_crop=points1x[(ycord1>MINy)
                   & (ycord1<MAXy)
                   & (xcord1>MINx)
                   & (xcord1<MAXx)]
points1y_crop=points1y[(ycord1>MINy)
                   & (ycord1<MAXy)
                   & (xcord1>MINx)
                   & (xcord1<MAXx)]
points2y_crop=points2y[(ycord2>MINy)
                   & (ycord2<MAXy)
                   & (xcord2>MINx)
                   & (xcord2<MAXx)]
points2x_crop=points2x[(ycord2>MINy)
                   & (ycord2<MAXy)
                   & (xcord2>MINx)
                   & (xcord2<MAXx)]
points1_new=points1y_crop,points1x_crop
points2_new=points2y_crop,points2x_crop

#turn pixel numbers into images 
#new size of images
#image 1 
no_ypixels1=len(np.where(points1x_crop==max(points1x_crop))[0])
no_xpixels1=len(np.where(points1y_crop==max(points1y_crop))[0])
image_data1_new=np.zeros((no_ypixels1,no_xpixels1))

#section of original images to keep
min_y1pixels,max_y1pixels=min(points1y_crop),max(points1y_crop)
min_x1pixels,max_x1pixels=min(points1x_crop),max(points1x_crop)

# and extract
image_data1_new[:,:]=image_data1[min_y1pixels:max_y1pixels,
                                   min_x1pixels:max_x1pixels]

#repeat previous 3 steps for second image
no_ypixels2=len(np.where(points2x_crop==max(points2x_crop))[0])
no_xpixels2=len(np.where(points2y_crop==max(points2y_crop))[0])
image_data2_new=np.zeros((no_ypixels2,no_xpixels2))
min_y2pixels,max_y2pixels=min(points2y_crop),max(points2y_crop)
min_x2pixels,max_x2pixels=min(points2x_crop),max(points2x_crop)
image_data2_new[:,:]=image_data2[min_y2pixels:max_y2pixels,
                                   min_x2pixels:max_x2pixels]

#now we have 2 images showing the same patch of sky!
'''
--------------------------------------------------------------------
'''
#SECTION 3 - MAKE IMAGES SAME SIZE (PIXEL SIZE)







'''
--------------------------------------------------------------------
'''
#SECTION 4 - COMPARE PIXEL VALUES

def convert(coords1y,coords1x,wcs2):
    points_converted=system2.wcs_world2pix(coordsx,coordsy,1)
    points2y,points2x=np.int_(points_converted[1]),np.int_(points_converted[0])
    points_converted=(points2y,points2x)
    return points_converted





'''
#y1,x1=image_data1.shape
#y2,x2=image_data2.shape
#
#co1y=coordinates1[0]
#co1x=coordinates1[1]
#co2y=coordinates2[0]
#co2x=coordinates2[1]
#if y1*x1 > y2*x2:
#    
#    minx,maxx=coordinates2[1].min(),coordinates2[1].max()
#    miny,maxy=coordinates2[0].min(),coordinates2[0].max()
#    co1x_new=co1x[np.where((co1x>minx) & (co1x<maxx) & (co1y>miny) & (co1y<maxy))]
#    co1y_new=co1y[np.where((co1x>minx) & (co1x<maxx) & (co1y>miny) & (co1y<maxy))]
#    places=np.where(coordinates1==(co1y_new,co1x_new))
#    
#else:
#    image_new=np.zeros_like(image_data2)
#    minx,maxx=coordinates1[1].min(),coordinates1[1].max()
#    miny,maxy=coordinates1[0].min(),coordinates1[0].max() 
#    
#    placex=np.where((co2x>minx) & (co2y<maxx) & (co2y>miny) & (co2y<maxy))
#    placey=np.where((co2x>minx) & (co2y<maxx) & (co2y>miny) & (co2y<maxy))
#    place=(placey,placex)
#    new=points2[0][placey],points2[1][placex]
#    
#    sizexmin,sizexmax=min(new[1]),max(new[1])
#    sizex=sizexmax-sizexmin
#    sizeymin,sizeymax=min(new[0]),max(new[0])
#    sizey=sizeymax-sizeymin
#    image_new=np.zeros((sizey,sizex))
#    image_new[new]=image_data2[new]
#    
#    co2x_new=co2x[placex]
#    co2y_new=co2y[placey]
#    
#
#    
#def crop(co1,co2,axis):
#    if axis=='y':
#        coord1=co1[0]
#        coord2=co2[0]
#    if axis=='x':
#        coord1=co1[1]
#        coord2=co2[2]
#    if max(coord1)>max(coord2):
#        coord1=coord1[coord1]
        
    
'''

