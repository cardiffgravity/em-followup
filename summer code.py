#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 09:56:42 2018

@author: lewis
"""
from __future__ import division
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
from matplotlib.colors import LogNorm
from reproject import reproject_interp
import os


'''
--------------------------------------------------------------------
'''

#SECTION 1 - LOAD IN IMAGES AND EXTRACT COORDINATES



def getimage(keyword):

    path1='/Users/lewisprole/Documents/University/year3/summer_project/images/%s/%s_1.fits'%(keyword,keyword)
    hdu_list1=fits.open(path1)
    hdu1 = fits.open(path1)[0]
    wcs1 = WCS(hdu1.header)
    image_data1=hdu1.data
    hdu_list1.close()
    
    path2='/Users/lewisprole/Documents/University/year3/summer_project/images/%s/%s_2.fits'%(keyword,keyword)
    hdu_list2=fits.open(path2)
    hdu2 = fits.open(path2)[0]
    wcs2 = WCS(hdu2.header)
    image_data2=hdu2.data
    hdu_list2.close()
    
    return hdu1,wcs1,image_data1,hdu2,wcs2,image_data2

#function which extracts coordinates 
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

def save_jpeg(name,Vmin,Vmax):    
    fits_data = fits.getdata(name)
    plt.figure(),plt.imshow(fits_data,vmin=Vmin,vmax=Vmax, cmap='Reds')
    imagename = name.replace('.fits', '.png')
    plt.savefig(imagename)
    plt.close('all')
    

path_images='/Users/lewisprole/Documents/University/year3/summer_project/images'
for filename in os.listdir(path_images):
    if filename=='.DS_Store':
        print('lol')
    else:
        hdu1,wcs1,image_data1,hdu2,wcs2,image_data2=getimage(filename)
        
        coordinates1,points1 = coords(image_data1,image_data2,wcs1,wcs2)
        coordinates2,points2 = coords(image_data2,image_data1,wcs2,wcs1)
        
        xrange1 = max(coordinates1[1])-min(coordinates1[1])
        xrange2 = max(coordinates2[1])-min(coordinates2[1])
        if xrange1>xrange2:
            new_image_data1, footprint = reproject_interp(hdu1, hdu2.header)
            new_image_data2=image_data2
        else:
            new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
            new_image_data1=image_data1
            
        med1=np.median(new_image_data1)
        min1=new_image_data1.min()
        max1=new_image_data1.max()
        range1=max1-min1
        
        med2=np.median(new_image_data2)
        min2=new_image_data2.min()
        max2=new_image_data2.max()
        range2=max2-min2
        

#        save_jpeg(new_image_data1,med1-range1/10,med1+range1/10,filename,1)
#        save_jpeg(new_image_data2,med1-range2/10,med2+range2/10,filename,1)
#        plt.figure(),plt.imshow(new_image_data1,vmin=med1-range1/10,vmax=med1+range1/10,cmap='Reds')
#        plt.figure(),plt.imshow(new_image_data2,vmin=med2-range2/10,vmax=med2+range2/10,cmap='Reds')
        subtraction=new_image_data1-new_image_data2
        hdu = fits.PrimaryHDU(new_image_data1)
        hdul = fits.HDUList([hdu])
        hdul.writeto('%s_1.fits'%filename)
        save_jpeg('%s_1.fits'%filename,med1-range1/10,med1+range1)
        os.remove('%s_1.fits'%filename)
        
        hdu = fits.PrimaryHDU(new_image_data2)
        hdul = fits.HDUList([hdu])
        hdul.writeto('%s_2.fits'%filename)
        save_jpeg('%s_2.fits'%filename,med2-range2/10,med2+range2)
        os.remove('%s_2.fits'%filename)

'''
#set up separate grids for x and y coordinates
#2D array giving the x/y coordinates in terms of pixels 
pixelsy1,pixelsx1 = np.zeros_like(image_data1),np.zeros_like(image_data1)
pixelsy1[points1]=points1[0]
pixelsx1[points1]=points1[1]
pixelsy2,pixelsx2 = np.zeros_like(image_data2),np.zeros_like(image_data2)
pixelsy2[points2]=points2[0]
pixelsx2[points2]=points2[1]

#2D array giving coordinates in terms of ICRS
ygrid1,xgrid1=np.zeros_like(image_data1),np.zeros_like(image_data1)
ygrid1[points1]=coordinates1[0]
xgrid1[points1]=coordinates1[1]

ygrid2,xgrid2=np.zeros_like(image_data2),np.zeros_like(image_data2)
ygrid2[points2]=coordinates2[0]
xgrid2[points2]=coordinates2[1]
'''
'''
--------------------------------------------------------------------
'''
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

#set up mask showing area you'd like to keep
mask1=np.zeros_like(xgrid1)
mask1[(xgrid1>MINx) & (xgrid1<MAXx) & (ygrid1>MINy) & (ygrid1<MAXy)]=1
rangex1=len(xgrid1[1])
rangey1=len(xgrid1[0])



mask2=np.zeros_like(xgrid2)
mask2[(xgrid2>MINx) & (xgrid2<MAXx) & (ygrid2>MINy) & (ygrid2<MAXy)]=1
rangex2=len(xgrid2[1])
rangey2=len(xgrid2[0])

#not a perfect rectangle of 2d data, need to clean up the edges
#find x size of square wanted 
lengthsx1=np.zeros(rangex1)
for i in range (rangex1):
    lengthsx1[i]=sum(mask1[:,i])
newsize_x1=max(lengthsx1)
#repeat for y measure
lengthsy1=np.zeros(rangey1)
for i in range (rangey1):
    lengthsy1[i]=sum(mask1[:,i])
newsize_y1=max(lengthsy1)
#create cropped image space
image_data1_crop=np.zeros((int(newsize_y1),int(newsize_x1)))

#repeat for image 2
#find x size of square wanted 
lengthsx2=np.zeros(rangex2)
for i in range (rangex2):
    lengthsx2[i]=sum(mask2[:,i])
newsize_x2=max(lengthsx2)
#repeat for y measure
lengthsy2=np.zeros(rangey2)
for i in range (rangey2):
    lengthsy2[i]=sum(mask2[i,:])
newsize_y2=max(lengthsy2)
#create cropped image space
image_data2_crop=np.zeros((int(newsize_y2),int(newsize_x2)))
crop_cord2=np.where(image_data2_crop==0)

#so now we have:
#a) grid of coordinates for whole image
#b) grid of pixel coordinates for whole image
#c) crop mask for whole image 
#d) new cropped image space
#e) arrays describing sum of mask in x/y directions 
#use e) to create grid of coordinates to be kept 
keep_mask2=np.ones_like(image_data2)
locx2=np.where(lengthsx2==0)
keep_mask2[:,locx2]=0
newsize_x2=len(lengthsx2[np.where(lengthsx2>0)])
locy2=np.where(lengthsy2==0)
keep_mask2[locy2,:]=0
newsize_y2=len(lengthsy2[np.where(lengthsy2>0)])
image_data2_crop=np.zeros((int(newsize_y2),int(newsize_x2)))
crop_cord2=np.where(image_data2_crop==0)
image_data2_crop[crop_cord2]=image_data2[np.where(keep_mask2==1)]

keep_mask1=np.ones_like(image_data1)
locx1=np.where(lengthsx1==0)
keep_mask1[:,locx1]=0
newsize_x1=len(lengthsx1[np.where(lengthsx1>0)])
locy1=np.where(lengthsy1==0)
keep_mask1[locy1,:]=0
newsize_y1=len(lengthsy1[np.where(lengthsy1>0)])
image_data1_crop=np.zeros((int(newsize_y1),int(newsize_x1)))
crop_cord1=np.where(image_data1_crop==0)
image_data1_crop[crop_cord1]=image_data1[np.where(keep_mask1==1)]

ax = plt.subplot(projection=wcs2)
ax.imshow(keep_mask2)
#for j in range (yrange1):
#    sy1=sum(ymask1[j,:])


    



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
no_ypixels1=len(np.where(points1x_crop==(np.median(points1x_crop)))[0])
no_xpixels1=len(np.where(points1y_crop==(np.median(points1y_crop)))[0])
image_data1_new=np.zeros((no_ypixels1,no_xpixels1))

#section of original images to keep
min_y1pixels,max_y1pixels=min(points1y_crop),max(points1y_crop)
min_x1pixels,max_x1pixels=min(points1x_crop),max(points1x_crop)

# and extract

image_data1_new[:,:]=image_data1[min_y1pixels:max_y1pixels+1,
                                   min_x1pixels:max_x1pixels+1]

#repeat previous 3 steps for second image
no_ypixels2=len(np.where(points2x_crop==(np.median(points2x_crop)))[0])
no_xpixels2=len(np.where(points2y_crop==(np.median(points2y_crop)))[0])
image_data2_new=np.zeros((no_ypixels2,no_xpixels2))
min_y2pixels,max_y2pixels=min(points2y_crop),max(points2y_crop)
min_x2pixels,max_x2pixels=min(points2x_crop),max(points2x_crop)
image_data2_new[:,:]=image_data2[min_y2pixels:max_y2pixels+1,
                                   min_x2pixels:max_x2pixels+1]

#now we have 2 images showing the same patch of sky!
subtraction=image_data1_new-image_data2_new
'''
'''
--------------------------------------------------------------------
'''
'''
#SECTION 3 - ROTATE IMAGE & plot
#pip install reproject
#from reproject import reporject_interp
#array,footprint=reproject_interp(hdu2,hud.header)
#ax1=plt.subplot(1,2,1,projection=WCS(hdu.header))




'''
'''
--------------------------------------------------------------------
'''
'''
#SECTION 4 - COMPARE PIXEL VALUES


#number of pixels in images 
N1=len(image_data1_crop[:,0])*len(image_data1_crop[0,:])
N2=len(image_data2_crop[:,0])*len(image_data2_crop[0,:])
ny1,nx1=len(image_data1_crop[:,0]),len(image_data1_crop[0,:])
ny2,nx2=len(image_data2_crop[:,0]),len(image_data2_crop[0,:])
if N1 > N2:
    ratioy=ny1/ny2
    ratiox=nx1/nx2
    image_data2_new_resize=ndimage.zoom(image_data2_crop,(ratioy,ratiox),order=0)
    image_data1_new_resize=image_data1_crop
if N2 > N1:
    ratioy=ny2/ny1
    ratiox=nx2/nx1
    image_data1_new_resize=ndimage.zoom(image_data1_crop,(ratioy,ratiox),order=0)
    image_data2_new_resize=image_data2_crop
    
subtraction=image_data2_new_resize-image_data1_new_resize
'''
'''
--------------------------------------------------------------------
'''
'''
#SECTION - SAVE AS JPEG


def fit2png(filename):    
    fits_data = fits.getdata(filename)
    plt.imshow(fits_data, cmap='Reds', norm=LogNorm())
    imagename = filename.replace('.fits', '.png')
    plt.savefig(imagename)
    plt.close('all')
    pass
'''
