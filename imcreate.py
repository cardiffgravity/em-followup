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
import matplotlib.colors as colors
from reproject import reproject_interp
import os
import csv
import shutil
import time
import threading
import skimage.draw

'''
--------------------------------------------------------------------

REQUIREMENTS
create a general folder that contains the following empty folders:
    
'images'
'LCO_images'
'fits_images'
'subtraction'
'relative'
'Zooniverse_upload'
'rejects'
'finished_images'

in 'images', create a folder per target to store target images
2 comparison images of single target per target folder
images must have same name keyword as folder
e.g. target folder called 'keyword', images within called 'keyword_1', 'keyword_2'
name target folders using actual target name 

to use: example
path_general='/Users/lewisprole/Documents/University/year3/summer_project'
CPR(path_general)
bright_diff(path_general)
sub(path_general)
out_save(path_general,'yes')

this is all that is neccessary

--------------------------------------------------------------------
'''

def remove_space(path_general,dir1):
    '''function that will remove any spaces in the target folder name
    input: path to general folder, name of images folder e.g 'images'
    '''
    print('remove_space')
    path=path_general+'/'+dir1
    for filename in sorted(os.listdir(path)):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')   
        else:
            
            name_new=filename.replace(' ','_')
            os.rename(path+'/%s'%(filename) , path+'/%s'%(name_new))
            
            

def combine_fold(path_general):
    '''function to move target folders from the LCO downloads folder to the main 
    images folder
    input: path to general folder
    '''
    print('combine_fold')
    files = os.listdir(path_general+'/LCO_images/raw')
    dest=path_general+'/images'
    for f in files:
        if f=='.DS_Store': #problematic folder that shows up sometimes
            print('.')   
        else:
            if os.path.exists(dest+'/'+f)==True:
                shutil.move(dest+'/'+f, path_general+'/rejects')
            shutil.move(path_general+'/LCO_images/raw/'+f, dest)
        

def reject_dir(path_general):
    '''function that checks if there are at least 2 imags per target
    if not the target folder is sent to reject folder
    input: path to general folder
    '''
    print('reject_dir')
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')   
        else:
            n=len(sorted(os.listdir(path_general+'/images/%s'%filename)))
            dest=path_general+'/rejects'
            if n<2:
                if os.path.exists(dest+'/'+filename)==True:
                    shutil.rmtree(dest+'/'+filename,ignore_errors=True)
                shutil.move(path_general+'/images/%s'%filename, dest)
    

def getimage(path_general,keyword):
    '''function that grabs fits image infomation from .fits file 
    located in the 'images' folder
    input: path to general folder, name of image target folder
    returns:hdu,wcs & image data of both images
    '''

    path=path_general+'/images/%s/%s_1.fits'%(keyword,keyword)
    hdu_list1=fits.open(path)
    hdu1 = fits.open(path)[0]
    wcs1 = WCS(hdu1.header)
    image_data1=hdu1.data
    hdu_list1.close()
    
    path=path_general+'/images/%s/%s_2.fits'%(keyword,keyword)
    hdu_list2=fits.open(path)
    hdu2 = fits.open(path)[0]
    wcs2 = WCS(hdu2.header)
    image_data2=hdu2.data
    hdu_list2.close()
    
    return hdu1,wcs1,image_data1,hdu2,wcs2,image_data2



#similar function used later for extracting new fits images 
def get_fitsimage(path_general,keyword):
    '''same function as getimage but taking images from
    the 'fits_images' folder instead of 'images'
    input: path to general folder, name of image target folder
    returns:hdu,wcs & image data of both images
    '''

    path=path_general+'/fits_images/%s_1.fits'%(keyword)
    hdu1 = fits.open(path)[0]
    wcs1 = WCS(hdu1.header)
    image_data1=hdu1.data
    
    path=path_general+'/fits_images/%s_2.fits'%(keyword)
    hdu2 = fits.open(path)[0]
    wcs2 = WCS(hdu2.header)
    image_data2=hdu2.data

    return hdu1,wcs1,image_data1,hdu2,wcs2,image_data2


    
def coords(data1,data2,system,system2): 
    '''function to grab galactic coordinates of image pixels
    input:image data from image1, same for image 2, wcs of image1, wcs of image2
    output:list of (ra,dec) coords, minimum ra, maximum ra, minimum dec, maximum dec
    '''
    x=data1.shape[1]
    y=data1.shape[0]
    grid=np.ones((y,x))
    points=np.where(grid==1)
    pointsy,pointsx=points[0],points[1]
    coords1=system.wcs_pix2world(pointsx,pointsy,1)
    coordsy,coordsx=coords1[1],coords1[0]
    c_join=(coordsy,coordsx)
    min_ra,max_ra=min(coordsy),max(coordsy)
    min_dec,max_dec=min(coordsx),max(coordsx)

    return c_join,points,min_ra,max_ra,min_dec,max_dec


def fz_remove(path_general):
    '''function which will uncompresses .fits.fz files into .fits files
    will check every target file in 'images' folder
    input:path to general folder
    '''
    print('fz_remove')
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            path=path_general+'/images/%s'%filename
            for name in os.listdir(path):
                if name.endswith('.fz')==True:
                    image=astropy.io.fits.getdata(path+'/'+name)
                    header=fits.getheader(path+'/'+name,ext=1)
                    hdu = fits.PrimaryHDU(image,header=header)
                    hdul = fits.HDUList([hdu])
                    hdul.writeto(path+'/%s'%name[:-3])   
                    os.remove(path+'/%s'%name)
        pass
             
    
      
def image_rename(path_general):
    '''function that will go into each target file in 'images' folder 
    and change the names of the images based on the target folder name
    e.g. target foler named 'M31',images renames 'M31_1.fits' & 'M31_2.fits'
    input: path to general folder
    '''
    print('image_rename')
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            path=path_general+'/images/%s'%filename
            n=len(os.listdir(path))
            names=os.listdir(path)
            
                
            try:
                names.remove('.DS_Store')
                n=n-1
            except: 
                pass
            for i in range (n):
                os.rename(path+'/%s'%(names[i]),path+'/%s%i.fits'%(filename,i+1))
            for i in range (n):
                os.rename(path+'/%s%i.fits'%(filename,i+1),path+'/%s_%i.fits'%(filename,i+1))
    


def save_png(name,Vmin,Vmax,path_general,wcs,types,contour):
    '''function that will save the .fits images as .png images 
    input: target name, min and max displayed values, path to general folder,
    wcs of image, subtraction image of target, relative subtraction image of target,
    title of figure, any contour data
    '''
    n=name.find('_')
    t_name=name[:n]
    if types=='subtraction':
        fits_data = fits.getdata(path_general+'/subtraction/%s'%name)
        c='viridis'
        title='%s: subtraction image'%t_name
    if types=='relative':
        fits_data = fits.getdata(path_general+'/relative/%s'%name)
        c='plasma'
        title='%s: percentage change'%t_name
    if types=='normal':
        fits_data = fits.getdata(path_general+'/fits_images/%s'%name)
        c='magma'
        n=name.find('.')
        t_name=name[:n]
        title='%s'%t_name
        
    plt.subplot(projection=wcs)
    plt.imshow(fits_data,vmin=Vmin,vmax=Vmax,cmap=c,norm=LogNorm())
    plt.grid(color='white', ls='solid')
    if Vmax-Vmin>=1:            
        plt.colorbar(label="log scale")
    else:
        plt.colorbar(ticks=[Vmin,Vmax],label="log scale")
    plt.title(title)
    imagename = name.replace('.fits', '.png')
    plt.savefig(path_general+'/Zooniverse_upload/%s'%imagename)
    plt.close('all')
    



 
def manifest(name,path_general):
    '''function to create manifest to be uploaded to zooniverse
    .csv telling zooniverse which images belong to the same subject set,
    saves it in the Zoonivere_upload folder
    input:target name, path to general folder
    '''
    if os.path.exists(path_general+"/Zooniverse_upload/Manifest.csv")==False:
        with open(path_general+"/Zooniverse_upload/Manifest.csv", "w+") as csvfile:
            write = csv.writer(csvfile)
            write.writerow(['image1','image2','subtraction','relative'])
            write.writerow(['%s_1.png'%name,'%s_2.png'%name,'%s_sub.png'%name,'%s_rel.png'%name])
    
    else:
        with open(path_general+"/Zooniverse_upload/Manifest.csv", "a") as csvfile:
            write = csv.writer(csvfile)
            write.writerow(['%s_1.png'%name,'%s_2.png'%name,'%s_sub.png'%name,'%s_rel.png'%name])
    pass
    


def iso_star(image_data,iterations):
    '''function that will isolate the bright objects in an image by iteratively
    masking out pixels above the median value of the remaining data.
    input:image data, number of iterations of masking loop
    output: binary image the same dimensions of the image data showing stars as value 1,
    binary image the same dimensions of the image data showing background as value 1.
    '''
    image_copy=copy.deepcopy(image_data)
    stars=np.zeros_like(image_data)
    background=np.zeros_like(image_data)
    nans=np.zeros_like(image_data)
    nans[~np.isfinite(image_data)]=1
    pool=image_copy[np.isfinite(image_copy)]
    med=np.median(pool)
    for i in range (iterations):
        pool=pool[np.where(pool>med)]
        med=np.median(pool)
    
    pool=pool[np.where(pool>med)]
    image_copy[np.where(nans==1)]=0
    stars[np.where(image_copy>med)]=1
    background[np.where(image_copy<0.5*med)]=1
    image_copy[np.where(nans==1)]=np.nan
    return stars,background



def background_estimate(image_data1,image_data2,use_for):
    '''function to estimate a value of the background of an image by splitting the image
    into small squares taking the average background value across all squares.
    function can be used in 2 ways:
    If you wish to use the background value as an offset to calabrate the data brightness
    (a.k.a. use_for='subtract') bright objects will be masked out and a mean average of the 
    remaining data is calculated which also considers negative values to give a base background value.
    If you wish to use the average value as a normalising constant to multiply the image by,
    (a.k.a. use_for='divide') only the brightest objects are masked out leaving a more representative data pool, 
    pixel values are also taken as their absolute value when calculating the mean.
    input: image_data1,image_data2,use of function
    output:binary images where 1 represents the background for images 1 & 2, the background value 
    for images 1 & 2, binary images where 1 represents masked out bright objects in images 1 & 2.
    '''

    y,x=image_data1.shape[0],image_data1.shape[1]
    #record nans to put back into image later 
    nan_mask=np.zeros_like(image_data1)
    nan_mask[~np.isfinite(image_data1)]=1
    #figure out how much of the image wont make it into a square
    remx=x%1500 
    remy=y%1500
    if 1500>x:
        x=1500
        remx=0
    if 1500>y:
        y=1500
        remy=0
    
    #split image into smaller squares for local background values
    #squares of 500 pixels
   
    background1=np.zeros_like(image_data1)
    stars1=np.zeros_like(image_data1)

    #remove this overlap from background calculation    
    background1[int(y-remy):,:]=np.nan
    background1[:,int(x-remx):]=np.nan
    
    for i in range (int(x/1500)):
        for j in range (int(y/1500)):
            square=copy.deepcopy(image_data1[1500*j:1500*j+1500,1500*i:1500*i+1500])               
            
            star_mask,bg_mask=iso_star(square,2)

            square[~np.isfinite(square)]=0 #avoiding "invalid value encountered in greater"
            if use_for=='subtract':         
                star_mask,bg_mask=iso_star(square,2) #isolated stars
                square[np.where(star_mask==1)]=np.nan #isolated background
#            if use_for =='divide':
#                star_mask,bg_mask=iso_star(square,1) #get rid of low values
#                square[np.where(star_mask==0)]=np.nan
                
            background1[1500*j:1500*j+1500,1500*i:1500*i+1500]=square
            stars1[1500*j:1500*j+1500,1500*i:1500*i+1500]=star_mask
    #re-introduce nans
    background1[np.where(nan_mask==1)]=np.nan
    plt.figure(),plt.imshow(background1)
    #mean background across whoel image
    if use_for=='divide':
        bg1=np.mean(np.absolute(background1[np.isfinite(background1)]))
    if use_for=='subtract':
        bg1=np.mean(background1[np.isfinite(background1)])
    

    #repeat for image 2 
    y,x=image_data2.shape[0],image_data2.shape[1]
    #record nans to put back into image later 
    nan_mask=np.zeros_like(image_data2)
    nan_mask[~np.isfinite(image_data2)]=1
    #figure out how much of the image wont make it into a square
    remx=x%1500 
    remy=y%1500
    if 1500>x:
        x=1500
        remx=0
    if 1500>y:
        y=1500
        remy=0

    background2=np.zeros_like(image_data2)
    stars2=np.zeros_like(image_data2) 
    
    background2[int(y-remy):,:]=np.nan
    background2[:,int(x-remx):]=np.nan       
    
    for i in range (int(x/1500)):
        for j in range (int(y/1500)):
            square=copy.deepcopy(image_data2[1500*j:1500*j+1500,1500*i:1500*i+1500])
            
            square[~np.isfinite(square)]=0 #avoiding "invalid value encountered in greater"
            if use_for=='subtract':         
                star_mask,bg_mask=iso_star(square,2) #isolated stars
                square[np.where(star_mask==1)]=np.nan #isolated background
#            if use_for =='divide':
#                star_mask,bg_mask=iso_star(square,1) #get rid of low values
#                square[np.where(star_mask==0)]=np.nan
                
                
            background2[1500*j:1500*j+1500,1500*i:1500*i+1500]=square
            stars2[1500*j:1500*j+1500,1500*i:1500*i+1500]=star_mask
                
    #re-introduce nans
    background2[np.where(nan_mask==1)]=np.nan
    
    if use_for=='divide':
        bg2=np.mean(np.absolute(background2[np.isfinite(background2)]))
    if use_for=='subtract':
        bg2=np.mean(background2[np.isfinite(background2)])
    
    return background1,background2,bg1,bg2,stars1,stars2

#function to get average ratio of brightness between images 
def bright_ratio(image_data1,image_data2):
    '''function to calculate the ratio of brightness between images by splitting
    the images into smaller squares and summing the total brightness of their pixels,
    the ratio of the sum in both images is calcuated and the mean value across all squares
    is taken.
    input:image_data1,image_data2
    output: average ratio of brightness
    '''
    RATIO = np.zeros_like(image_data1)

    x,y=RATIO.shape[1],RATIO.shape[0]
    
    nan_mask=np.zeros_like(image_data1)
    nan_mask[~np.isfinite(image_data1)]=1
    nan_mask[~np.isfinite(image_data2)]=1
    
    image_data1[np.where(nan_mask==1)]=0
    image_data2[np.where(nan_mask==1)]=0
    

    remx=x%1500
    remy=y%1500
    if 1500>x:
        x=1500
        remx=0
    if 1500>y:
        y=1500
        remy=0
    
    for i in range (int(x/1500)):
        for j in range (int(y/1500)):
            square1=copy.deepcopy(image_data1[1500*j:1500*j+1500,1500*i:1500*i+1500])
            
            square2=copy.deepcopy(image_data2[1500*j:1500*j+1500,1500*i:1500*i+1500])
   

            stars,bg=iso_star(square1,3)
            
            
            pool1=square1[np.where(stars==1)]
            pool2=square2[np.where(stars==1)]

            
            pool1=pool1[np.isfinite(pool1)]
            pool2=pool2[np.isfinite(pool2)]
        
            
            
            
            brightness1=sum(pool1)
            brightness2=sum(pool2)
            print(brightness1,brightness2)
            ratio=(brightness1/brightness2)
      
            RATIO[1500*j:1500*j+1500,1500*i:1500*i+1500]=ratio

            
    RATIO[int(y-remy):,:]=np.nan
    RATIO[:,int(x-remx):]=np.nan
    
    
    av_ratio=np.mean(RATIO[np.isfinite(RATIO)])
        
    return av_ratio
            
        

#function that crops & rotates images, and re-sizes pixels
#so that both images have same format
#save as fits images in folder 'fits_images'
def CPR(path_general):
    '''function that takes 2 images of the same object and crops them, changes the
    pixel sizes and roates the images so they are matching, making them easily compared.
    Once the images are edited they are saved in the fits_images folder
    input: path to general folder 
    '''
    print('CPR')
    #for each file in 'images' folder, load data & wcs, get coordinates of pixels
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
            
        else:
            print(filename)
            hdu1,wcs1,image_data1,hdu2,wcs2,image_data2=getimage(path_general,filename)        
            coordinates1,points1,min_ra1,max_ra1,min_dec1,max_dec1 = coords(image_data1,image_data2,wcs1,wcs2)
            coordinates2,points2,min_ra2,max_ra2,min_dec2,max_dec2 = coords(image_data2,image_data1,wcs2,wcs1)
            
            
            
            #rough estimate of which image is bigger (only repormat the larger image)
            xrange1 = max(coordinates1[1])-min(coordinates1[1])
            xrange2 = max(coordinates2[1])-min(coordinates2[1])
            yrange1 = max(coordinates1[0])-min(coordinates1[0])
            yrange2 = max(coordinates2[0])-min(coordinates2[0])
            
            if xrange1>=xrange2 : 
                if  yrange1>= yrange2:
                    new_image_data1, footprint = reproject_interp(hdu1, hdu2.header)
                    new_image_data2=image_data2 #reproject reformats image to same format
                    wcs=wcs2
                else:
                    if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                        new_image_data1, footprint = reproject_interp(hdu1, hdu2.header)
                        new_image_data2=image_data2 #reproject reformats image to same format
                        wcs=wcs2
                    else:
                        new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
                        new_image_data1=image_data1
                        wcs=wcs1
                        
            if xrange2>xrange1 :
                if  yrange2>=yrange1:
                    new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
                    new_image_data1=image_data1
                    wcs=wcs1
                
                else:
                    if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                        new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
                        new_image_data1=image_data1 #reproject reformats image to same format
                        wcs=wcs1
                    
                    else:
                        new_image_data1, footprint = reproject_interp(hdu1, hdu2.header)
                        new_image_data2=image_data2 #reproject reformats image to same format
                        wcs=wcs2
                        
                    
            
            if xrange2==xrange1:
                if  yrange2>yrange1:
                    new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
                    new_image_data1=image_data1
                    wcs=wcs1
                else:
                    if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                        new_image_data2, footprint = reproject_interp(hdu2, hdu1.header)
                        new_image_data1=image_data1 #reproject reformats image to same format
                        wcs=wcs1
                        

            
            
            if os.path.exists(path_general+'/fits_images/%s_1.fits'%filename)==True:
                os.remove(path_general+'/fits_images/%s_1.fits'%filename)
                os.remove(path_general+'/fits_images/%s_2.fits'%filename)
            header=wcs.to_header()
            hdu = fits.PrimaryHDU(new_image_data1,header=header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/fits_images/%s_1.fits'%filename)
            
            hdu = fits.PrimaryHDU(new_image_data2,header=header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/fits_images/%s_2.fits'%filename)
            
            
            
def bright_diff(path_general):
    '''function to calebrate the image brightness of both images to match eachother
    1)finds the background offset and adds it to image 2
    2)normalises both images by dividing by their new background value
    3)calculates ratio between image brightnesses and applies it to image2
    4)replaces the pre-existing files in 'fits_images' folder
    input:path to general folder
    '''
    
    print('bright_diff')
    for filename in sorted(os.listdir(path_general+'/images')):
        print(filename)
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            #get images
            hdu1,wcs1,new_image_data1,hdu2,wcs2,new_image_data2=get_fitsimage(path_general,filename)
            #blur
            new_image_data1=ndimage.filters.gaussian_filter(new_image_data1,3)
            new_image_data2=ndimage.filters.gaussian_filter(new_image_data2,3)
            #remember the pixels that are NaNs for later
            nan_mask=np.zeros_like(new_image_data1)
            nan_mask[~np.isfinite(new_image_data1)]=1
            nan_mask[~np.isfinite(new_image_data2)]=1
            
            
            #estimate background
            background1,background2,bg1,bg2,stars1,stars2=background_estimate(new_image_data1,new_image_data2,'subtract')
            
            #subtract background to centre stars around 0
            diff=bg1-bg2
#            new_image_data1=new_image_data1-bg1
            new_image_data2=new_image_data2+diff
            
            #get new background
            background1,background2,bg1,bg2,stars1,stars2=background_estimate(new_image_data1,new_image_data2,'divide')
            
            #normalise
            new_image_data1=new_image_data1/bg1
            new_image_data2=new_image_data2/bg2
            
            #calabrate scales
            ratio=bright_ratio(new_image_data1,new_image_data2)
            new_image_data2=new_image_data2*ratio
            


            #remove old fits files and replace with new processed fits files
            #give edited iamge wcs of other image
            os.remove(path_general+'/fits_images/%s_2.fits'%filename) #remove old file
            os.remove(path_general+'/fits_images/%s_1.fits'%filename)
            wcs=wcs2
            header=wcs.to_header()
            hdu = fits.PrimaryHDU(new_image_data2,header=header) #write new file
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/fits_images/%s_2.fits'%filename)
            
            header=wcs.to_header()
            hdu = fits.PrimaryHDU(new_image_data1,header=header) #write new file
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/fits_images/%s_1.fits'%filename)


def circles(path_general,keyword):    
    '''function that will create a cirlce around a given ra/dec coordinate
    input: path to general folder, target folder name
    output: binary images with cirlces as 1s to be used as contour data on a plot
    FUNCTION NOT FINISHED, ONLY CREATES CIRCLE AT CENTRE COORDINATE ATM
    INTENDED FOR IDENTIFYING VARIABLE STARS WITHIN THE IMAGE
    NEEDS RA/DEC ENTRIES FROM VARIABLE STAR ARCHIVE
    '''
    hdu1,wcs1,image_data1,hdu2,wcs2,image_data2=get_fitsimage(path_general,keyword)
    mask=np.zeros_like(image_data1)
    pixels=wcs1.wcs.crpix
#    pixels=wcs1.wcs_world2pix(ra,dec,1)
    px,py=int(pixels[1]),int(pixels[0])
    mask[px,py]=1
    points=np.where(mask==1)
    if len(points)==0:
        pass
    if len(points[0])==1:
        y,x=points[0],points[1]
        cy,cx=skimage.draw.circle(int(y),int(x),50)
        mask[cy,cx]=1
    if len(points[0])==2:
        for i in points: 
            y,x=i[0],i[1]
            cy,cx=skimage.draw.circle(int(y),int(x),50)
            mask[cy,cx]=1
    return mask
            

def sub(path_general):
    '''function to make subtraction and relative subtraction images 
    saves images in 'subtraction' & 'relative' folders respectively
    input:path to general folder
    '''
    print('sub')
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            hdu1,wcs1,new_image_data1,hdu2,wcs2,new_image_data2=get_fitsimage(path_general,filename)
            
#            #blur
            new_image_data1=ndimage.filters.gaussian_filter(new_image_data1,3)
            new_image_data2=ndimage.filters.gaussian_filter(new_image_data2,3)
            
            #record nans
            nan_mask=np.zeros_like(new_image_data1)
            nan_mask[~np.isfinite(new_image_data1)]=1
            nan_mask[~np.isfinite(new_image_data2)]=1
            
            #tempory remove nans
            new_image_data1[np.where(nan_mask==1)]=0
            new_image_data2[np.where(nan_mask==1)]=0
            
            #difference image
            subtraction=np.absolute(new_image_data1-new_image_data2)
            
            #relative image difference (percentage)
            relative=np.absolute(subtraction/new_image_data1)*100   
            
            #put the nans back
            subtraction[np.where(nan_mask==1)]=np.nan
            relative[np.where(nan_mask==1)]=np.nan
            
            #blur
            subtraction=ndimage.filters.gaussian_filter(subtraction,3)
            relative=ndimage.filters.gaussian_filter(relative,3)
            
            
            if os.path.exists(path_general+'/subtraction/%s_sub.fits'%filename)==True:
                os.remove(path_general+'/subtraction/%s_sub.fits'%filename)
            hdu = fits.PrimaryHDU(subtraction) #write new file
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/subtraction/%s_sub.fits'%filename)
            
            if os.path.exists(path_general+'/relative/%s_rel.fits'%filename)==True:
                os.remove(path_general+'/relative/%s_rel.fits'%filename)
            hdu = fits.PrimaryHDU(relative) #write new file
            hdul = fits.HDUList([hdu])
            hdul.writeto(path_general+'/relative/%s_rel.fits'%filename)
            

def out_save(path_general,subtraction,relative):
    '''function to calculate sensible figure brightness limmits and plot the images
    and save figures as .png images in the 'Zooniverse_upload' folder.
    input: path to general folder, to plot the regular images (1&2) subtraction=='no'
    and relative=='no',for subtraction=='yes' there imust be an existing subtraction image in the
    'subtraction' folder to plot the subtraction image, same for relative=='yes' to plot the relative image.
    note: both subtraction=='yes' AND relative=='yes' will not work.
    '''
    print('out_save')
    if os.path.exists(path_general+'/Zooniverse_upload/Manifest.csv')==True:
        os.remove(path_general+'/Zooniverse_upload/Manifest.csv')
    for filename in sorted(os.listdir(path_general+'/images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            print(filename)
            hdu1,wcs1,new_image_data1,hdu2,wcs2,new_image_data2=get_fitsimage(path_general,filename)
#            circle=circles(path_general,filename)
            
            #pool of real data (no NaNs)
            pool1=new_image_data1[np.isfinite(new_image_data1)]
            pool1=pool1[~np.isnan(pool1)]
            min1,max1=np.percentile(pool1,[90,99])
            
            
            pool2=new_image_data2[np.isfinite(new_image_data2)]
            pool2=pool2[~np.isnan(pool2)]
            min2,max2=np.percentile(pool2,[90,99])
            
            MAX=max(max1,max2)
            MIN=max(min1,min2)           
            if max1==MAX:
                MIN=min1
            else:
                MIN=min2
            
            
            
            if subtraction=='yes':
                path=path_general+'/subtraction/%s_sub.fits'%(filename)
                hdu_list=fits.open(path)
                hdu = fits.open(path)[0]
                subs=hdu.data
                
                hdu_list.close()
                pools=subs[np.isfinite(subs)]
                pools=pools[~np.isnan(pools)]
                mins,maxs=np.percentile(pools,[95,100])

                save_png('%s_sub.fits'%filename,mins,maxs,path_general,wcs1,'subtraction','no')
                os.remove(path_general+'/subtraction/%s_sub.fits'%filename)
                plt.close('all')
                
            if relative=='yes':  
                path=path_general+'/relative/%s_rel.fits'%(filename)
                hdu_list=fits.open(path)
                hdu = fits.open(path)[0]
                subs=hdu.data
                
                hdu_list.close()
                pools=subs[np.isfinite(subs)]
                pools=pools[~np.isnan(pools)]
                mins,maxs=np.percentile(pools,[95,100])
                save_png('%s_rel.fits'%filename,mins,maxs,path_general,wcs1,'relative','no')
                os.remove(path_general+'/relative/%s_rel.fits'%filename)
                plt.close('all')
            
            plt.close('all')
            save_png('%s_1.fits'%filename,MIN,MAX,path_general,wcs1,'normal','no')
#            os.remove(path_general+'/fits_images/%s_1.fits'%filename)
            plt.close('all')
            
            save_png('%s_2.fits'%filename,MIN,MAX,path_general,wcs2,'normal','no')
#            os.remove(path_general+'/fits_images/%s_2.fits'%filename)
            plt.close('all')

            manifest(filename,path_general)


def fold_check(path_general):
    '''function to check how many image target folders there are in the 
    'images' folder and run the image-processing functions if there are >30
    input: path to general folder
    '''
    path_folder=path_general+'/images'
    n=len(os.walk(path_folder).next()[1])
    if n>=30:
        print('images are ready for proccessing')
        remove_space(path_general,'images')
        remove_space(path_general,'LCO_images/raw')
        combine_fold(path_general)
        reject_dir(path_general)
        fz_remove(path_general)
        image_rename(path_general)
        CPR(path_general)
        bright_diff(path_general)
        sub(path_general)
        out_save(path_general,'yes','yes')
        files = os.listdir(path_folder)
        
        for f in files:
            if f=='.DS_Store': #problematic folder that shows up sometimes
                print('.')   
            else:            
                shutil.move(path_folder+'/'+f, path_general+'/finished_images')        
    else:
        print('insufficient images: min 30 objects')
    pass


    



#path_general='/Users/lewisprole/Documents/University/year3/summer_project'




