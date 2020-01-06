#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 09:56:42 2018

@author: annac
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import signal
from astropy.io import fits
import astropy.stats
import numpy.ma
import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D
#from skimage import img_as_bool, io, color, morphology
#from skimage.morphology import skeletonize
#from skimage import measure
#from numpy import random
#from skimage import measure
import copy
from astropy.wcs.wcs import WCS
#from astropy.utils.data import get_pkg_data_filename
#from astropy.coordinates import ICRS, Galactic
#from astropy import units as u
#import matplotlib.ticker as mticker
#import math
#from matplotlib.colors import LogNorm
#import matplotlib.colors as colors
from reproject import reproject_interp
import os
import csv
import shutil
#import time
#import threading
import skimage.draw
from xml.etree import ElementTree
import requests
#from scipy.fftpack import ifftn
#from scipy.signal import medfilt
#import pathlib
from scipy.fftpack import fft2, fft, fftshift
from scipy.stats import norm
import scipy as sp

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

path_general='/Users/annac/Documents/Summer Project'
remove_punctuation(path_general,'images')
remove_punctuation(path_general,'LCO_images/raw')
combine_fold(path_general)
reject_dir(path_general)
fz_remove(path_general)
image_rename(path_general)
CPR(path_general)
bright_diff(path_general)
sub(path_general)
out_save(path_general,'yes','yes')


--------------------------------------------------------------------
'''

def remove_punctuation(path_general,dir1):
    '''function that will remove any spaces in the target folder name
    input: path to general folder, name of images folder e.g 'images'
    '''
    print('remove_punctuation')
    path=path_general+'/'+dir1
    for filename in sorted(os.listdir(path)):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            name_new=filename.replace(' ','_')
            name_new=name_new.replace('"','')
            name_new=name_new.replace("'",'')
            os.rename(path+'/%s'%(filename),path+'/%s'%(name_new))



def combine_fold(path_general):
    '''function to move target folders from the LCO downloads folder to the main
    images folder
    input: path to general folder
    '''
    print('combine_fold')
    files=os.listdir(path_general+'/LCO_images/raw')
    dest=path_general+'/images'
    for f in files:
        if f=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            if os.path.exists(dest+'/'+f)==True:
                shutil.rmtree(dest+'/'+f,ignore_errors=True)
                #shutil.move(dest+'/'+f,path_general+'/rejects/'+f)
            shutil.move(path_general+'/LCO_images/raw/'+f,dest)


def reject_dir(path_general):
    '''function that checks if there are at least 2 images per target
    if not the target folder is sent to reject folder
    input: path to general folder
    '''
    print('reject_dir')
    for filename in sorted(os.listdir(path_general+'/images')):
        #for file in sorted(os.listdir(path_general+'/archive_images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            n1=len(sorted(os.listdir(path_general+'/images/%s'%filename)))
            n2=len(sorted(os.listdir(path_general+'/archive_images/%s'%filename)))
            dest=path_general+'/rejects'
            n=n1+n2
            if n<2:
                if os.path.exists(dest+'/'+filename)==True:
                    shutil.rmtree(dest+'/'+filename,ignore_errors=True)
                shutil.move(path_general+'/images/%s'%filename,dest)
                shutil.move(path_general+'/archive_images/%s'%filename,dest)


def radial_mean(data,coords):
    '''function to find radial mean
    input: image data and coordinates of centre of radius
    '''
    y,x=np.indices((data.shape))
    r=np.sqrt((x-coords[0])**2+(y-coords[1])**2)
    rbin=(30*r/r.max()).astype(np.int)
    radialmean=ndimage.mean(data,labels=rbin,index=np.arange(1,rbin.max()+1))
            
    return radialmean
        
        
def resolution(path_general,obj):
    '''function to match up resolutions between LCO and archive images
    input: path to general folder
    '''
    print('resolution')
    for f in os.listdir(path_general+'/archive_images/%s'%obj):
        files=path_general+'/archive_images/%s/%s'%(obj,f)
        intdata=fits.getdata(files,ext=0)
        header=fits.getheader(files,ext=0)
        intdata=np.require(intdata,dtype=np.float32)
        f1=fft2(intdata)
        f1abs=np.log10(np.abs(f1))
        plt.figure()
        plt.imshow(f1abs)
        plt.figure()
        plt.plot(f1abs[0,0:150])
        plt.xscale('log')
        plt.yscale('log')
        plt.figure()
        plt.imshow(intdata)
        
        data=intdata[0:150,0:150]
        coords=(0,0)
        rad=radial_mean(data,coords)

        plt.figure()
        plt.plot(rad)
        plt.title('%s archive image radial mean'%f[:-5])
        
    file=path_general+'/images/%s/%s_3.fits'%(obj,obj)
    intdata=fits.getdata(file,ext=0)
    header=fits.getheader(file,ext=0)
    intdata=np.require(intdata,dtype=np.float32)
    f1=fft2(intdata)
    f1abs=np.log10(np.abs(f1))
    plt.figure()
    plt.imshow(f1abs)
    plt.figure()
    plt.plot(f1abs[0,0:150])
    plt.xscale('log')
    plt.yscale('log')
    plt.figure()
    plt.imshow(intdata)
        
    data=intdata[0:150,0:150]
    coords=(0,0)
    rad=radial_mean(data,coords)

    plt.figure()
    plt.plot(rad)
    plt.title('%s LCO image radial mean'%obj)
        
    '''
    #mu=np.ndarray.mean(intdata)
    #sig=np.ndarray.std(intdata)
    data=np.reshape(intdata,((np.shape(intdata)[0])*(np.shape(intdata)[1])))

    #for i in range(len(intdata)):
    #   for j in range(len(intdata[i])):
    # Fit a normal distribution to the data:
    mu,std=norm.fit(data)
        
    # Plot the histogram.
    plt.figure()
    plt.hist(data,bins=200)
        
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)
        
    plt.show()
    '''

def getimage(path_general,keyword):                        
    '''function that grabs fits image infomation from .fits file    
    located in the 'images' folder                                  
    input: path to general folder, name of image target folder      
    returns:hdu,wcs & image data of image                           
    '''                                                             
                                                                    
    path=path_general+'/images/%s/%s_3.fits'%(keyword,keyword)            
    #hdu_list1=fits.open(path)                                      
    hdu1=fits.open(path)[0]                                       
    wcs1=WCS(hdu1.header)                                         
    image_data1=hdu1.data                                           
    #hdu_list1.close()                                              
                                                                    
    return hdu1,wcs1,image_data1                                    
                                                                    
                                                                    
#similar function used later for extracting new fits images         
def get_fitsimage(path_general,keyword):                   
    '''same function as getimage but taking images from             
    the 'fits_images' folder instead of 'images'                    
    input: path to general folder, name of image target folder      
    returns:hdu,wcs & image data of image                           
    '''                                                             
                                                                    
    path=path_general+'/fits_images/%s/%s_1.fits'%(keyword,keyword)       
    hdu1=fits.open(path)[0]                                       
    wcs1=WCS(hdu1.header)                                         
    image_data1=hdu1.data                                           
                                                                    
    return hdu1,wcs1,image_data1                                    


#similar function for the archive images
def get_archiveimage(path_general,keyword,survey):
    '''same function as getimage but taking images from
    the 'archive_image' folder instead of 'images'
    input: path to general folder, name of image target folder
    returns:hdu,wcs & image data of image
    '''
    for i in survey:
        path=path_general+'/archive_images/%s/%s_%s.fits'%(keyword,keyword,i)
        #hdu_list2=fits.open(path)
        hdu2=fits.open(path)[0]
        wcs2=WCS(hdu2.header)
        image_data2=hdu2.data
        #hdu_list2.close()
            
    return hdu2,wcs2,image_data2



def get_archivefits(path_general,keyword,filename,survey):
    '''same function as getimage but taking images from
    the 'archive_fits' folder instead of 'images'
    input: path to general folder, name of image target folder
    returns:hdu,wcs & image data of image
    '''
    for i in survey:
        path=path_general+'/archive_fits/%s/%s_%s.fits'%(keyword,keyword,i)
        hdu2=fits.open(path)[0]
        wcs2=WCS(hdu2.header)
        image_data2=hdu2.data
        
    return hdu2,wcs2,image_data2



def coords(data1,data2,system1,system2):
    '''function to grab galactic coordinates of image pixels
    input:image data from image1, same for image 2, wcs of image1, wcs of image2
    output:list of (ra,dec) coords, minimum ra, maximum ra, minimum dec, maximum dec
    '''
    x=data1.shape[1]
    y=data1.shape[0]
    grid=np.ones((y,x))
    points=np.where(grid==1)
    pointsy,pointsx=points[0],points[1]
    coords1=system1.wcs_pix2world(pointsx,pointsy,1)
    coordsy,coordsx=coords1[1],coords1[0]
    c_join=(coordsy,coordsx)
    min_ra,max_ra=min(coordsy),max(coordsy)
    min_dec,max_dec=min(coordsx),max(coordsx)

    return c_join,points,min_ra,max_ra,min_dec,max_dec


def fz_remove(path_general): #may not need?
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
                    image=fits.getdata(path+'/'+name)
                    header=fits.getheader(path+'/'+name,ext=1)
                    hdu=fits.PrimaryHDU(image,header=header)
                    hdul=fits.HDUList([hdu])
                    hdul.writeto(path+'/%s'%name[:-3])
                    os.remove(path+'/%s'%name)
        pass



def image_rename(path_general):
    '''function that will go into each target file in 'images' folder
    and change the names of the images based on the target folder name
    e.g. target folder named 'M31',images renames 'M31_1.fits' & 'M31_2.fits'
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
            for i in range(n):
                os.rename(path+'/%s'%(names[i]),path+'/%s%i.fits'%(filename,i+1))
            for i in range (n):
                os.rename(path+'/%s%i.fits'%(filename,i+1),path+'/%s_%i.fits'%(filename,i+1))    



def save_png(name,obj,Vmin,Vmax,path_general,wcs,types,contour):
    '''function that will save the .fits images as .png images
    input: target name, min and max displayed values, path to general folder,
    wcs of image, subtraction image of target, relative subtraction image of target,
    title of figure, any contour data
    '''
    n=name.find('.')
    t_name=name[:n]
    if types=='subtraction':
        fits_data=fits.getdata(path_general+'/subtraction/%s/%s'%(obj,name))
        c='viridis'
        title='%s: subtraction image'%t_name[:-3]
        #norm=LogNorm()
    if types=='relative':
        fits_data=fits.getdata(path_general+'/relative/%s/%s'%(obj,name))
        #c='plasma'
        c='viridis'
        title='%s: percentage change'%t_name[:-3]
        #norm=LogNorm()
    if types=='normal':
        fits_data=fits.getdata(path_general+'/fits_images/%s/%s'%(obj,name))
        #c='magma'
        c='viridis'
        #n=name.find('.')
        #t_name=name[:n]
        title='%s'%t_name
        #norm=colors.SymLogNorm(Vmax/2)
    if types=='archive':
        fits_data=fits.getdata(path_general+'/archive_fits/%s/%s'%(obj,name))
        #c='magma'
        c='viridis'
        #n=name.find('.')
        #t_name=name[:n]
        title='%s'%t_name
        #norm=colors.SymLogNorm(Vmax/2)

    plt.subplot(projection=wcs)
    plt.imshow(fits_data,vmin=Vmin,vmax=Vmax,cmap=c)#,norm=norm)
    plt.grid(color='white',ls='solid')
    plt.colorbar()
    #if Vmax-Vmin>=1:
        #plt.colorbar(label="log scale")
    #else:
        #plt.colorbar(ticks=[Vmin,Vmax],label="log scale")
    plt.title(title)
    imagename=name.replace('.fits','.png')
    plt.savefig(path_general+'/Zooniverse_upload/%s'%imagename)
    plt.close('all')



def manifest(name,path_general,survey): #need for archive this isnt right
    '''function to create manifest to be uploaded to zooniverse
    .csv telling zooniverse which images belong to the same subject set,
    saves it in the Zoonivere_upload folder
    input:target name, path to general folder
    '''
    if os.path.exists(path_general+"/Zooniverse_upload/Manifest.csv")==False:
        with open(path_general+"/Zooniverse_upload/Manifest.csv", "w+") as csvfile:
            write=csv.writer(csvfile)
            for i in survey:
                write.writerow(['image1','image2','subtraction','relative'])
                write.writerow(['%s.png'%name,'%s_%s.png'%(name,i),'%s_sub.png'%name,'%s_rel.png'%name])

    else:
        with open(path_general+"/Zooniverse_upload/Manifest.csv", "a") as csvfile:
            write=csv.writer(csvfile)
            for i in survey:
                write.writerow(['%s.png'%name,'%s_%s.png'%(name,i),'%s_sub.png'%name,'%s_rel.png'%name])
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
    for i in range(iterations):
        pool=pool[np.where(pool>med)]
        med=np.median(pool)

    image_copy[np.where(nans==1)]=0
    stars[np.where(image_copy>med)]=1
    background[np.where(image_copy<0.5*med)]=1
    image_copy[np.where(nans==1)]=np.nan
    return stars,background



def med_filt(diff_data,square_size):
    '''function that applies a median filter to an image
    median pixel value of a large square of dimension 'square_size' is calculated
    and replaces all pixel values in the square, the square moves by (square_size/5)
    each turn, overlapping 5 times.
    input:image,square size.
    output:median filter of image.
    '''

    y,x=diff_data.shape[0],diff_data.shape[1]
    mfilter=np.zeros_like(diff_data)

    s=int(square_size)
    r=int(s/5)
    for i in range (int(x/r-s/r+1)):
        for j in range (int(y/r-s/r+1)):
            square=diff_data[j*r:j*r+s,i*r:i*r+s]

            pool=square[np.isfinite(square)]

            med=np.median(pool)


            mfilter[int(j*r):int(j*r+s),int(i*r):int(i*r+s)]=med

    #number of pixels left at the the end that didn't make it into a square
    xl=x%s
    yl=y%s
    #make them nans
    mfilter[int(-yl):,:]=np.nan
    mfilter[:,int(-xl):]=np.nan

    return mfilter

def offset_filter(image_data1,image_data2):
    '''function to create a filter of the varying offset between 2 images
    to be added to image2, by applying a median filter to the difference image.
    The stars are removed to consider background values only.
    input:image1,image2
    output:image2 after the offset filter has been added.
    '''
    #get star positions
    stars1,b=iso_star(image_data1,2)
    stars2,b=iso_star(image_data2,2)
    #shared star map for both images
    stars=np.zeros_like(image_data1)
    stars[np.where(stars1==1)]=1
    stars[np.where(stars2==1)]=1
    #map image differences
    diff_filter=image_data1-image_data2     
    #remove star differences
    diff_filter[np.where(stars==1)]=np.nan
    #blur with median filter
    diff_filter=med_filt(diff_filter,500)
    diff_filter=ndimage.filters.gaussian_filter(diff_filter,10)
    #apply to new_image_data2
    image_data2=image_data2+diff_filter
    return image_data2


def ratio_filter(image_data1,image_data2):
    '''function to create a filter of the varying ratio between 2 images
    to be multiplied to image2, by applying a median filter to the ratio image.
    input:image1,image2
    output:image2 after the offset filter has been multiplied.
    '''
    stars1,b=iso_star(image_data1,2)
    stars2,b=iso_star(image_data2,2)
    #shared star map for both images
    stars=np.zeros_like(image_data1)
    stars[np.where(stars1==1)]=1
    stars[np.where(stars2==1)]=1
    #map image differences
    diff_filter=image_data1/image_data2
    diff_filter[np.where(stars==0)]=np.nan
    #blur with median filter

    diff_filter=med_filt(diff_filter,500)
    diff_filter=ndimage.filters.gaussian_filter(diff_filter,10)
    image_data2=image_data2*diff_filter
    return image_data2


def background_estimate(image_data1,size):
    '''function to create a background fiter by removing the stars and
    generating nosie based on the surrounding background to fill in the stars.
    Noise is generated in small overlapping boxes to represent local background.
    input: image, size of the sample box for noise generation.
    '''

    y,x=image_data1.shape[0],image_data1.shape[1]
    background1=np.zeros_like(image_data1)

    #record nans to put back into image later
    nan_mask=np.zeros_like(image_data1)
    nan_mask[~np.isfinite(image_data1)]=1

    #remove stars
    stars1,b=iso_star(image_data1,1)
    image=copy.deepcopy(image_data1)
    image[np.where(stars1==1)]=np.nan
    #figure out how much of the image wont make it into a square

    #square overlap
    s=size
    r=s/5
    X=int(x/r-s/r+1)
    Y=int(y/r-s/r+1)
    #figuring out how much of the image doesn't get covered at the end
    remx=x%r
    remy=y%r
    #accounting for if the square is bigger than the image
    if X<=0:
        X=1
        remx=0
    if Y<=0:
        Y=1
        remy=0

    #split image into smaller squares
    for i in range(int(X)):
        for j in range(int(Y)):

            square=copy.deepcopy(image[int(r*j):int(r*j+s),int(r*i):int(r*i+s)])
            ys,xs=square.shape[0],square.shape[1]
            stars=stars1[int(r*j):int(r*j+s),int(r*i):int(r*i+s)]
            #remove nans
            pool=square[np.isfinite(square)]

            #create noise baed on surrounding values
            noise=numpy.random.normal(np.median(pool),np.max(pool),(ys,xs))
            #fill in the stars with background nosie
            square[np.where(stars==1)]=noise[np.where(stars==1)]


            background1[int(r*j+2*r):int(r*j+3*r),int(r*i+2*r):int(r*i+3*r)]=square[int(2*r):int(3*r),int(2*r):int(3*r)]


    background1[:int(2*r),:]=np.nan
    background1[:,:int(2*r)]=np.nan

    background1[:,int(-remx-2*r):]=np.nan
    background1[int(-remy-2*r):,:]=np.nan


    return background1



def get_variables(fromRA,toRA,fromDec,toDec):
    '''function to find any variable stars that might be in the image
    input: minimum and maximum ra, minimum and maximum dec
    output: list of ra and decs
    '''

    payload={'fromra':fromRA,
           'tora':toRA,
           'fromdec':fromDec,
           'todec':toDec}

    response=requests.get('https://www.aavso.org/vsx/index.php?view=api.list',params=payload)
    tree=ElementTree.fromstring(response.content)
    ras=[]
    decs=[]
    for variable in tree.findall('VSXObject'):
        ra=variable.find('RA2000').text
        dec=variable.find('Declination2000').text
        ras.append(ra)
        decs.append(dec)
    return ras,decs


def circles(path_general,keyword,survey):
    '''function that will create a circle around a given ra/dec coordinate
    input: path to general folder, target folder name
    output: binary images with circles as 1s to be used as contour data on a plot
    '''
    for i in os.listdir(path_general+'/archive_images/%s'%keyword):
        hdu1,wcs1,image_data1=get_fitsimage(path_general,keyword)
        hdu2,wcs2,image_data2=get_archivefits(path_general,keyword,i,survey)
        mask=np.zeros_like(image_data1)
    
        c_join,points,min_ra,max_ra,min_dec,max_dec=coords(image_data1,image_data2,wcs1,wcs2)
    
        ras,decs=get_variables(min_ra,max_ra,min_dec,max_dec)
    
    
        if len(ras)>0:
            ras=np.array(ras,dtype='float64')
            decs=np.array(decs,dtype='float64')
            pixels=wcs1.wcs_world2pix(ras,decs,1)
            px,py=np.int_(pixels[1]),np.int_(pixels[0])
            mask[px,py]=1
            points=np.where(mask==1)
            if len(points)==0:
                pass
            if len(points[0])==1:
                y,x=points[0],points[1]
                cy,cx=skimage.draw.circle(int(y),int(x),50)
                mask[cy,cx]=1
            if len(points[0])>1:
                for i in points:
                    y,x=i[0],i[1]
                    cy,cx=skimage.draw.circle(int(y),int(x),50)
                    mask[cy,cx]=1
    return mask



#function that crops & rotates images, and re-sizes pixels
#so that both images have same format
#save as fits images in folder 'fits_images'
def CPR(path_general,survey):
    '''function that takes 2 images of the same object and crops them, changes the
    pixel sizes and rotates the images so they are matching, making them easily compared.
    Once the images are edited they are saved in the fits_images folder
    input: path to general folder
    '''
    print('CPR')
    #for each file in 'images' folder, load data & wcs, get coordinates of pixels
    for filename in sorted(os.listdir(path_general+'/images')):
        #for file in sorted(os.listdir(path_general+'/archive_images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')

        else:
            print(filename)
            for i in os.listdir(path_general+'/archive_images/%s'%filename):
                hdu1,wcs1,image_data1=getimage(path_general,filename)
                hdu2,wcs2,image_data2=get_archiveimage(path_general,filename,survey)
                coordinates1,points1,min_ra1,max_ra1,min_dec1,max_dec1=coords(image_data1,image_data2,wcs1,wcs2)
                coordinates2,points2,min_ra2,max_ra2,min_dec2,max_dec2=coords(image_data2,image_data1,wcs2,wcs1)
        
        
                #rough estimate of which image is bigger (only reformat the larger image)
                xrange1=max(coordinates1[1])-min(coordinates1[1])
                xrange2=max(coordinates2[1])-min(coordinates2[1])
                yrange1=max(coordinates1[0])-min(coordinates1[0])
                yrange2=max(coordinates2[0])-min(coordinates2[0])

                if xrange1>=xrange2:
                    if  yrange1>=yrange2:
                        new_image_data1,footprint=reproject_interp(hdu1,hdu2.header)
                        new_image_data2=image_data2 #reproject reformats image to same format
                        wcs=wcs2
                    else:
                        if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                            new_image_data1,footprint=reproject_interp(hdu1,hdu2.header)
                            new_image_data2=image_data2 #reproject reformats image to same format
                            wcs=wcs2
                        else:
                            new_image_data2,footprint=reproject_interp(hdu2,hdu1.header)
                            new_image_data1=image_data1
                            wcs=wcs1
        
                if xrange2>xrange1 :
                    if  yrange2>=yrange1:
                        new_image_data2,footprint=reproject_interp(hdu2,hdu1.header)
                        new_image_data1=image_data1
                        wcs=wcs1
        
                    else:
                        if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                            new_image_data2,footprint=reproject_interp(hdu2,hdu1.header)
                            new_image_data1=image_data1 #reproject reformats image to same format
                            wcs=wcs1
        
                        else:
                            new_image_data1,footprint=reproject_interp(hdu1,hdu2.header)
                            new_image_data2=image_data2 #reproject reformats image to same format
                            wcs=wcs2
        
        
        
                if xrange2==xrange1:
                    if  yrange2>yrange1:
                        new_image_data2,footprint=reproject_interp(hdu2,hdu1.header)
                        new_image_data1=image_data1
                        wcs=wcs1
                    else:
                        if np.abs(xrange1-xrange2)>np.abs(yrange1-yrange2):
                            new_image_data2,footprint=reproject_interp(hdu2,hdu1.header)
                            new_image_data1=image_data1 #reproject reformats image to same format
                            wcs=wcs1


                #create object folder inside data folder if not existent
                if not os.path.exists(path_general+'/fits_images/%s'%filename):
                    os.mkdir(path_general+'/fits_images/%s'%filename)
                #moving to fits_images
                if os.path.exists(path_general+'/fits_images/%s/%s_1.fits'%(filename,filename))==True:
                    os.remove(path_general+'/fits_images/%s/%s_1.fits'%(filename,filename))
                header=wcs.to_header()
                hdu=fits.PrimaryHDU(new_image_data1,header=header)
                hdul=fits.HDUList([hdu])
                hdul.writeto(path_general+'/fits_images/%s/%s_1.fits'%(filename,filename),overwrite=True)
                hdul.close()
        
                #create object folder inside data folder if not existent
                if not os.path.exists(path_general+'/archive_fits/%s'%filename):
                    os.mkdir(path_general+'/archive_fits/%s'%filename)
                #moving to archive_fits
                if os.path.exists(path_general+'/archive_fits/%s/%s'%(filename,i))==True:
                    os.remove(path_general+'/archive_fits/%s/%s'%(filename,i))
                header=wcs.to_header()
                hdu=fits.PrimaryHDU(new_image_data2,header=header)
                hdul=fits.HDUList([hdu])
                hdul.writeto(path_general+'/archive_fits/%s/%s'%(filename,i),overwrite=True)
                hdul.close()
            


def bright_diff(path_general,survey):
    '''function to calebrate the image brightness of both images to match eachother
    1)finds the background offset and adds it to image 2
    2)normalises both images by dividing by their new background value
    3)calculates ratio between image brightnesses and applies it to image2
    4)replaces the pre-existing files in 'fits_images' folder
    input:path to general folder
    '''

    print('bright_diff')
    for filename in sorted(os.listdir(path_general+'/fits_images')):
        #for file in sorted(os.listdir(path_general+'/archive_fits')):
        print(filename)
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            #get images
            for i in os.listdir(path_general+'/archive_fits/%s'%filename):
                hdu1,wcs1,new_image_data1=get_fitsimage(path_general,filename)
                hdu2,wcs2,new_image_data2=get_archivefits(path_general,filename,i,survey)
    
    
    
                #normalise
                new_image_data1=ndimage.filters.gaussian_filter(new_image_data1,3)
                d1=new_image_data1[np.isfinite(new_image_data1)]
                new_image_data1=(new_image_data1-np.mean(d1))/np.std(d1)
                    
                new_image_data2=ndimage.filters.gaussian_filter(new_image_data2,3)
                d2=new_image_data2[np.isfinite(new_image_data2)]
                new_image_data2=(new_image_data2-np.mean(d2))/np.std(d2)
                    
                #apply offset filter
                new_image_data2=offset_filter(new_image_data1,new_image_data2) #makes data2 into nan
                    
                    
                '''
                #remove old fits files and replace with new processed fits files
                #give edited iamge wcs of other image
                with Image.open(path_general+'/fits_images/%s_2.fits'%filename):
                    os.remove(path_general+'/fits_images/%s_2.fits'%filename) #remove old file
                with Image.open(path_general+'/fits_images/%s_1.fits'%filename):
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
                '''




def sub(path_general,survey):
    '''function to make subtraction and relative subtraction images
    saves images in 'subtraction' & 'relative' folders respectively
    input:path to general folder
    '''
    print('sub')
    for filename in sorted(os.listdir(path_general+'/images')):
        #for file in sorted(os.listdir(path_general+'/archive_images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            for i in os.listdir(path_general+'/archive_images/%s'%filename):
                hdu1,wcs1,new_image_data1=get_fitsimage(path_general,filename)
                hdu2,wcs2,new_image_data2=get_archivefits(path_general,filename,i,survey)
    
    
                #record nans
                nan_mask=np.zeros_like(new_image_data1)
                nan_mask[~np.isfinite(new_image_data1)]=1
                nan_mask[~np.isfinite(new_image_data2)]=1
    
    
                #difference image
                subtraction=np.absolute(new_image_data1-new_image_data2)
                    
                #relative image difference (percentage)
                relative=np.absolute(subtraction/new_image_data1)*100
    
                #blur
                relative=ndimage.filters.gaussian_filter(relative,3)
                subtraction=ndimage.filters.gaussian_filter(subtraction,3)
    
    
                #create object folder inside data folder if not existent
                if not os.path.exists(path_general+'/subtraction/%s'%filename):
                    os.mkdir(path_general+'/subtraction/%s'%filename)
                if os.path.exists(path_general+'/subtraction/%s/%s_sub.fits'%(filename,i[:-5]))==True:
                    os.remove(path_general+'/subtraction/%s/%s_sub.fits'%(filename,i[:-5]))
                hdu=fits.PrimaryHDU(subtraction) #write new file
                hdul=fits.HDUList([hdu])
                hdul.writeto(path_general+'/subtraction/%s/%s_sub.fits'%(filename,i[:-5]),overwrite=True)
                #hdul.close()
    
                #create object folder inside data folder if not existent
                if not os.path.exists(path_general+'/relative/%s'%filename):
                    os.mkdir(path_general+'/relative/%s'%filename)
                if os.path.exists(path_general+'/relative/%s/%s_rel.fits'%(filename,i[:-5]))==True:
                    os.remove(path_general+'/relative/%s/%s_rel.fits'%(filename,i[:-5]))
                hdu=fits.PrimaryHDU(relative) #write new file
                hdul=fits.HDUList([hdu])
                hdul.writeto(path_general+'/relative/%s/%s_rel.fits'%(filename,i[:-5]),overwrite=True)
                #hdul.close()


def out_save(path_general,subtraction,relative,survey):
    '''function to calculate sensible figure brightness limits and plot the images
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
        #for file in sorted(os.listdir(path_general+'/archive_images')):
        if filename=='.DS_Store': #problematic folder that shows up sometimes
            print('.')
        else:
            print(filename)
            for i in os.listdir(path_general+'/archive_images/%s'%filename):
                hdu1,wcs1,new_image_data1=get_fitsimage(path_general,filename)
                hdu2,wcs2,new_image_data2=get_archivefits(path_general,filename,i,survey)
    
                circle=circles(path_general,filename,survey)
    
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
                    path=path_general+'/subtraction/%s/%s_sub.fits'%(filename,i[:-5])
                    hdu_list=fits.open(path)
                    hdu=fits.open(path)[0]
                    subs=hdu.data
    
                    hdu_list.close()
                    pools=subs[np.isfinite(subs)]
                    pools=pools[~np.isnan(pools)]
                    mins,maxs=np.percentile(pools,[95,100])
    
                    save_png('%s_sub.fits'%i[:-5],'%s'%filename,mins,maxs,path_general,wcs1,'subtraction',circle)
                    #os.remove(path_general+'/subtraction/%s_sub.fits'%filename)
    
                if relative=='yes':
                    path=path_general+'/relative/%s/%s_rel.fits'%(filename,i[:-5])
                    hdu_list=fits.open(path)
                    hdu=fits.open(path)[0]
                    subs=hdu.data
                        
                    hdu_list.close()
                    pools=subs[np.isfinite(subs)]
                    pools=pools[~np.isnan(pools)]
                    mins,maxs=np.percentile(pools,[95,99])
                    
                    save_png('%s_rel.fits'%i[:-5],'%s'%filename,mins,maxs,path_general,wcs1,'relative','no')
                    #os.remove(path_general+'/relative/%s_rel.fits'%filename)
    
            save_png('%s_1.fits'%filename,'%s'%filename,MIN,MAX,path_general,wcs1,'normal','no')
            #os.remove(path_general+'/fits_images/%s_1.fits'%filename)
    
            for i in survey:
                save_png('%s_%s.fits'%(filename,i),'%s'%filename,MIN,MAX,path_general,wcs2,'archive','no')
                #os.remove(path_general+'/fits_images/%s_2.fits'%filename)
    
            manifest(filename,path_general,survey)
                

def fold_check(path_general,survey,obj):    
    '''function to check how many image target folders there are in the
    'images' folder and run the image-processing functions if there are >1.
    function will move original fits images to the 'finished_images' folder after completion.
    input: path to general folder
    '''
    path_folder=path_general+'/images'
    path_folder_LCO=path_general+'/LCO_images' #???????????
    path_folder_archives=path_general+'/archive_images'
    n1=len(os.walk(path_folder).__next__()[1])
    n2=len(os.walk(path_folder_LCO).__next__()[1])
    n3=len(os.walk(path_folder_archives).__next__()[1])
    n=n1+n2+n3
    '''
    print('images are ready for processing')
    remove_punctuation(path_general,'images')
    remove_punctuation(path_general,'LCO_images/raw')
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
    '''
    
    if n>=1:
        print('images are ready for processing')
        remove_punctuation(path_general,'images')
        remove_punctuation(path_general,'LCO_images/raw')
        remove_punctuation(path_general,'archive_images')
        combine_fold(path_general)
        reject_dir(path_general)
        fz_remove(path_general)
        image_rename(path_general)
        resolution(path_general,obj)
        CPR(path_general,survey)
        bright_diff(path_general,survey)
        sub(path_general,survey)
        out_save(path_general,'yes','yes',survey)
        files=os.listdir(path_folder)

        for f in files:
            if f=='.DS_Store': #problematic folder that shows up sometimes
                print('.')
            else:
                if os.path.exists(path_general+'/finished_images/%s'%f)==True:
                    shutil.rmtree(path_general+'/finished_images/%s'%f,ignore_errors=True)
                shutil.move(path_folder+'/'+f,path_general+'/finished_images')
    else:
        print('insufficient images: min 1 objects')
    pass
    







