# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 16:09:13 2019

@author: annac
"""
import os
import shutil
from astropy.io import fits
from astropy.wcs.wcs import WCS
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from angles import arcs2d
import matplotlib.pyplot as plt


def remove_punctuation(path_general,dir1):
    '''function that will remove any spaces in the target folder name
    input: path to general folder, name of images folder e.g 'images'
    '''
    print('remove_punctuation')
    path=path_general+'/'+dir1
    for filename in sorted(os.listdir(path)):
        name_new=filename.replace(' ','_')
        name_new=name_new.replace('"','')
        name_new=name_new.replace("'",'')
        os.rename(path+'/%s'%(filename),path+'/%s'%(name_new))

        
def reject_dir(path_general):
    '''function that checks if there are at least 10 images per target
    if not the target folder is sent to reject folder
    input: path to general folder
    '''
    print('reject_dir')
    for filename in sorted(os.listdir(path_general+'/lightcurve_data')):
        n=len(sorted(os.listdir(path_general+'/lightcurve_data/%s'%filename)))
        #n2=len(sorted(os.listdir(path_general+'/archive_lightcurve/%s'%filename)))
        dest=path_general+'/rejects'
        #n=n1+n2
        if n<10:
            if os.path.exists(dest+'/'+filename)==True:
                shutil.rmtree(dest+'/'+filename,ignore_errors=True)
            shutil.move(path_general+'/lightcurve_data/%s'%filename,dest)
            #shutil.move(path_general+'/archive_lightcurve/%s'%filename,dest)

'''
def fz_remove(path_general):
    ''''''function which will uncompresses .fits.fz files into .fits files
    will check every target file in 'lightcurve_data' folder
    input:path to general folder
    ''''''
    print('fz_remove')
    for filename in sorted(os.listdir(path_general+'/lightcurve_data')):
        path=path_general+'/lightcurve_data/%s'%filename
        for name in os.listdir(path):
            if name.endswith('.fz')==True:
                image=fits.getdata(path+'/'+name)
                header=fits.getheader(path+'/'+name,ext=1)
                hdu=fits.PrimaryHDU(image,header=header)
                hdul=fits.HDUList([hdu])
                hdul.writeto(path+'/%s'%name[:-3])
                os.remove(path+'/%s'%name)
        pass
''' 

 
def image_rename(path_general):
    '''function that will go into each target file in 'lightcurve_data' folder
    and change the names of the images based on the target folder name
    e.g. target folder named 'M31',images renames 'M31_1.fits' & 'M31_2.fits'
    input: path to general folder
    '''
    print('image_rename')
    for filename in sorted(os.listdir(path_general+'/lightcurve_data')):
        path=path_general+'/lightcurve_data/%s'%filename
        n=len(os.listdir(path))
        names=os.listdir(path)

        for i in range(n):
            os.rename(path+'/%s'%(names[i]),path+'/%s%i.fits.fz'%(filename,i+1))
        for i in range (n):
            os.rename(path+'/%s%i.fits.fz'%(filename,i+1),path+'/%s_%i.fits.fz'%(filename,i+1))


        '''
        for i in range(len(os.listdir(path))):
            files=path+os.listdir(path)[i]
            file=fits.open(files[i])
            image=files[i].data
            header=fits.getheader(files[i],ext=0)
        '''            

            
def data_ext(path_general):
    '''function that will extract the data needed from each file and store in arrays
    will plot lightcurves for flux and magnitude
    input: path to general folder
    '''
    print(data_ext)
    for filename in sorted(os.listdir(path_general+'/lightcurve_data')):
        path=path_general+'/lightcurve_data/%s/'%filename
        MJD=np.zeros(len(os.listdir(path))) #modified julian date
        flux=np.zeros(len(os.listdir(path)))
        mag=np.zeros(len(os.listdir(path)))
        #plt.figure(figsize=(16,16)).add_subplot(1,1,1,projection=wcs)
        for i in range(len(os.listdir(path))):
            
            file=path+'%s_%s.fits.fz'%(filename,i+1) #each file in folder
            image1=fits.getdata(file,ext=1) #data
            header1=fits.getheader(file,ext=1) #header
            image2=fits.getdata(file,ext=2) #more data
            header2=fits.getheader(file,ext=2) #another header with flux and peak
            MJD[i]=header1['MJD-OBS']
            RA=header1['CAT-RA'] #RA of object
            Dec=header1['CAT-DEC'] #Dec of object
            c=SkyCoord(RA,Dec,unit=(u.hourangle,u.deg))
            decRA=c.ra.deg #decimal
            decDec=c.dec.deg #decimal
            wcs=WCS(header1).celestial
            
            #convert counts to flux/mag hopefully
            time=header1['EXPTIME']
            photflam=image2['FLUX'][34] #or flux? need photflam :(            
            totflux=image1/time #should be (data*photflam)/exptime
            totmag=-2.5*(np.log10(np.abs(totflux))) #+ photzpt zeropoint
            
            '''
            #find section of image that has object
            centreRA=header1['CRVAL1']
            centreDec=header1['CRVAL2']
            diffRA=np.abs(RA-centreRA)
            diffDec=np.abs(Dec-centreDec)
            res=header1['PIXSCALE'] #resolution arcsec/pixel
            res=arcs2d(res)
            a=diffRA/res
            b=diffDec/res
            
            centrex=header1['CRPIX1']
            centrey=header1['CRPIX2']
            x1=int(-(-(centrex-a)//1))
            x2=int(-(-(centrex+a)//1))
            y1=int(-(-(centrey-b)//1))
            y2=int(-(-(centrey+b)//1))
            section=image1[x1:x2,y1:y2]
            plt.figure()
            plt.imshow(section)
            plt.title('sss17a')
            '''
            
            #pixels of object
            
            x,y=wcs.all_world2pix(decRA,decDec,1)
            xpix=int(-(-x)//1)
            ypix=int(-(-y)//1)
            val=image1[xpix,ypix]
            
            #plt.subplot(4,4,i+1)
            #plt.imshow(image1[xpix-100:xpix+100,ypix-100:ypix+100])
            #plt.title('sss17a: {}'.format(i))
            #average flux and magnitude for object
            flux[i]=val/time
            #flux[i]=image2['FLUX'][0]
            mag[i]=-2.5*(np.log10(np.abs(flux[i])))
            xcoord=image2['X']
            ycoord=image2['Y']
            ra=image2['RA']
            dec=image2['DEC']
            x2,y2=wcs.all_world2pix(ra,dec,1)
            
            #plt.subplot(4,4,i+1)
            plt.figure().add_subplot(1,1,1,projection=wcs)
            plt.imshow(image1)#[xpix-100:xpix+100,ypix-100:ypix+100])
            plt.plot(x,y,'bx')
            plt.plot(xcoord,ycoord,'mx')
            plt.xlim(x-100,x+100)
            plt.ylim(y-100,y+100)
            #plt.plot(RA,Dec,'yx')
            #plt.plot(x2,y2,'gx')
            plt.title('sss17a: {}'.format(i))
            
            for i in range(len(xcoord)):
                plt.annotate(i,(xcoord[i],ycoord[i]))
            
            
        
        #plot lightcurve
        plt.figure()
        plt.plot(MJD,flux,'b.',label='Flux')
        #plt.plot(MJD,mag,'m.',label='Magnitude')
        plt.title('Lightcurve for sss17a binary neutron star merger')
        plt.xlabel('Modified Julian date')
        plt.ylabel('Value')
        plt.legend(loc='best')            
                
          
def lightcurve(path_general):
    #remove_punctuation(path_general,'lightcurve_data')
    #remove_punctuation(path_general,'archive_lightcurve')
    reject_dir(path_general)
    #fz_remove(path_general)
    image_rename(path_general)
    data_ext(path_general)
    
    