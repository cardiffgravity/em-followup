# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 16:27:13 2019

@author: annac
"""

from astroquery.sdss import SDSS
from astroquery.skyview import SkyView
import astropy.coordinates as coords
from astropy import units as u
from astropy.io import fits
from astropy.wcs.wcs import WCS
import os
import shutil
#import datetime

'''
how to perform an object cross-ID with SDSS
start with the position of a source found in another survey (LCO archive)
search within a 5 arcsecond radius for optical counterparts in SDSS
keyword argument spectro requires matches to have spectroscopy not just photometry

SkyView offers a cutout service for a number of imaging surveys
to see the list of surveys, use the list_surveys method
'''

def archives(path_general,obj):#,verbose=False):
    path=path_general+'/lightcurve_data/%s/'%obj
    file=path+os.listdir(path)[0]
    
    image=fits.getdata(file)
    header=fits.getheader(file,ext=1)
    hdu=fits.PrimaryHDU(data=image,header=header)
    hdul=fits.HDUList([hdu])
    '''
    hdul.writeto(file[:-3])
    d=int(datetime.datetime.now().timestamp())
    hdul.writeto(file.replace('.fz','').replace('.fits','{}.fits'.format(d)),overwrite=True)
    
    if verbose:
        print('header:',header)
    '''
    RA=hdul[0].header['CRVAL1']
    Dec=hdul[0].header['CRVAL2']
    print('RA=%.3f'%RA,'Dec=%.3f'%Dec) #RA and Dec of LCO image
    wcs=WCS(hdul[0].header).celestial
    
    #SDSS
    #pos = coords.SkyCoord('20h53m14.4s +43d44m9.27s', frame='icrs')
    pos=coords.SkyCoord(RA*u.deg,Dec*u.deg,frame='icrs')
    #xid = SDSS.query_region(pos, spectro=True)
    xid=SDSS.query_region(pos,10*u.arcsec)#, spectro=True)
    print(xid)
    #print(SDSS.AVAILABLE_TEMPLATES)
    #template = SDSS.get_spectral_template('qso') #download spectral templates from SDSS

    #create object folder inside data folder if not existent
    if not os.path.exists(path_general+'/archive_lightcurve/%s'%obj):
        os.mkdir(path_general+'/archive_lightcurve/%s'%obj)
    #else:
     #   shutil.rmtree(path_general+'/archive_images/%s'%obj)#,ignore_errors=True)
    
    if xid!=None:
        #sp = SDSS.get_spectra(matches=xid) #only works for spectro=True
        im=SDSS.get_images(matches=xid)#, band='u,g,r,i,z') #download spectra and images for our match
        if os.path.exists(path_general+'/archive_lightcurve/%s/%s_SDSS.fits'%(obj,obj)):
            os.remove(path_general+'/archive_lightcurve/%s/%s_SDSS.fits'%(obj,obj))
        fits.HDUList.writeto(im[0],path_general+'/archive_lightcurve/%s/%s_SDSS.fits'%(obj,obj),overwrite=True)
    
    #SkyView
    #SkyView.list_surveys()
    survey=['SDSSg','SDSSi','SDSSr','SDSSu','SDSSz','DSS']
    paths=[]
    try:
        paths=SkyView.get_images(pos,survey=survey) #searching for and downloading files
    except:
        pass

    if len(paths)!=0:
        for i in range(len(paths)):
            fits.HDUList.writeto(paths[i],path_general+'/archive_lightcurve/%s/%s_%s.fits'%(obj,obj,survey[i]),overwrite=True)
    
    #if no archive images program exits
    if xid==None and len(paths)==0:
                
        print('No archive images found.')
        #exit()
    
    if xid!=None and len(paths)!=0:
        survey=[]    
        #creating list of surveys for use in further functions
        if xid!=None and len(paths)!=0:
            survey=['SDSS','SDSSg','SDSSi','SDSSr','SDSSu','SDSSz','DSS']
        elif xid==None and len(paths)!=0:
            survey=['SDSSg','SDSSi','SDSSr','SDSSu','SDSSz','DSS']
        else:
            survey=['SDSS']
        
        print(survey)
        return survey