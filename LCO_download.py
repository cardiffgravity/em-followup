# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:44:02 2019

@author: annac
"""

import os
import time
import sys
import calendar
import requests
import numpy as np


def download_frames(sdate,edate,headers,prop,datafolder,obj,filt,tel):
    """Download files
      This function downloads all the frames for a given range of dates, querying
      50 frames at a time (i.e., if 150 frames have to be downloaded, the process 
      is repeated 3 times, each time downloading 50 frames). This number 
      assumes connections can be as bad as to be able to download only ~1 Mb per 
      minute (each get request shares file urls that last 48 hours only), assuming 
      60 MB frames (worst case scenarios).
 
      It returns the number of total identified frames for the given range and the 
      number of frames downloaded (which is equal to the number of identified frames 
      if no data for that time range was detected on the system).
      Args:
          sdate(time.time): Search for data collected on this date or later
          edate(time.time): Search for data collected before this date
          headers(dict): Authentication token from the LCO archive
          prop(list): List of proposal IDs to search for
          obj(list): List of objects to search for
          filt(list): List of filters used
          tel(list): List of telescope sizes used
          datafolder (string): Directory to put the data
      Returns:
          tuple: list of files found on the archive, list of files actually downloaded
    """
    nidentified=0
    ndownloaded=0
    response=requests.get('https://archive-api.lco.global/frames/?'+
                            'limit=50&'+
                            'RLEVEL=91&'+
                            'start='+sdate+'&'+
                            'end='+edate+'&'+
                            'PROPID='+prop+'&'+
                            'OBJECT='+obj+'&'+
                            'FILTER='+filt+'&'+
                            'TELID='+tel,
                            headers={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/75.0.3770.142 Safari/537.36'}).json()
    print(response)
    frames=response['results']
    print(frames)
    if len(frames)!=0:
        print('\t > Frames identified for the '+sdate+'/'+edate+' period. Checking frames...')
        while True:
            for frame in frames:
                print(frame)
                nidentified+=1
                # Get date of current image frame:
                date=frame['OBJECT']

                # Create new folder with the date if not already there:
                outpath=os.path.join(datafolder,date)
                if not os.path.exists(outpath):
                    os.mkdir(outpath)

                # Check if file is already on disk and that is not a _cat.fits. If not there
                # and is not a _cat.fits, download the file:
                if not os.path.exists(os.path.join(outpath,frame['filename'])) and\
                   '_cat.fits'!=frame['filename'][-9:]:
                    spectra=False
                    if spectra and not frame['filename'].endswith('.tar.gz'):
                        continue
                    print('\t   + File '+frame['filename']+' not found in '+outpath)
                    print('\t     Downloading ...')
                    with open(os.path.join(outpath,frame['filename']),'wb') as f:
                        f.write(requests.get(frame['url']).content)
                    ndownloaded+=1
            if response.get('next'):
                response=requests.get(response['next'],headers=headers).json()
                frames=response['results']
            else:
                break
    return nidentified,ndownloaded


def get_headers_from_token(username,password):
    """
      This function gets an authentication token from the LCO archive.
      Args:
          username (string): User name for LCO archive
          password (string): Password for LCO archive
      Returns:
          dict: LCO authentication token
    """
    # Get LCOGT token:
    response=requests.post('https://archive-api.lco.global/api-token-auth/',
                             data={'username': username,
                                   'password': password}
                             ).json()

    token=response.get('token')

    # Store the Authorization header
    headers={'Authorization': 'Token '+token}
    return headers

def download(path_general,sdate,edate,proposalID,datafolder,obj,filt,tel):
        
    starting_date=sdate
    ending_date=edate
    PROPID=proposalID
    dfolder=datafolder
    OBJECT=obj
    FILTER=filt
    TELID=tel

    print('\n\t ----------------------------------------------')
    print('\t                lcogtDD v.1.2.\n')
    print('\t Author: Nestor Espinoza (nespino@astro.puc.cl)')
    print('\t                         (github@nespinoza)')
    print('\t w/ contributions from: BJ Fulton (bfulton@caltech.edu)')
    print('\t                                  (github@bjfultn)')
    print('\t ----------------------------------------------\n')
    # Check that user input is ok:
    if starting_date is None:
        print('\t lgogtDD input error: Please, insert a starting date from which')
        print('\t                      to download data from. Usage example:\n')
        print('\t                        python download_data -sdate 2016-04-01')
        print('\n')
        sys.exit()

    # Get current date (in order to explore it, we need to leave
    # ending_date=ending_date + 1 day:
    if ending_date is None:
        ending_date=time.strftime("%Y-%m-%d")
        print('\t > Checking data from {} to {}...\n'.format(starting_date,ending_date))
        c_y,c_m,c_d=ending_date.split('-')
        if int(c_d)+1<=calendar.monthrange(int(c_y),int(c_m))[-1]:
            ending_date=c_y+'-'+c_m+'-'+str(int(c_d)+1)
        elif int(c_m)+1<=12:
            ending_date=c_y+'-'+str(int(c_m)+1)+'-01'
        else:
            ending_date=str(int(c_y)+1)+'-01-01'
    else:
        print('\t > Checking data from {} to {}...\n'.format(starting_date,ending_date))
    # Get data from user file:
    f=open(path_general+'/data.txt','r')
    username=(f.readline().split('=')[-1]).split()[0]
    password=(f.readline().split('=')[-1]).split()[0]
    datafolder=(f.readline().split('=')[-1]).split()[0]    #?
    proposals=(f.readline().split('=')[-1]).split(',')     #?
    objects=(f.readline().split('=')[-1]).split(',')       #?
    filters=(f.readline().split('=')[-1]).split(',')       #?
    telescopes=(f.readline().split('=')[-1]).split(',')    #?

    if PROPID is not None:
        proposals=PROPID.split(',')

    if dfolder is not None:
        datafolder=dfolder
        
    if OBJECT is not None:
        objects=OBJECT.split(',')
        
    if FILTER is not None:
        filters=FILTER.split(',')
        
    if TELID is not None:
        telescopes=TELID.split(',')

    print('\t > Proposals from which data will be fetched: {}'.format(' '.join(proposals)))
    for i in range(len(proposals)):
        proposals[i]=proposals[i].split()[0]
    print('\t > Objects in proposals: {}'.format(' '.join(objects)))
    for i in range(len(objects)):
        objects[i]=objects[i].split()[0]
    print('\t > Filters used for observations: {}'.format(' '.join(filters)))
    for i in range(len(filters)):
        filters[i]=filters[i].split()[0]
    print('\t > Telescope diameter used: {}'.format(' '.join(telescopes)))
    for i in range(len(telescopes)):
        telescopes[i]=telescopes[i].split()[0]
    f.close()

    headers=get_headers_from_token(username,password)

    # Get frame names from starting to ending date:
    for prop in proposals:
        c_y,c_m,c_d=starting_date.split('-')
        e_y,e_m,e_d=np.array(ending_date.split('-')).astype('int')
        while True:
            sdate=c_y+'-'+c_m+'-'+c_d
            if int(c_d)+1<=calendar.monthrange(int(c_y),int(c_m))[-1]:
                edate=c_y+'-'+c_m+'-'+str(int(c_d)+1)
            elif int(c_m)+1<=12:
                edate=c_y+'-'+str(int(c_m)+1)+'-01'
            else:
                edate=str(int(c_y)+1)+'-01-01'

            # Download frames in the defined time ranges:
            nidentified,ndownloaded=download_frames(sdate,edate,headers,prop,datafolder,obj,filt,tel)
            if nidentified!=0:
                print('\t   Final count: '+str(nidentified)+' identified frames, downloaded '+
                      str(ndownloaded)+' new ones.')

            # Get next year, month and day to look for. If it matches the user-defined
            # or current date, then we are done:
            c_y,c_m,c_d=edate.split('-')
            if int(c_y)==e_y and int(c_m)==e_m and int(c_d)==e_d:
                break

    print('\n\t Done!\n')