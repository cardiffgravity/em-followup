#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:00:29 2018

@author: annac
"""
#import numpy as np
#import pandas as pd
#import csv
import imcreate
import download_data_LCO
#from Request_LCO import *
#import time
import threading
import Archives
import argparse

''' HOW TO USE

1) directories

create a general folder that contains the following empty folders:
    
'images'
'LCO_images'
'fits_images'
'subtraction'
'relative'
'Zooniverse_upload'
'rejects'
'finished_images'

For pre-downloaded images:
in 'images', create a folder per target to store target images
2 comparison images of single target per target folder
images must have same name keyword as folder
e.g. target folder called 'keyword', images within called 'keyword_1.fits', 'keyword_2.fits'
name target folders using actual target name 

2) making requests

enter username, password, RA, DEC, proposal ID, object magnitude and exposure time 
as exampled below.

3) downloading data

-define variable path_general as directory leading to general folder containing 'images','LCO_images', etc.
-also define directory 'datafolder' leading to LCO_images folder
-enter LCO proposal ID and start/end date of download period
-dates in yyyy-mm-dd format
-set the period of how often the number of downloaded images is checked.
-create txt file called userdata.txt in general folder
write:
LCOGT archive username, password, datafolder, and the name of the proposals
datafolder (I think) is the name of where you want to save files 
names of the proposals and objects separated by commas 

e.g.
username = 
password = 
datafolder = 
proposals = 
objects =

Once the LCO request has been made, the download code will repeat periodically until there are >30 targets 
in the images & LCO_images folders combined, at which point the processing code runs. Set the 'period' variable in seconds.
Once the processing code runs and saves the finished images to Zooniverse upload, it removes the originals 
from the images folder.


'''
''' ENTER USER DATA FOR LCO REQUEST
-------------------------------------------------------------
'''
'''
username= "username"
password = "password"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "proposalID"
object_ID = "obj"
magnitude = 9
expt=30
period =10000
path_general='path/to/general_folder'
'''

''' REQUEST LCO IMAGES TO BE TAKEN
-------------------------------------------------------------
'''
'''
#run request
userrequest = Request(RA, Dec,magnitude,PROPOSAL_ID,username,password)
response = requests.post('https://observe.lco.global/api/userrequests/',
    headers = userrequest.retrieve_token(username, password),
    json = userrequest.final_for_submission(RA, Dec, magnitude, PROPOSAL_ID, object_ID))
try:
    response.raise_for_status()
except requests.exceptions.HTTPError as exc:
    print('Request failed: {}'.format(response.content))
    raise exc
    
userrequest_dict = response.json()  # The API will return the newly submitted userrequest as json

# Print out the url on the portal where we can view the submitted request
print('View this observing request: https://observe.lco.global/userrequests/{}/'.format(userrequest_dict['id']))
'''

''' ENTER DOWNLOAD/PROCESSING DETAILS
-------------------------------------------------------------
''' 
'''
parser=argparse.ArgumentParser(prog='Run code',description='Running code for em followup')
parser.add_argument('-s','--start_date',type=str,default='2019-05-01',dest='startdate',help='Start date of observations')
parser.add_argument('-e','--end_date',type=str,default='2019-07-30',dest='enddate',help='End date of observations')
parser.add_argument('-p','--proposal_ID',type=str,default='FTPEPO2017AB-001',dest='proposalID',help='Proposal ID of observations')
parser.add_argument('-o','--object_name',type=str,default='SN2019gmh',dest='object',help='Object observed')
parser.add_argument('-d','--dl',action='store_true',dest='dodl',default=True,help='Set to download files from LCO archive')
args=parser.parse_args()
sdate=args.startdate
edate=args.enddate
proposalID=args.proposalID
obj=args.object
dl=args.dodl

#print('Start',sdate,'End',edate,'Proposal',proposalID,'Object',obj)
'''

#sdate=input('Start date of observations YYYY-MM-DD\n')
#edate=input('End date of observations YYYY-MM-DD\n')
#proposalID=input('Proposal ID\n')
#obj=input('Object name\n')

sdate='2019-01-01'
edate='2019-08-27'
proposalID='FTPEPO2017AB-001'#"LCOEPO2018A-004"#'LCOEPO2019A-005,FTPEPO2014A-003'
obj='SN2019gmh'#'M51'#"SN2019gmh"
path_general="/Users/annac/Documents/summer_project"
datafolder="/Users/annac/Documents/summer_project/LCO_images"
#period=10000


''' DOWNLOAD AND PROCESS IMAGES 
-------------------------------------------------------------
'''

#function downloads all LCO images taken with the propID in given time frame
#checks to see if there are 2 targets in image folder
#runs image processing code if >2 images 
def repeat(path_general):
    spectra=False
    #if dl:
    download_data_LCO.download(path_general,sdate,edate,proposalID,datafolder,obj)
    survey=Archives.archives(path_general,obj)
    #imcreate.fold_check(path_general,survey,obj)
    #threading.Timer(10000,repeat(path_general)).start()

repeat(path_general)
