#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:00:29 2018

@author: lewisprole
"""
import numpy as np
import pandas as pd
import csv
import imcreate
import download_data_LCO
from Request_LCO import *
import time
import threading

 

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
names of the proposals separated by commas 

e.g.
username = 
password = 
datafolder = 
proposals = 

Once the LCO request has been made, the download code will repeat periodically untill there are >30 targets 
in the images & LCO_images folders combined, at which poiint the processing code runs. Set the 'period' variable in seconds.
Once the processing code runs and saves the finished images to Zooniverse upload, it removes the originals 
from the images folder.

'''

''' ENTER USER DATA FOR LCO REQUEST
-------------------------------------------------------------
'''

username= "username"
password = "password"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "proposalID"
magnitude = 9
expt=30
period =10000
path_general='path/to/general_folder'


''' REQUEST LCO IMAGES TO BE TAKEN
-------------------------------------------------------------
'''

#run request
userrequest = Request(RA, Dec,magnitude,PROPOSAL_ID,username,password)
response = requests.post(
    'https://observe.lco.global/api/userrequests/',
    headers = userrequest.retrieve_token(username, password),
    json = userrequest.final_for_submission(RA, Dec, magnitude, PROPOSAL_ID)  
)
try:
    response.raise_for_status()
except requests.exceptions.HTTPError as exc:
    print('Request failed: {}'.format(response.content))
    raise exc
    
userrequest_dict = response.json()  # The API will return the newly submitted userrequest as json

# Print out the url on the portal where we can view the submitted request
print('View this observing request: https://observe.lco.global/userrequests/{}/'.format(userrequest_dict['id']))


''' ENTER DOWNLOAD/PROCESSING DETAILS
-------------------------------------------------------------
e.g. 
sdate='2018-06-19'
edate='2018-08-19'
proposalID="LCOEPO2018A-004"
path_general="/Users/lewisprole/Documents/University/year3/summer_project"
datafolder="/Users/lewisprole/Documents/University/year3/summer_project/LCO_images"
'''


''' DOWNLOAD AND PROCESS IMAGES 
-------------------------------------------------------------
'''

#function downloads all LCO images taken with the propID in given time frame
#checks to see if there are 30 targets in image folder
#runs image processing code if >30 images 
def repeat(period,path_general):
    spectra=False
    download_data_LCO.download(path_general,sdate,edate,proposalID,datafolder)
    imcreate.fold_check(path_general)
    threading.Timer(period, repeat).start()

repeat(period,path_general)
    


