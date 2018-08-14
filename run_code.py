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

create general folder containing the followuing:
an 'images' folder for pre-downloaded images.
an 'LCO_images' folder where new LCO images are downloaded.
create 'Zooniverse_upload' folder where processed images are saved.
create empty folders called 'fits_images','subtraction','rejects','finished_images'


2) making requests

can request 'single' or 'multiple' targets
-if single, set file_directory to 0.
-if multiple, set RA,Dec,mag to 0
-multiple objects in csv file with columns RA,Dec,mag

3) downloading data

-within 'images' folder, create a new folder per target
-2 images of the target per target file
-naming e.g. if target folder named 'object', images in folder must be named 'object_1' and 'object_2'
-call target folders whatever you like, but no duplicates
-define path_general directory leading to general folder containing 'images'/'LCO_images'
-also define directory 'datafolder' leading to LCO_images folder
-enter LCO proposal ID and start/end date of download period

Once the LCO request has been made, the download code will repeat periodically untill there are >30 targets 
in the images folder, at which poiint the processing code runs. Set the 'period' variable in seconds.
Once the processing code runs and saves the finished images to Zooniverse upload, it removes the originals 
from the images folder.

'''

''' ENTER USER DATA FOR LCO REQUEST
-------------------------------------------------------------
'''

username= "lewisprole"
password = "Walkingdead66"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "LCOEPO2018A-004"
magnitude = 9
expt=30


''' REQUEST LCO IMAGES TO BE TAKEN
-------------------------------------------------------------
'''

#run request
#userrequest = Request(RA, Dec,magnitude,PROPOSAL_ID,username,password)
#response = requests.post(
#    'https://observe.lco.global/api/userrequests/',
#    headers = userrequest.retrieve_token(username, password),
#    json = userrequest.final_for_submission(RA, Dec, magnitude, PROPOSAL_ID)  
#)
#try:
#    response.raise_for_status()
#except requests.exceptions.HTTPError as exc:
#    print('Request failed: {}'.format(response.content))
#    raise exc
#    
#userrequest_dict = response.json()  # The API will return the newly submitted userrequest as json
#
## Print out the url on the portal where we can view the submitted request
#print('View this observing request: https://observe.lco.global/userrequests/{}/'.format(userrequest_dict['id']))


''' SET UP DIRECTORY TO FOLDER CONTAINING 'IMAGES' FOLDER
-------------------------------------------------------------
'''

#enter dOWnload/processing details

sdate='2018-06-19'
edate='2018-08-19'
proposalID="LCOEPO2018A-004"
path_general="/Users/lewisprole/Documents/University/year3/summer_project"
#path_general='/Libraries/Documents/summer_project/em_followup_master'
datafolder="/Users/lewisprole/Documents/University/year3/summer_project/LCO_images"
spectra=False


''' DOWNLOAD AND PROCESS IMAGES 
-------------------------------------------------------------
'''

#function downloads all LCO images taken with the propID in given time frame
#checks to see if there are 30 targets in image folder
#runs image processing code if >30 images 
def repeat(period,path_general):
    download_data_LCO.download(path_general,sdate,edate,proposalID,datafolder)
    imcreate.fold_check(path_general)
    threading.Timer(period, repeat).start()


    


