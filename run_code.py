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
#import CODE #Tilly's file extraction code 

'''
object type can be 'single' or 'multiple'
-if single, set file_directory to 0.
-if multiple, set RA,Dec,mag to 0
-multiple objects in csv file with columns RA,Dec,mag
'''


username= "lewisprole"
password = "Walkingdead66"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "LCOEPO2018A-004"
magnitude = 9
expt=30

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

#enter dwonload/processing details
path_general='/Users/lewisprole/Documents/University/year3/summer_project'
sdate='2018-06-19'
edate='2018-08-19'
proposalID="LCOEPO2018A-004"
path_general="/Users/lewisprole/Documents/University/year3/summer_project"
datafolder="/Users/lewisprole/Documents/University/year3/summer_project/LCO_images"
spectra=False

#download iamges
download_data_LCO.download(sdate,edate,proposalID,datafolder)
#process images
imcreate.fold_check(path_general,path_general+'/images')