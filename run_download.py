# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:40:47 2019

@author: annac
"""

'''
create txt file called data.txt
write:
LCOGT archive username, password, datafolder, name of the proposals,
objects, filters used, telescope sizes
datafolder is the name of where you want to save files 
names of the proposals, objects, filters and telescopes separated by commas 

e.g.
username= 
password= 
datafolder= 
proposals= 
objects=
filters=
telescopes= 
'''

import LCO_download
import threading
import lightcurve_archives
import argparse
import lightcurve

''' ENTER DOWNLOAD/PROCESSING DETAILS
-------------------------------------------------------------
''' 
'''
parser=argparse.ArgumentParser(prog='Run code',description='Running code for sss17a lightcurve')
parser.add_argument('-s','--start_date',type=str,default='2017-08-17',dest='startdate',help='Start date of observations')
parser.add_argument('-e','--end_date',type=str,default='2017-09-06',dest='enddate',help='End date of observations')
parser.add_argument('-p','--proposal_ID',type=str,default='LCO2017AB-001',dest='proposalID',help='Proposal ID(s) of observations, separate by a comma')
parser.add_argument('-o','--object_name',type=str,default='sss17a',dest='object',help='Object(s) observed, separate by a comma')
parser.add_argument('-f','--filter_used',type=str,default='gp',dest='filter',help='Filter(s) used for observations, separate by a comma')
parser.add_argument('-t','--telescope_used',type=str,default='1m0a',dest='telescope',help='Telescope(s) size used for observations, separate by a comma')
parser.add_argument('-d','--dl',action='store_true',dest='dodl',default=True,help='Set to download files from LCO archive')
args=parser.parse_args()
sdate=args.startdate
edate=args.enddate
proposalID=args.proposalID
obj=args.object
filt=args.filter
tel=args.telescope
dl=args.dodl

#print('Start',sdate,'End',edate,'Proposal',proposalID,'Object',obj,'Filter',filt,'Telescope',tel)
'''

sdate='2017-08-17'
edate='2017-09-06'
proposalID='LCO2017AB-001'
obj='sss17a'
filt='gp'
tel='1m0a'
path_general="/Users/annac/Documents/summer_project"
datafolder="/Users/annac/Documents/summer_project/lightcurve_data"

''' DOWNLOAD AND PROCESS IMAGES 
-------------------------------------------------------------
'''

#function downloads all LCO images taken with the propID in given time frame
#checks to see if there are 10 targets in image folder
#runs lightcurve code if >10 images 
def repeat(path_general):
    spectra=False
    #if dl:
    LCO_download.download(path_general,sdate,edate,proposalID,datafolder,obj,filt,tel)
    #lightcurve_archives.archives(path_general,obj)
    lightcurve.lightcurve(path_general)
    #threading.Timer(10000,repeat(path_general)).start()

repeat(path_general)
