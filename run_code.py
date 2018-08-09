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
#import CODE #Tilly's file extraction code 

'''
object type can be 'single' or 'multiple'
-if single, set file_directory to 0.
-if multiple, set RA,Dec,mag to 0
-multiple objects in csv file with columns RA,Dec,mag
'''


#simple function allowing code to be run for single target
#^ or csv file containing multiple targets 
#file_directory='fake_targets.csv'
#
#def run(object_type,file_directory,RA,Dec,mag):
#    
#    if object_type=='single':
#        print(RA,Dec,mag)
##        CODE(RA,Dec,mag)
#        pass
#    
#    if object_type=='multiple':
#        data=pd.read_csv(file_directory)
#        for i in range (len(data)):
#            RA=data['RA'][i]
#            Dec=data['Dec'][i]
#            mag=data['mag'][i]
#            print(RA,Dec,mag)
##            CODE(RA,Dec,mag)
#        pass

path_general='/Users/lewisprole/Documents/University/year3/summer_project'
sdate='2018-06-19'
edate='2018-08-19'
proposalID="LCOEPO2018A-004"
path_general="/Users/lewisprole/Documents/University/year3/summer_project"
datafolder="/Users/lewisprole/Documents/University/year3/summer_project/LCO_images"
spectra=False

download_data_LCO.download(sdate,edate,proposalID,datafolder)
#imcreate.fold_check(path_general,path_general+'/images')