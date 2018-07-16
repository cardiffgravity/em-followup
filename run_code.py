#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:00:29 2018

@author: lewisprole
"""
import numpy as np
import pandas as pd
import csv
import CODE #Tilly's file extraction code 

'''
Lewis' notes
object type can be 'single' or 'multiple'
-if single, set file_directory to 0.
-if multiple, set RA,Dec,mag to 0
-multiple objects in csv file with columns RA,Dec,mag
'''


#simple function allowing code to be run for single target
#^ or csv file containing multiple targets 
file_directory='fake_targets.csv'

def run(object_type,file_directory,RA,Dec,mag):
    
    if object_type=='single':
        CODE(RA,Dec,mag)
        pass
    
    if object_type=='multiple':
        data=pd.read_csv(file_directory)
        for i in range (len(data)):
            RA=data['RA'][i]
            Dec=data['Dec'][i]
            mag=data['mag'][i]
            CODE(RA,Dec,mag)
        pass