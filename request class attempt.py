# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:29:58 2018

@author: Tilly
"""
import requests
from datetime import datetime
from datetime import timedelta
import astropy.units as u
import astropy.coordinates as coord
from astroquery.simbad import Simbad

class Request():
    
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_time = (datetime.now() + timedelta(days=7)).strftime("%Y-%m-%d %H:%M:%S")
        
    windows = [{
            'start': start_time,
            'end': end_time
            }]
    constraints = {
            'max_airmass': 1.6,
            'min_lunar_distance': 30
            }
    
    def __init__(self, RA, Dec, magnitude):
        self.RA = RA
        self.Dec = Dec
        self.magnitude = magnitude
    
    def coord_to_name(self, RA, Dec):
        """
        RA and Dec need to be input in decimal degrees for LCO, this function
        chooses the name of the nearest object from SIMBAD and uses that for 
        the name and title of the requested observation.
        """
        result_table = Simbad.query_region(coord.SkyCoord(RA, Dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d1m0s')
        name = str(result_table['MAIN_ID'][0])[2:-1]
    
        return name 
    
    def exposure_time(self, magnitude):
        """ takes estimated magnitude for target object and
        calculates the exposure time suitable for the 0.4m 
        LCO telescope"""
        pass

        
#    
#    def add_target(self, target):
#        self(requests)(target) = target
#        return self

test = Request(34,23,10)
help(Request)
