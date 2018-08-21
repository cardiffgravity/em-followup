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
import os

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
    location = {
        'telescope_class': '0m4'
    }
    
    def __init__(self, RA, Dec, magnitude, PROPOSAL_ID, username, password, scope_filter = 'v'):
        self.RA = RA
        self.Dec = Dec
        self.magnitude = magnitude
        self.scope_filter = scope_filter
        self.PROPOSAL_ID = PROPOSAL_ID
        self.username = username
        self.password = password
    
    def coord_to_name(self, RA, Dec):
        """ 
        RA and Dec need to be input in decimal degrees for LCO, this function
        chooses the name of the nearest object from SIMBAD and uses that for 
        the name and title of the requested observation.
        
        Args:
            RA (float): Right Ascention in decimal degrees
            Dec (float): Declination in decimal degrees
        Returns:
            string: first associated name with the given coordinates
        """
        
        result_table = Simbad.query_region(coord.SkyCoord(RA, Dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d1m0s')
        name = str(result_table['MAIN_ID'][0])
    
        return name 
    
    def exposure_time(self, magnitude):
        """ 
        takes estimated magnitude for target object and
        calculates the exposure time suitable for the 0.4m 
        LCO telescope
        
        Args:
            magnitude (float): the target magnitude
        Returns:
            exposure time (float) specific to the 0.4m telescope
        """
        
        expt = 10**(0.375*(magnitude-6))
        
        return expt
        
    
    def add_target_name(self, RA, Dec):
        """ 
        takes target coordinates in and uses the coords_to_name 
        function to search the name of the target and adds the name 
        into the target dictionary
        
        Args:
            RA (float): Right Ascention in decimal degrees
            Dec (float): Declination in decimal degrees
        Returns:
            dict: target element of LCO request dict now including the RA and
            Dec of the object and its name
        """
        
        target_name = self.coord_to_name(RA, Dec)
        target = {
                'name': target_name,
                'type': 'SIDEREAL',
                'ra': RA,
                'dec': Dec,
                'epoch': 2000
                }
        return target
    
    def add_exposure(self, magnitude):
        """
        Uses the exposure_time function and inputs value into
        molecules dictionary
        Args:
            magnitude (float): the target magnitude
        Returns:
            dict: molecules element of LCO request dict including the newly
            calculated exposure time
        """
        
        expt = self.exposure_time(magnitude)
        
        molecules = [
                {
                        'type': 'EXPOSE',
                        'instrument_name': '0M4-SCICAM-SBIG',
                        'filter': self.scope_filter,
                        'exposure_time': expt,
                        'exposure_count': 1,
                        }
                ]
        
        return molecules

    def final_for_submission(self,RA,Dec,magnitude,PROPOSAL_ID):
        """
        This function produces the dictionary to be submitted to LCO.
        Args:
            RA (float): Right Ascention in decimal degrees
            Dec (float): Declination in decimal degrees
            magnitude (float): the target magnitude
            PROPOSAL_ID (string): LCO proposal ID
        Returns:
            dict: full request for LCO
        """
        
        userrequest = {
            'group_id': self.coord_to_name(RA, Dec),  
            'proposal': PROPOSAL_ID,
            'ipp_value': 1.05,
            'operator': 'SINGLE',
            'observation_type': 'NORMAL',
            'requests': [{
                'target': self.add_target_name(RA, Dec),
                'molecules': self.add_exposure(magnitude),
                'windows': self.windows,
                'location': self.location,
                'constraints': self.constraints
            }]
        }
        return userrequest
    
    def retrieve_token(self, username, password):
        
        """
        This function gets an authentication token from the LCO archive.
        Args:
              username (string): User name for LCO archive
              password (string): Password for LCO archive
        Returns:
              dict: LCO authentication token
        """
        
        response = requests.post(
                    'https://observe.lco.global/api/api-token-auth/',
                    data = {
                            'username': username,
                            'password': password
                            }
                    ).json()
        token = response.get('token')
        
        authorisation = {'Authorization': 'Token ' + token}
        

        return authorisation
        




