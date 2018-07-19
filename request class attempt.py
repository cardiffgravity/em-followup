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
        """ RA and Dec need to be input in decimal degrees for LCO, this function
        chooses the name of the nearest object from SIMBAD and uses that for 
        the name and title of the requested observation."""
        
        result_table = Simbad.query_region(coord.SkyCoord(RA, Dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d1m0s')
        name = str(result_table['MAIN_ID'][0])[2:-1]
    
        return name 
    
    def exposure_time(self, magnitude):
        """ takes estimated magnitude for target object and
        calculates the exposure time suitable for the 0.4m 
        LCO telescope"""
        pass
    
    def add_target_name(self, RA, Dec):
        """ takes target coordinates in and uses the coords_to_name 
        function to search the name of the target and adds the name 
        into the target dictionary"""
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
        uses the exposure_time function and inputs value into
        molecules dictionary
        """
        
        expt = self.exposure_time(magnitude)
        
        expt = 30
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
        # Get LCOGT token:
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
        


##############################################################################


#RUN THE CODE


username= "username"
password = "password"
RA = 299.5917
Dec = 35.2017
PROPOSAL_ID = "PROPID"
magnitude = 9

instance = Request(RA, Dec,magnitude,PROPOSAL_ID,username,password)

response = requests.post(
    'https://observe.lco.global/api/userrequests/',
    headers=instance.retrieve_token(username, password),
    json=instance.final_for_submission(RA,Dec,magnitude,PROPOSAL_ID)  
)

try:
    response.raise_for_status()
except requests.exceptions.HTTPError as exc:
    print('Request failed: {}'.format(response.content))
    raise exc
    
userrequest_dict = response.json()  # The API will return the newly submitted userrequest as json

# Print out the url on the portal where we can view the submitted request
print('View this observing request: https://observe.lco.global/userrequests/{}/'.format(userrequest_dict['id']))

