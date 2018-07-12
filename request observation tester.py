import requests
from datetime import datetime
from datetime import timedelta
import astropy.units as u
import astropy.coordinates as coord
from astroquery.simbad import Simbad

'''
submit_request.py
Submit an imaging UserRequest to be scheduled for observing.
'''

API_TOKEN = 'API TOKEN'  # API token obtained from https://observe.lco.global/accounts/profile/
PROPOSAL_ID = 'PROPOSAL ID'  # Proposal IDs may be found here: https://observe.lco.global/proposals/

# The target of the observation


RA = input('Right Ascension: \n')
Dec = input('Declination: \n')

def coords2name(RA, Dec):
    """
    RA and Dec need to be input in decimal degrees for LCO, this function
    chooses the name of the nearest object from SIMBAD and uses that for 
    the name and title of the requested observation.
    """
    result_table = Simbad.query_region(coord.SkyCoord(RA, Dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d1m0s')
    name = str(result_table['MAIN_ID'][0])[2:-1]
    
    return name 

target_name = coords2name(RA, Dec)
print(target_name)
 
target = {
    'name': target_name,
    'type': 'SIDEREAL',
    'ra': RA, # want this to be user input 
    'dec': Dec, # want this to be user input
    'epoch': 2000
}

# The configurations for this request. 

expt = 30 # some function goes here. Waiting on some info from Edward.

molecules = [
    {
        'type': 'EXPOSE',
        'instrument_name': '0M4-SCICAM-SBIG',
        'filter': 'v',
        'exposure_time': expt,
        'exposure_count': 1,
    }
]

# The windows in which this request should be considered for observing. 
# The windows are set automatically to be the time and date now and 7 days from now
start = datetime.now()
end = start + timedelta(days=7)

windows = [{
    'start': start.strftime("%Y-%m-%d %H:%M:%S"),
    'end': end.strftime("%Y-%m-%d %H:%M:%S")
}]

# The telescope class that should be used for this observation
location = {
    'telescope_class': '0m4'
}

# Additional constraints to be added to this request
constraints = {
    'max_airmass': 1.6,
    'min_lunar_distance': 30
}


# The full userrequest, with additional meta-data
userrequest = {
    'group_id': target_name,  # The title
    'proposal': PROPOSAL_ID,
    'ipp_value': 1.05,
    'operator': 'SINGLE',
    'observation_type': 'NORMAL',
    'requests': [{
        'target': target,
        'molecules': molecules,
        'windows': windows,
        'location': location,
        'constraints': constraints
    }]
}

# Now that we have a fully formed UserRequest, we can submit it to the api.
response = requests.post(
    'https://observe.lco.global/api/userrequests/',
    headers={'Authorization': 'Token {}'.format(API_TOKEN)},
    json=userrequest  # Make sure you use json!
)

# Make sure this api call was successful
try:
    response.raise_for_status()
except requests.exceptions.HTTPError as exc:
    print('Request failed: {}'.format(response.content))
    raise exc

userrequest_dict = response.json()  # The API will return the newly submitted userrequest as json

# Print out the url on the portal where we can view the submitted request
print('View this observing request: https://observe.lco.global/userrequests/{}/'.format(userrequest_dict['id']))

# Gets an update on request
request_info = requests.get('https://observe.lco.global/api/userrequests/{}/'.format(userrequest_dict['id'])).json()
print(request_info.get('state'))