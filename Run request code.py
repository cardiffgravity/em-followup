# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 09:57:30 2018

@author: Tilly
"""

from Request_LCO import *


##############################################################################


#RUN THE CODE


username= "username"
password = "password"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "LCOEPO2018A-004"
magnitude = 9
expt=30


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


#function to add targets to a .txt file
#or create new file if none exist
path_to_file = '/Users/lewisprole/Documents/GitHub/em-followup'
def txtfile(path_to_file,RA,Dec):
    if os.path.exists(path_to_file+'/target_list.txt')==False:
        f = open("target_list.txt","w+")
        f.write('name expt filter RA Dec\n')
    else:
        f = open("target_list.txt","a+")
    name = userrequest.coord_to_name(RA, Dec)       
    f.write('%s %f v %s %s \n'%(name, expt, RA, Dec))
    f.close()
    pass

txtfile(path_to_file,RA,Dec)