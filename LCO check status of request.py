# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 09:29:20 2018

@author: Tilly
"""

import time
import requests

request_id = 677078

def check_status(request_id):
    """
    uses the request ID to return the state of the request
    
    Example
    input:
    request_id = 677078
    
    output:
    'COMPLETE', 'PENDING', 'CANCELED' or 'EXPIRED'
    """
    request_info = requests.get('https://observe.lco.global/api/userrequests/{}/'.format(request_id)).json()
    return request_info.get('state')


while True:
    status = check_status(request_id)
    if status == 'COMPLETE':
        #RUN CODE TO DOWNLOAD THOSE FILES
        break
    
    elif status == 'PENDING':
        continue #to check every 6 hours
    
    else:
        print('REQUEST NO LONGER ACTIVE')
