# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:25:38 2018

@author: Tilly
"""
import requests
import numpy as np
from xml.etree import ElementTree

fromRA = 9.1
fromDec = -4.1
toRA = 10.8
toDec = 3.2

def get_variables(fromRA, toRA, fromDec, toDec):
    '''
    
    '''

    payload = {'fromra':fromRA,
           'tora':toRA, 
           'fromdec':fromDec,
           'todec':toDec}

    response = requests.get('https://www.aavso.org/vsx/index.php?view=api.list', params = payload)

    tree = ElementTree.fromstring(response.content)
    
    ras=[]
    decs=[]
    for variable in tree.findall('VSXObject'):
        ra = variable.find('RA2000').text
        dec = variable.find('Declination2000').text
        ras.append(ra)
        decs.append(dec)
    return ras, decs
    

















































