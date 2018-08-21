# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:25:38 2018

@author: Tilly
"""
import requests

from xml.etree import ElementTree
fromRA = 9.1
fromDec = -4.1
toRA = 10.8
toDec = 3.2

payload = {'fromra':fromRA,
           'tora':toRA, 
           'fromdec':fromDec,
           'todec':toDec}

response = requests.get('https://www.aavso.org/vsx/index.php?view=api.list', params = payload)

tree = ElementTree.fromstring(response.content)

for variable in tree.findall('VSXObject'):
    ra = variable.find('RA2000').text
    dec = variable.find('Declination2000').text
    print(ra, dec)
    

















































