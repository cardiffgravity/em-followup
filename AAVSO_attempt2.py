# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:25:38 2018

@author: Tilly
"""
import requests



    
fromRA = 123
fromDec = 123
toRA = 150
toDec = 140
payload = {'fromra':fromRA,'tora':toRA, 'fromdec':fromDec,'todec':toDec}
web = requests.request('GET','https://www.aavso.org/vsx/index.php?view=api.list', params = payload)


print(web)













































