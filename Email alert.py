# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:37:10 2018

@author: Tilly
"""

import win32com.client as win32

outlook = win32.Dispatch('outlook.application')
mail = outlook.CreateItem(0)
mail.To = 'email address'
mail.Subject = 'EM GW alert'
mail.Body = 'time to upload data set to zooniverse'

mail.Send()
