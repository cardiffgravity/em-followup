# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 09:22:22 2018

@author: Tilly
"""

def exposure_time(magnitude):
    expt = 10**(0.375*magnitude - 6)
    return expt

