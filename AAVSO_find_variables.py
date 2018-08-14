# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 09:25:00 2018

@author: Tilly
"""

import requests



#def get_details(ra,dec,filter_band='V',field_of_view=5):
    '''
    this function finds the chart id of all NON varying stars in the field...
    silly aavso
    '''
#    result = []
#    vsp_template = 'https://www.aavso.org/apps/vsp/api/chart/?format=json&fov={}&maglimit=22&ra={}&dec={}'
#    r = requests.get(vsp_template.format(field_of_view, ra, dec))
#    chart_id = r.json()['chartid']
#    
#    for star in r.json()['photometry']:
#        comparison = {}
#        comparison['auid'] = star['auid']
#        comparison['ra'] = star['ra']
#        comparison['dec'] = star['dec']
#        for band in star['bands']:
#            if band['band'] == filter_band:
#                comparison['vmag'] = band['mag']
#                comparison['error'] = band['error']
#        result.append(comparison)
#    return result, chart_id



'''
trying to write a code to go throught the vsx instead of vsp but different input and
output parameters, much less user friendly, not sure if this is even gonna work
'''


def something(RA, Dec, filter_band = 'V', fov = 7.5):
    vsx_template = 'https://www.aavso.org/vsx/index.php?view=results.get&coords={}&format=d&size={}&unit=2'
    r = requests.get(vsx_template.format())
