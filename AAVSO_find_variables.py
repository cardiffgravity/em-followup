# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 09:25:00 2018

@author: Tilly
"""

import requests


def get_chart_id(ra,dec,filter_band='V',field_of_view=5):
    vsp_template = 'https://www.aavso.org/apps/vsp/api/chart/?format=json&fov={}&maglimit=18.5&ra={}&dec={}'
    r = requests.get(vsp_template.format(field_of_view, ra, dec))
    chart_id = r.json()['chartid']
    return chart_id


