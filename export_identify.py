# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 15:00:36 2018

@author: Tilly
"""

import pandas as pd
import json
import numpy as np

data = pd.read_csv('identify-variables-classifications.csv')
useful = data['annotations']

instance1 = useful[0]

instance1 = json.loads(instance1)

instance1T0 = instance1[2]

values = instance1T0['value']

x = y = np.zeros(len(values))

for i in range (len(values)):
    x[i] = values[i]['x']
    y[i] = values[i]['y']
