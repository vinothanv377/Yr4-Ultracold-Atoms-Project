# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:32:30 2021

@author: VXV762
"""

import numpy as np
Nt = 100
therm_time = 5
print(np.linspace(0, Nt, int((Nt)/(therm_time)), endpoint=False))
y = np.linspace(0, Nt, int((Nt)/(therm_time)), endpoint=False)
print(np.size(y))