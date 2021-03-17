# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:32:30 2021

@author: VXV762
"""

import numpy as np
import random
import numpy as np
#from sklearn.preprocessing import normalize

x1 = np.random.rand(3)
#x = [random.random(),random.random(), random.random()]
norm1 = x1 / np.linalg.norm(x1)
print(norm1[0]**2+norm1[1]**2+norm1[2]**2)
print(x1)
#norm2 = normalize(x[:,np.newaxis], axis=0).ravel()
#print (np.all(norm1 == norm2))
# True
y = [1,2,3,4,5,6]
print(np.average(y))