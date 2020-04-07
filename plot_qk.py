# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 19:03:03 2020

@author: Lyy
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

df = pd.read_excel("qk.xlsx")
x = np.linspace(0, 149, 150)
y = x * -6 + 1213.26
plt.plot(df['dens'], df['VPHPL'], 'ro')
plt.plot(x, y, 'b-')
plt.show()