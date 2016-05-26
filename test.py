# -*- coding: utf-8 -*-
"""
Created on Fri May 27 06:16:05 2016

@author: rainson
"""

import kwant
from scipy.sparse import coo_matrix as coom

sys =  kwant.builder.Builder()
lat = kwant.lattice.square()
sys[[lat(0,0), lat(0,1), lat(1,0), lat(1,1)]] = 0.0
sys[lat.neighbors()] = 1.0
kwant.plot(sys)
#sys = sys.finalized()
