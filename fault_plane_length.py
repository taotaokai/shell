#!/usr/bin/env python
# -*- coding: utf-8 -*-

# empirical fault plane size estimated form seismic moment magnitude
import sys
import numpy as np

#====== specify model
Mw = float(sys.argv[1])

M0 = 10**(1.5*Mw + 9.1) # Kanamoori

# density (kg/m^3)
rho = 2.6*10**3

# shear wave velocity (m/s)
vs = 3.2*10**3

# shear modulus (Pa)
mu = rho*vs**2

# fault slip (m)
# empirical relation between slip and Mw
# strike-slip
d = 10**(0.68*Mw - 4.59)

A = M0/mu/d

# width = length/2
length = (2.0*A)**0.5

print("length ~ ", length)
