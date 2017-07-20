#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Taper functions
"""
import numpy as np

#====== exponential taper
def exp_taper(x, xc, a=0.01, e=0.0):
    """exponential taper at two ends stop,pass, pass,stop
        xc: array-like
            stop,pass[,pass,stop]
        x: scalar or array-like
            sample point(s)
        a: taper amplitude at stop frequencies.
        e: constant amplitude in stop band
    """
    if not 0 < a < 1:
        raise ValueError('attenuation value must be within 0 and 1.')
    b = np.log(a)
    x = np.array(x)
    y = np.ones(len(x))
    nc = len(xc)
    if nc == 2: # one-sided taper
        l = xc[1] - xc[0] # taper width
        if l == 0:
            raise ValueError('pass and stop values cannot be the same.')
        elif l > 0: # tapered at left side
            idx = x < xc[1]
            y[idx] = np.exp(b*(xc[1]-x[idx])**2/l**2)*(1.0-e) + e
            idx = xc[1] >= x; y[idx] = 1.0
        else: # tapered at right side
            idx = xc[1] <= x
            y[idx] = np.exp(b*(x[idx]-xc[1])**2/l**2)*(1.0-e) + e
            idx = x < xc[1]; y[idx] = 1.0
    elif nc == 4: # two-sided taper
        if not (xc[0]<xc[1]<xc[2]<xc[3]):
            raise ValueError('4 cutoff values must be in increasing order.')
        else:
            idx = x < xc[1]
            y[idx] = np.exp(b*(xc[1]-x[idx])**2/l**2)*(1.0-e) + e
            
            idx = (xc[1] <= x) & (x <= xc[2]); y[idx] = 1.0

            idx = xc[2] < x
            y[idx] = np.exp(b*(x[idx]-xc[2])**2/l**2)*(1.0-e) + e
    else:
        raise ValueError('number of cutoff values must be either 2 or 4.')

    return y


#====== cosine taper
def cosine_taper(x, xc):
    """cosine taper at two ends stop,pass, pass,stop
        xc: (array-like)
            stop,pass[,pass,stop]
        x: scalar or array-like
            sample points
    """
    nc = len(xc)
    if np.isscalar(x):
        x = np.array([x, ])
    else:
        x = np.array(x)
    y = np.ones(len(x))

    if nc == 2: # sided taper
        l = xc[1] - xc[0] # taper width
        if l == 0:
            raise ValueError('pass and stop values cannot be the same.')
        elif l > 0: # tapered at left side
            idx = x < xc[0]; y[idx] = 0.0
            idx = (xc[0] <= x) & (x <= xc[1])
            y[idx] = 0.5 - 0.5*np.cos(np.pi*(x[idx] - xc[0])/l)
            idx = x > xc[1]; y[idx] = 1.0
        else: # tapered at right side
            idx = x > xc[0]; y[idx] = 0.0
            idx = (xc[1] <= x) & (x <= xc[0])
            y[idx] = 0.5 + 0.5*np.cos(np.pi*(x[idx] - xc[1])/l)
            idx = x < xc[1]; y[idx] = 1.0
    elif nc == 4: # two sided taper
        if not (xc[0]<xc[1]<xc[2]<xc[3]):
            raise ValueError('4 cutoff values must be in increasing order.')
        else:
            idx = x <= xc[0]; y[idx] = 0.0

            idx = (xc[0] < x) & (x < xc[1])
            y[idx] = 0.5 - 0.5*np.cos(np.pi*(x[idx] - xc[0])/(xc[1]-xc[0]))
            
            idx = (xc[1] <= x) & (x <= xc[2]); y[idx] = 1.0

            idx = (xc[2] < x) & (x < xc[3])
            y[idx] = 0.5 + 0.5*np.cos(np.pi*(x[idx] - xc[2])/(xc[3]-xc[2]))

            idx = x > xc[3]; y[idx] = 0.0
    else:
        raise ValueError('number of cutoff values must be either 2 or 4.')

    # restore return value to scalar when input x is a scalar
    if len(y) == 1: y = y[0]

    return y