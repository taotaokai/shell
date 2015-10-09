#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Modify RESP file to reverse the order of the FIR coefficients.

"""
import sys

#====== read command line args
resp_file = sys.argv[1]

#resp_file = "RESP.AH.ANQ.00.BHZ"

#====== read RESP file
with open(resp_file,'r') as f:
    lines = [l.replace('\n','')
            for l in f.readlines() if not(l.startswith('#'))]

#====== parse RESP file
nline = len(lines)
i = 0
while i < nline:
    l = lines[i]
    if not(l.startswith('B054F07')):
        print l 
        i += 1
    else:
        print l
        n_numer = int(l.split()[4])

        i += 1
        l = lines[i]
        print l
        n_denom = int(l.split()[4])
        if n_denom != 0:
            print n_denom
            raise ValueError('#denominators!=0')
        
        i += 1
        if n_numer > 0: 
            numerators = [x.split() for x in lines[i:(i+n_numer)]]
            # reverse the order
            k = 0
            for d in reversed(numerators):
                print '%s %5d %15s %15s' %(d[0], k, d[2], d[3])
                k += 1
            i += n_numer