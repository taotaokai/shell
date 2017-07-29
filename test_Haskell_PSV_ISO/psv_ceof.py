#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
#
import numpy as np
#
import matplotlib
#matplotlib.use("pdf")
import matplotlib.pyplot as plt
#
import Haskell_PSV_ISO as haskell

#====== parameters
def convert_arg_line_to_args(arg_line): 
  if not arg_line.startswith("#"):
    for arg in arg_line.split():
      if not arg.strip():
        continue
      yield arg

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  fromfile_prefix_chars='@',
  description='P-SV reflection/transmission coefficients at a plane interface.',
  epilog="z axis points downwards, half space 1 on negative axis and half space 2 on positive axis",
  )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument("--m1", metavar="Vel", type=float, nargs=3, default=[6, 3.5, 2.7],
    help="vp,vs,rho below the interface in km/s,km/s,g/cm^3")

parser.add_argument("--m2", metavar="Vel", type=float, nargs=3, default=[8, 4.5, 3.3],
    help="vp,vs,rho below the interface in km/s,km/s,g/cm^3")

parser.add_argument("-t", '--wave-type', type=str, default='p1',
    help="incident wave type, p1,s1: downgoing P/S, p2,s2: upgoing P/S")

parser.add_argument("-i", '--incident-angle', type=float, default=20.0,
    help="incident angle in degree ")

args = parser.parse_args()

vp1 = args.m1[0]
vs1 = args.m1[1]
ro1 = args.m1[2]

vp2 = args.m2[0]
vs2 = args.m2[1]
ro2 = args.m2[2]

wave_type = args.wave_type
incident_angle = args.incident_angle

#====== Mode vector
# determine ray parameter
incang = np.deg2rad(incident_angle)

if wave_type == 'p1':
  rayp = np.sin(incang)/vp1
elif wave_type == 's1':
  rayp = np.sin(incang)/vs1
elif wave_type == 'p2':
  rayp = np.sin(incang)/vp2
elif wave_type == 's2':
  rayp = np.sin(incang)/vs2

M1, Minv1, Q1 = haskell.mode_psv_iso(vp1,vs1,ro1,rayp)
M2, Minv2, Q2 = haskell.mode_psv_iso(vp2,vs2,ro2,rayp)

print(M1, Minv1, Q1)
