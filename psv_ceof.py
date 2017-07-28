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
  epilog="z axis points downwards, half space 1 on negative z-axis and half space 2 on positive z-axis",
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

incarg = np.complex(incang)

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

#print(rayp)
#print(M1, Minv1, Q1)
#print(M2, Minv2, Q2)

#====== solve equation that matches the displacements at the interface
if wave_type == 'p1':
  Minv1_Pd2 = np.dot(Minv1, M2[:,0])
  #print(Minv1_Pd2)
  Minv1_Sd2 = np.dot(Minv1, M2[:,2])
  #print(Minv1_Sd2)
  a = Minv1_Pd2[0]
  b = Minv1_Sd2[0]
  c = Minv1_Pd2[2]
  d = Minv1_Sd2[2]
  #print(a,b,c,d)
  det = a*d - c*b
  Tp1p2 = d/det
  Tp1s2 = -c/det
  m1 = Tp1p2*Minv1_Pd2 + Tp1s2*Minv1_Sd2
  Rp1p1 = m1[1]
  Rp1s1 = m1[3]
  print("Tp1p2 Tp1s2 Rp1p1 Rp1s1")
  print(Tp1p2,Tp1s2,Rp1p1,Rp1s1)
elif wave_type == 's1':
  Minv1_Pd2 = np.dot(Minv1, M2[:,0])
  Minv1_Sd2 = np.dot(Minv1, M2[:,2])
  a = Minv1_Pd2[0]
  b = Minv1_Sd2[0]
  c = Minv1_Pd2[2]
  d = Minv1_Sd2[2]
  det = a*d - c*b
  Ts1p2 = -b/det
  Ts1s2 = a/det
  m1 = Tp1p2*Minv1_Pd2 + Tp1s2*Minv1_Sd2
  Rs1p1 = m1[1]
  Rs1s1 = m1[3]
  print("Ts1p2 Ts1s2 Rs1p1 Rs1s1")
  print(Ts1p2,Ts1s2,Rs1p1,Rs1s1)
elif wave_type == 'p2':
  Minv2_Pu1 = np.dot(Minv2, M1[:,1])
  Minv2_Su1 = np.dot(Minv2, M1[:,3])
  a = Minv1_Pu2[1]
  b = Minv1_Su2[1]
  c = Minv1_Pu2[3]
  d = Minv1_Su2[3]
  det = a*d - c*b
  Tp2p1 = d/det
  Tp2s1 = -c/det
  m2 = Tp2p1*Minv2_Pu1 + Tp2s1*Minv2_Su1
  Rp2p2 = m2[0]
  Rp2s2 = m2[2]
  print("Tp2p1 Tp2s1 Rp2p2 Rp2s2")
  print(Tp2p1,Tp2s1,Rp2p2,Rp2s2)
elif wave_type == 's2':
  Minv2_Pu1 = np.dot(Minv2, M1[:,1])
  Minv2_Su1 = np.dot(Minv2, M1[:,3])
  a = Minv1_Pu2[1]
  b = Minv1_Su2[1]
  c = Minv1_Pu2[3]
  d = Minv1_Su2[3]
  det = a*d - c*b
  Ts2p1 = -b/det
  Ts2s1 = a/det
  m2 = Ts2p1*Minv2_Pu1 + Ts2s1*Minv2_Su1
  Rs2p2 = m2[0]
  Rs2s2 = m2[2]
  print("Ts2p1 Ts2s1 Rs2p2 Rs2s2")
  print(Ts2p1,Ts2s1,Rs2p2,Rs2s2)
