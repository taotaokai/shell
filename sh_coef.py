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
import Haskell_isotropy as haskell

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
  description='SH reflection/transmission coefficients at a plane interface.',
  epilog="z axis points downwards, half space 1 on negative z-axis and half space 2 on positive z-axis",
  )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument("--m1", metavar="Vel", type=float, nargs=2, default=[3.5, 2.7],
    help="vs[km/s],rho[g/cm^3] above the interface (negtive z)")

parser.add_argument("--m2", metavar="Vel", type=float, nargs=2, default=[4.5, 3.3],
    help="vs[km/s],rho[g/cm^3] below the interface (positive z)")

parser.add_argument("-t", '--wave-type', type=str, default='s1', choices=['s1','s2'],
    help="incident wave type, s1/s2: down/up-going S")

group = parser.add_mutually_exclusive_group()

group.add_argument("-i", '--incident-angle', type=float, default=20.0,
    help="incident angle in degree")

group.add_argument("-p", '--ray-parameter', type=float, default=None,
    help="ray parameter in s/km")

args = parser.parse_args()

vs1 = args.m1[0]
ro1 = args.m1[1]

vs2 = args.m2[0]
ro2 = args.m2[1]

wave_type = args.wave_type
incident_angle = args.incident_angle
ray_parameter = args.ray_parameter

if incident_angle > 90.0 or incident_angle < 0.0:
  raise ValueError("-i --incident-angle must have a value between 0.0 and 90.0.")

if ray_parameter:
  if ray_parameter < 0.0:
    raise ValueError("-p --ray-parameter must be greater than 0.0.")

#====== Mode vector
# determine ray parameter
if not ray_parameter:
  incang = np.deg2rad(incident_angle)
  if wave_type == 's1':
    ray_parameter = np.sin(incang)/vs1
  elif wave_type == 's2':
    ray_parameter = np.sin(incang)/vs2

M1, Minv1, Q1 = haskell.mode_sh_iso(vs1,ro1,ray_parameter)
M2, Minv2, Q2 = haskell.mode_sh_iso(vs2,ro2,ray_parameter)

#====== solve equation that matches the displacements at the interface
if wave_type == 's1':
  # MinvSu2 * (Sd1 + R*Su1) = 0 
  MinvSu2_Sd1 = np.dot(Minv2[1,:], M1[:,0])
  MinvSu2_Su1 = np.dot(Minv2[1,:], M1[:,1])
  Rs1s1 = -1 * MinvSu2_Sd1/MinvSu2_Su1

  MinvSd2_Sd1 = np.dot(Minv2[0,:], M1[:,0])
  MinvSd2_Su1 = np.dot(Minv2[0,:], M1[:,1])
  Ts1s2 = MinvSd2_Sd1 + Rs1s1*MinvSd2_Su1

  print("Ts1s2 Rs1s1")
  print(Ts1s2, Rs1s1)
  print(np.abs(Ts1s2), np.abs(Rs1s1))
  print(np.angle(Ts1s2,deg=True),np.angle(Rs1s1,deg=True))

elif wave_type == 's2':
  # MinvSd1 * (Su2 + R*Sd2) = 0 
  MinvSd1_Sd2 = np.dot(Minv1[0,:], M2[:,0])
  MinvSd1_Su2 = np.dot(Minv1[0,:], M2[:,1])
  Rs2s2 = -1 * MinvSd1_Su2/MinvSd1_Sd2

  MinvSu1_Sd2 = np.dot(Minv1[1,:], M2[:,0])
  MinvSu1_Su2 = np.dot(Minv1[1,:], M2[:,1])
  Ts2s1 = MinvSu1_Su2 + Rs2s2*MinvSu1_Sd2

  print("Ts2s1 Rs2s2")
  print(Ts2s1, Rs2s2)
  print(np.abs(Ts2s1), np.abs(Rs2s2))
  print(np.angle(Ts2s1,deg=True),np.angle(Rs2s2,deg=True))
