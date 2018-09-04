#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HISTORY

2017-07-20 created
"""
import sys
import warnings
import argparse
#
import numpy as np
#
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
#
from obspy import read, Trace, Stream, UTCDateTime
from obspy.io.sac import SACTrace
#import pyproj
# 
#from taper import cosine_taper
from lanczos_interp1 import lanczos_interp1

def is_equal(lst):
  return not lst or [lst[0]]*len(lst) == lst

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
  description='moving window cross-correlation using power-normalized coherence method',
  #epilog="It is assumed that the zero time in sac record corresponds to the cross-correlation zero lag time.",
  )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument("--Cnn", type=str, default=None,
  help="correlation file between north-north component.")

parser.add_argument("--Cee", type=str, default=None,
 help="correlation file between east-east component.")

parser.add_argument("--Cne", type=str, default=None,
 help="correlation file between north-east component.")
 
parser.add_argument("--start_time", type=float, default=0.0, 
  help="starttime of the data to process in second.")

parser.add_argument("--end_time", type=float, default=5, 
  help="endtime of the data to process in second.")

#parser.add_argument("--min_freq", type=float, default=0.01, 
#  help="minimum frequency in Hz")
#
#parser.add_argument("--max_freq", type=float, default=1, 
#  help="maximum frequency in Hz")

parser.add_argument("--aniso_range", metavar='ALPHA', type=float, nargs=3, default=[0.0, 0.5, 0.01],
    help="search range of anisotropy strength (V_fast/V_slow -1) and grid interval")

parser.add_argument("--azimuth_range", metavar='THETA', type=float, nargs=3, default=[0.0, 180, 1],
    help="search range of fast direction measured clockwise from the north and grid interval (in degrees).")

#parser.add_argument("--taper_percentage", type=float, default=0.1,
#    help="Decimal percentage of taper at both ends (ranging from 0. to 0.5)")

parser.add_argument("--out_fig", type=str, default=None,
    help="output figure")

args = parser.parse_args()

#====== read data
print("=======================================")
print("Read in data")
print("=======================================\n")

st = Stream()
st += read(args.Cnn)
st += read(args.Cee)
st += read(args.Cne)

print(st)

if not is_equal([(tr.stats.starttime, tr.stats.delta, tr.stats.npts) for tr in st]):
  print("[ERROR] Not equal time samples!")
  sys.exit(-1)

#====== compute power normalized cross-coherence matrix
print("=======================================")
print("Process each correlation time window")
print("=======================================\n")

#-- sum the two sides of the correlation function  
tr = st[0]
npts = tr.stats.npts
dt = tr.stats.delta
if npts%2 != 1:
  print("[ERROR] npts (%d) must be an odd number!"%(npts))
npos = int((npts-1)/2)

Cnn = st[0].data[npos::-1] + st[0].data[npos:]
Cee = st[1].data[npos::-1] + st[1].data[npos:]
Cne = st[2].data[npos::-1] + st[2].data[npos:]

#times = dt*np.arange(npos+1)
#plt.plot(times, Cnn, times, Cee, times, Cne)
#plt.legend(['Cnn','Cee','Cne'])
#plt.show()

#-- cut time window
idx0 = int(args.start_time/dt)
idx1 = int(args.end_time/dt)
times = np.arange(idx0,idx1) * dt

#print(idx0,idx1)

Cnn_cut = Cnn[idx0:idx1]
Cee_cut = Cee[idx0:idx1]
Cne_cut = Cne[idx0:idx1]

#plt.plot(times, Cnn_cut, times, Cee_cut, times, Cne_cut)
#plt.legend(['Cnn_cut','Cee_cut','Cne_cut'])
#plt.show()

#-- 
alpha_axis = np.arange(args.aniso_range[0], args.aniso_range[1], args.aniso_range[2])
theta_axis = np.arange(args.azimuth_range[0], args.azimuth_range[1], args.azimuth_range[2])

nalpha = len(alpha_axis)
ntheta = len(theta_axis)

#print(alpha_axis, theta_axis)

CC = np.zeros((nalpha, ntheta))

for ialpha in range(nalpha):
  alpha = alpha_axis[ialpha]
  times_stretch = times * (1+alpha)
  Cnn_stretch = lanczos_interp1(Cnn, dt, times_stretch, na=40)
  Cee_stretch = lanczos_interp1(Cee, dt, times_stretch, na=40)
  Cne_stretch = lanczos_interp1(Cne, dt, times_stretch, na=40)

  for itheta in range(ntheta):
    theta = theta_axis[itheta]
    theta_rad = np.deg2rad(theta)
    stheta2 = np.sin(theta_rad)**2
    ctheta2 = np.cos(theta_rad)**2
    sctheta = np.sin(theta_rad)*np.cos(theta_rad)
    Cff = Cee_cut*stheta2 + 2*Cne_cut*sctheta + Cnn_cut*ctheta2 # assumed fast direction
    Css = Cee_stretch*ctheta2 - 2*Cne_stretch*sctheta + Cnn_stretch*stheta2 # assumed slow direction (theta+pi/2)

    CC[ialpha,itheta] = np.sum(Cff*Css)/(np.sum(Cff**2)*np.sum(Css**2))**0.5


indmax = np.unravel_index(np.argmax(CC, axis=None), CC.shape)

alpha_max = alpha_axis[indmax[0]]
theta_max = theta_axis[indmax[1]]

print(alpha_max, theta_max)

fig = plt.figure() # US Letter

h_cs = plt.contourf(theta_axis, alpha_axis, CC)
plt.plot(theta_max, alpha_max, 'ro')
#h_cs2 = plt.contour(theta_axis, alpha_axis, CC)
plt.xlabel('fast direction (degree)')
plt.ylabel('anisotropy strength (Vfast/Vslow-1)')
plt.title('correlation between fast and slow direction')

cbar = plt.colorbar(h_cs)
cbar.ax.set_ylabel('similarity')
# Add the contour line levels to the colorbar
#cbar.add_lines(h_cs2)

#plt.show()
#outfig = "%s/%s_polarization.pdf"%(out_dir,stnm_rec)
fig.savefig(args.out_fig, format="pdf")
