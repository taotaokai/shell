#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HISTORY

2017-07-20 created
"""
import sys
import argparse
#
import numpy as np
#
import matplotlib
#matplotlib.use("pdf")
import matplotlib.pyplot as plt
#
from obspy import read #, Trace, UTCDateTime
import pyproj
# 
from taper import cosine_taper
from lanczos_interp1 import lanczos_interp1

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
  description='Harmonic analysis of Rayleigh wave horizontal polarization deviation angle using noise cross-correlation data.',
  epilog="It is assumed that the zero time in sac record corresponds to the cross-correlation zero lag time.",
  )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument("-d", "--data-dir", type=str, default="./",
    help="diretory of sac files")

parser.add_argument("-l", "--station-list", type=argparse.FileType('r'), default="station.lst",
    help="text file containing lines of station info: stnm stla stlo")
 
parser.add_argument("-r", "--receiver", type=str, default=None, required=True,
    help="station name as the receiver")

parser.add_argument("-m", "--moveout-velocity", metavar='VEL', type=float, nargs=2, default=[3.0, 4.0],
    help="moveout velocities to determine the begin/end time of the polarization time window, in km/s")

parser.add_argument("-t", "--taper-width", type=float, default=5,
    help="taper width to extend both sides of the time window determined from moveout velocities, in seconds")

parser.add_argument("-w", "--window-length", type=float, default=10,
    help="time window length at zero distance, in seconds")

parser.add_argument("-c", "--component-name", metavar="COMP", type=str, nargs=3, default=['ZE','ZN','ZZ'],
    help="three component names to be appended to the sac files in the sac_list")

parser.add_argument("-f", "--filter", metavar="FREQ", type=float, nargs=2, default=[0.01, 1],
    help="bandpass filter frequency range, in Hz")

parser.add_argument("-g", "--dist-range", metavar="DIST", type=float, nargs=2, default=[0, 200],
    help="epi-distance range to be used, in km")

parser.add_argument("--sampling-interval", type=float, default=0.1,
    help="time interval for resampling data, in second")

parser.add_argument("--plot", action='store_true',
    help="whether to plot individial seismograms")

args = parser.parse_args()

station_list = args.station_list
data_dir = args.data_dir
stnm_rec = args.receiver

min_moveout_vel = args.moveout_velocity[0]
max_moveout_vel = args.moveout_velocity[1]
taper_width = args.taper_width
window_length = args.window_length

cmpnm_list = args.component_name

freqmin = args.filter[0]
freqmax = args.filter[1]

min_dist = args.dist_range[0]
max_dist = args.dist_range[1]

dt = args.sampling_interval

flag_plot = args.plot

nyqfreq = 0.5/dt
if freqmax > 0.9*nyqfreq:
  raise Exception("sampling interval too large!")

#====== read station list
#with open(station_list, 'r') as f:
lines = [ l.split() for l in station_list if not l.startswith('#') ]

stnm_list = [ l[0] for l in lines ]
stla_list = [ l[1] for l in lines ]
stlo_list = [ l[2] for l in lines ]

irec = stnm_list.index(stnm_rec)
stlo_rec = stlo_list[irec]
stla_rec = stla_list[irec]

#======
min_time = min_dist/max_moveout_vel - taper_width
max_time = max_dist/min_moveout_vel + window_length + taper_width

if min_time < 0: 
  min_time = 0

nt = int(np.ceil((max_time - min_time)/dt))
times = min_time + np.arange(nt)*dt

geod = pyproj.Geod(ellps="WGS84")
 
az_list = []
polar_angle_list = []
stnm_src_list = []

for istn in range(len(stnm_list)):

  if istn == irec: continue
  stnm_src = stnm_list[istn]
  stlo_src = stlo_list[istn]
  stla_src = stla_list[istn]

  print("====== ", stnm_src)

  #------ read sac files
  for i in range(len(cmpnm_list)):
    sacfile = "%s/COR_%s.%s.%s.SAC"%(data_dir,stnm_src,stnm_rec,cmpnm_list[i])
    if i == 0:
      st = read(sacfile)
    else:
      st += read(sacfile)
  
  #------ filter data
  st.detrend(type='linear')
  st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=10, zerophase=True)
  
  #------ calculate main phase arrival time
  az, baz, dist = geod.inv(stlo_src, stla_src, stlo_rec, stla_rec)
  dist_km = dist/1000.0
  print(dist_km)

  if dist_km < min_dist or dist_km > max_dist:
    print("[WARN] SKIP source station %s out of required distance range." % (stnm_src))
    continue

  #------ resample data
  ENZ = np.zeros((3,nt))
  for i in range(len(cmpnm_list)):
    tr = st[i]
    tb = tr.stats.sac['b'] # assume zero time corresponds to zero lag correlation
    ENZ[i,:] = lanczos_interp1(tr.data, tr.stats.delta, times-tb, na=40)
  
  #------ polarizatiion analysis
  min_time1 = dist_km/max_moveout_vel
  min_time0 = min_time1 - taper_width
  if min_time0 < 0:
    min_time0 = 0

  max_time1 = dist_km/min_moveout_vel + window_length
  max_time0 = max_time1 + taper_width

  twin_use = np.array([min_time0, min_time1, max_time1, max_time0])
  print(twin_use)

  min_record_time = max([tr.stats.sac['b'] for tr in st])
  max_record_time = min([tr.stats.sac['e'] for tr in st])

  if min(twin_use) < min_record_time or max(twin_use) > max_record_time:
    print("[WARN] SKIP source station %s does not record the analysis time window." % (stnm_src))
    continue

  taper_polar = cosine_taper(times, twin_use)
  #covMat = np.dot(ENZ*taper_polar, np.transpose(ENZ*taper_polar))
  #w, v = np.linalg.eig(covMat)

  ENZ_taper = ENZ*taper_polar

  covMat = np.dot(ENZ_taper[0:2,:], np.transpose(ENZ_taper[0:2,:]))
  w, v = np.linalg.eig(covMat)

  #------ calculate rotation angle between polarization (E,N) and radial direction
  radial = np.array([np.sin(np.deg2rad(baz-180.0)), np.cos(np.deg2rad(baz-180.0)), 0.0])

  if w[1] > w[0]:
    polar1 = np.array([v[0,1], v[1,1], 0.0])
  else:
    polar1 = np.array([v[0,0], v[1,0], 0.0])

  vcross = np.cross(radial, polar1)
  vdot = np.dot(radial, polar1)
  if vdot >= 0.0:
    ang = -1*np.arccos(np.abs(vdot))*np.sign(vcross[2])
  else:
    ang = np.arccos(np.abs(vdot))*np.sign(vcross[2])

  az_list.append((baz-180.0)%360.0)
  polar_angle_list.append(np.rad2deg(ang))
  stnm_src_list.append(stnm_src)

  if flag_plot:
    fig = plt.figure(figsize=(5,9)) # US Letter

    ax = fig.add_subplot(211)
    ax.plot(times, ENZ[0,:], '--', times, ENZ_taper[0,:], '-')
    ax.plot(times, ENZ[1,:], '--', times, ENZ_taper[1,:], '-')
    ax.plot(times, ENZ[2,:], '--', times, ENZ_taper[2,:], '-')
    ax.legend(['E','N','Z'])

    ax = fig.add_subplot(212, projection='3d')

    ax.plot(ENZ_taper[0,:],ENZ_taper[1,:],ENZ_taper[2,:])

    ax.plot(ENZ_taper[0,:],ENZ_taper[1,:],-1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[2,:])))
    ax.plot(ENZ_taper[0,:],1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[1,:])), ENZ_taper[2,:])
    ax.plot(-1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[0,:])),ENZ_taper[1,:],ENZ_taper[2,:])
    
    maxamp_xy = np.max((ENZ_taper[0,:]**2 + ENZ_taper[1,:]**2)**0.5)
    radial_x = np.array([-1, 1])*np.sin(np.deg2rad(baz))*maxamp_xy
    radial_y = np.array([-1, 1])*np.cos(np.deg2rad(baz))*maxamp_xy
    offset_z = -1.5*np.max(np.abs(ENZ_taper[2,:]))
    ax.plot(radial_x, radial_y,np.array([1,1])*offset_z)

    ax.set_xlabel('E')
    ax.set_ylabel('N')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')

    plt.title(stnm_src)

    plt.show()


#====== calculate Fourier series expansion
az = np.array(az_list)
polar_angle = np.array(polar_angle_list)

# sort az
idx = np.argsort(az)
az = az[idx]
polar_angle = polar_angle[idx]
#print(az)
#print(polar_angle)

# get integration intervals for each az samples
naz = len(az)
az_ext = np.zeros(naz+2)
az_ext[0] = az[-1]-360.0
az_ext[-1] = az[0]+360.0
az_ext[1:naz+1] = az
daz = (az_ext[2:naz+2] - az_ext[0:naz])/2.0
#print(daz)
#print(np.sum(daz))

# Fourier coefficients
nc = 10
c = np.zeros(nc)*1j

az_rad = np.deg2rad(az)
daz_rad = np.deg2rad(daz)
for i in range(nc):
  c[i] = 1/np.pi*np.sum(polar_angle*np.exp(-1j*i*az_rad)*daz_rad)
  if i == 0:
    c[i] = c[i]/2.0

amp = np.abs(c)
phi = np.angle(c,deg=True)

#print(amp)
#print(phi)

#====== plot Fourier coeff.
fig = plt.figure(figsize=(11, 8.5)) # US Letter

plt.subplot(211)
deg = np.arange(0,nc)
plt.plot(deg, amp, 'r-o')
for i in range(len(phi)):
  plt.text(deg[i]+0.1, amp[i], "%.1f"%(amp[i]))
plt.xlim([0, nc-1])
plt.ylabel("Amplitude (deg)")
plt.title("Fourier series of deviation angle of Rayleigh polarization: %s"%(stnm_rec))

plt.subplot(212)
deg1 = deg.copy(); deg1[0] = 1
plt.plot(deg, -phi/deg1, 'r-o')
for i in range(len(phi)):
  plt.text(deg[i]+0.1, -phi[i]/deg1[i], "%.1f"%(-phi[i]/deg1[i]))
plt.xlim([0, nc-1])
plt.ylim([-180, 180])
plt.ylabel("Angle (deg)")
plt.xlabel("Harmonic degree")

#plt.show()
outfig = "%s_fourier.pdf"%(stnm_rec)
fig.savefig(outfig, format="pdf")

#======  plot deviation angle
fig = plt.figure(figsize=(11, 8.5)) # US Letter

plt.plot(az, polar_angle, 'ro')
for i in range(len(az_list)):
  plt.text(az[i]+1, polar_angle[i], stnm_src_list[i])

az1 = np.linspace(0,2*np.pi,100)
polar_angle1 = np.zeros(az1.shape)
for i in range(3):
  plt.plot(np.rad2deg(az1), np.real(c[i]*np.exp(1j*i*az1)), '--')
  polar_angle1 += np.real(c[i]*np.exp(1j*i*az1))

#plt.plot(np.rad2deg(az1), polar_angle1, 'k-')

plt.xlabel('Radial direction (deg)')
plt.ylabel('Deviation angle of Rayleigh polarization (clockwise from radial, deg)')
plt.title('receiver station: %s'%(stnm_rec))
plt.xlim([0,360])
#plt.show()

#if not outfig:
outfig = "%s_polarization.pdf"%(stnm_rec)
fig.savefig(outfig, format="pdf")
