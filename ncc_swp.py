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
matplotlib.use("pdf")
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
  description='Harmonic analysis of Rayleigh wave horizontal polarization deviation angle using noise cross-correlation data. \
  The deviation angle is approximated to N degree by d(az) ~ sum(An*cos(n*(az - phi_n)), n=0,N)',
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

parser.add_argument("--harmonic-degree", type=int, default=10,
    help="number of harmonic degree in Fourier series")

parser.add_argument("--bootstrap-num", type=int, default=100,
    help="number of bootstrap samples, in second")

parser.add_argument("--bootstrap-ratio", type=float, default=0.8,
    help="ratio of data used in each bootstrap sampling, should be between 0.5 and 1.0")

parser.add_argument("-o", "--out-dir", type=str, default='./',
    help="output directory for result figures")

parser.add_argument("--plot", action='store_true',
    help="whether to plot individial seismograms")

parser.add_argument("--figure-title", type=str, default=None,
    help="title in output figure")

#parser.add_argument("--out-file", type=str, default=None,
#    help="output file for Fourier series")

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

bootstrap_num = args.bootstrap_num
bootstrap_ratio = args.bootstrap_ratio
if bootstrap_num <= 1:
  bootstrap_num = 1
  bootstrap_ratio = 0.8 

n_harmonic_degree = args.harmonic_degree

out_dir = args.out_dir
flag_plot = args.plot
figure_title = args.figure_title

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
#min_time = min_dist/max_moveout_vel - taper_width
#max_time = max_dist/min_moveout_vel + window_length + taper_width
#
#if min_time < 0: 
#  min_time = 0
#
#nt = int(np.ceil((max_time - min_time)/dt))
#times = min_time + np.arange(nt)*dt

geod = pyproj.Geod(ellps="WGS84")
 
az_list = []
polar_angle_list = []
stnm_src_list = []
ellipticity_list = []

for istn in range(len(stnm_list)):

  if istn == irec: continue
  stnm_src = stnm_list[istn]
  stlo_src = stlo_list[istn]
  stla_src = stla_list[istn]

  print("------ ", stnm_src)

  #------ read sac files
  flag_read_sac_ok = True
  for i in range(len(cmpnm_list)):
    #sacfile = "%s/COR_%s.%s.%s.SAC"%(data_dir,stnm_src,stnm_rec,cmpnm_list[i])
    #sacfile = "%s/COR_%s.%s.%s.sac"%(data_dir,stnm_src,stnm_rec,cmpnm_list[i])
    sacfile = "%s/%s.%s.%s"%(data_dir,stnm_src,stnm_rec,cmpnm_list[i])
    try:
      if i == 0:
        st = read(sacfile)
      else:
        st += read(sacfile)
    except:
      flag_read_sac_ok = False
      break
  if not flag_read_sac_ok:
    print("[WARN] SKIP %s, error read sac files."%(stnm_src))
    continue
  
  #------ filter data
  st.detrend(type='linear')
  st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=10, zerophase=True)
  
  #------ calculate station distance
  az, baz, dist = geod.inv(stlo_src, stla_src, stlo_rec, stla_rec)
  dist_km = dist/1000.0
  #print(dist_km)

  if dist_km < min_dist or dist_km > max_dist:
    print("[WARN] SKIP source station %s out of required distance range." % (stnm_src))
    continue

  #------ get time window
  # zero time corresponds to zero lag correlation
  min_time1 = dist_km/max_moveout_vel
  min_time0 = min_time1 - taper_width
  if min_time0 < 0:
    min_time0 = 0

  max_time1 = dist_km/min_moveout_vel + window_length
  max_time0 = max_time1 + taper_width

  twin_use = np.array([min_time0, min_time1, max_time1, max_time0])
  #print(twin_use)

  min_record_time = max([tr.stats.sac['b'] for tr in st])
  max_record_time = min([tr.stats.sac['e'] for tr in st])

  if min(twin_use) < min_record_time or max(twin_use) > max_record_time:
    print("[WARN] SKIP source station %s does not record the analysis time window." % (stnm_src))
    continue

  #------ resample data
  times = np.arange(min_record_time, max_record_time, dt)
  nt = len(times)
  ENZ = np.zeros((3,nt))
  for i in range(len(cmpnm_list)):
    tr = st[i]
    tb = tr.stats.sac['b'] # assume zero time corresponds to zero lag correlation
    ENZ[i,:] = lanczos_interp1(tr.data, tr.stats.delta, times-tb, na=40)
  
  #------ polarizatiion analysis
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
    polar2 = np.array([v[0,0], v[1,0], 0.0])
    ellipticity = w[0]/w[1]
  else:
    polar1 = np.array([v[0,0], v[1,0], 0.0])
    polar2 = np.array([v[0,1], v[1,1], 0.0])
    ellipticity = w[1]/w[0]

  vcross = np.cross(polar1, radial)
  vdot = np.dot(radial, polar1)
  if vdot >= 0.0:
    ang = np.arccos(np.abs(vdot))*np.sign(vcross[2])
  else:
    ang = -1*np.arccos(np.abs(vdot))*np.sign(vcross[2])

  az_list.append((baz-180.0)%360.0)

  polar_angle_list.append(np.rad2deg(ang))
  stnm_src_list.append(stnm_src)
  ellipticity_list.append(ellipticity)

  if flag_plot:
    fig = plt.figure(figsize=(5,9)) # US Letter

    ax = fig.add_subplot(211)
    ax.plot(times, ENZ[0,:], 'r-', lw=0.2)
    ax.plot(times, ENZ[1,:], 'b-', lw=0.2)
    ax.plot(times, ENZ[2,:], 'k-', lw=0.2)
    line_e, = ax.plot(times, ENZ_taper[0,:], 'r-', lw=1)
    line_n, = ax.plot(times, ENZ_taper[1,:], 'b-', lw=1)
    line_z, = ax.plot(times, ENZ_taper[2,:], 'k-', lw=1)
    ax.legend([line_e, line_n, line_z], ['E','N','Z'])

    ax = fig.add_subplot(212, projection='3d')

    ax.plot(ENZ_taper[0,:],ENZ_taper[1,:],ENZ_taper[2,:])

    ax.plot(ENZ_taper[0,:],ENZ_taper[1,:],-1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[2,:])))
    ax.plot(ENZ_taper[0,:],1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[1,:])), ENZ_taper[2,:])
    ax.plot(-1.5*np.ones(nt)*np.max(np.abs(ENZ_taper[0,:])),ENZ_taper[1,:],ENZ_taper[2,:])
    
    maxamp_xy = np.max((ENZ_taper[0,:]**2 + ENZ_taper[1,:]**2)**0.5)
    offset_z = -1.5*np.max(np.abs(ENZ_taper[2,:]))

    radial_x = np.array([-1, 1])*np.sin(np.deg2rad(baz))*maxamp_xy
    radial_y = np.array([-1, 1])*np.cos(np.deg2rad(baz))*maxamp_xy
    ax.plot(radial_x, radial_y, np.array([1,1])*offset_z, 'k')

    radial_x = np.array([-1, 1])*polar1[0]*maxamp_xy
    radial_y = np.array([-1, 1])*polar1[1]*maxamp_xy
    ax.plot(radial_x, radial_y, np.array([1,1])*offset_z, 'r')

    radial_x = np.array([-1, 1])*polar2[0]*maxamp_xy*ellipticity
    radial_y = np.array([-1, 1])*polar2[1]*maxamp_xy*ellipticity
    ax.plot(radial_x, radial_y, np.array([1,1])*offset_z, 'b')

    ax.set_xlabel('E')
    ax.set_ylabel('N')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')

    plt.title(stnm_src)

    plt.show()

#====== calculate Fourier series expansion

if not az_list:
  print("[ERROR] no usable data, exit!")
  sys.exit()

az = np.array(az_list)
polar_angle = np.array(polar_angle_list)
ellipticity = np.array(ellipticity_list)
stnm_src = np.array(stnm_src_list)

idx = np.argsort(az)
az = az[idx]
polar_angle = polar_angle[idx]
stnm_src = stnm_src[idx]
ellipticity = ellipticity[idx]
weight = 1.0 - ellipticity

nstn = len(az)
# max harmonic degree
#n_harmonic_degree = 10
# d(az) ~ sum(amp_n*cos(n*(az - phi_n)), n=0,N)
amp = np.zeros((bootstrap_num, n_harmonic_degree))
phi = np.zeros((bootstrap_num, n_harmonic_degree))

for isample in range(bootstrap_num):
  print("bootstrap No. %d"%(isample))
  # randomly select part of the data
  idx = np.random.rand(nstn)<=bootstrap_ratio
  az_sample = az[idx]
  polar_angle_sample = polar_angle[idx]
  stnm_src_sample = stnm_src[idx]
  ellipticity_sample = ellipticity[idx]
  weight_sample = weight[idx]
  # sort list against az
  idx = np.argsort(az_sample)
  az_sample = az_sample[idx]
  polar_angle_sample = polar_angle_sample[idx]
  stnm_src_sample = stnm_src_sample[idx]
  ellipticity_sample = ellipticity_sample[idx]
  # get integration intervals for each az samples
  naz = len(az_sample)
  az_ext = np.zeros(naz+2)
  az_ext[0] = az_sample[-1]-360.0
  az_ext[-1] = az_sample[0]+360.0
  az_ext[1:naz+1] = az_sample
  daz = (az_ext[2:naz+2] - az_ext[0:naz])/2.0
  #print(daz)
  #print(np.sum(daz))
  # Fourier coefficients
  c = np.zeros(n_harmonic_degree)*1j
  az_rad = np.deg2rad(az_sample)
  daz_rad = np.deg2rad(daz)
  sum_weight_daz_rad = np.sum(weight_sample*daz_rad)
  for i in range(n_harmonic_degree):
    c[i] = 2.0*np.sum(polar_angle_sample*np.exp(-1j*i*az_rad)*daz_rad*weight_sample)/sum_weight_daz_rad
    if i == 0: c[i] = c[i]/2.0
  amp[isample,:] = np.abs(c)
  phi[isample,:] = np.angle(c,deg=True)
  phi[isample,0] = (-phi[isample,0])%360
  phi[isample,1::] = (-phi[isample,1::])%360/np.arange(1,n_harmonic_degree)

# statistics
amp_mean = np.mean(amp, axis=0)
phi_mean = np.mean(phi, axis=0)
amp_std = np.std(amp, axis=0)
phi_std = np.std(phi, axis=0)

# write out results
out_file = "%s/%s_fourier.txt"%(out_dir,stnm_rec)
with open(out_file, "w") as f:
  # write out average/std
  f.write("#Statistics\n")
  f.write("#degree amplitude(deg) std cosine_phase_angle(deg) std\n")
  for ideg in range(n_harmonic_degree):
    f.write("%02d  %+8.1f  %+8.1f  %+8.1f  %+8.1f\n"%(
      ideg, amp_mean[ideg], amp_std[ideg], phi_mean[ideg], phi_std[ideg]))
  # write out each results from each sample
  for isample in range(bootstrap_num):
    f.write("#sample %d\n"%(isample))
    f.write("#degree amplitude(deg) cosine_phase_angle(deg)\n")
    for ideg in range(n_harmonic_degree):
      f.write("%02d  %+8.1f  %+8.1f\n"%(ideg, amp[isample,ideg], phi[isample,ideg]))

#====== plot Fourier coeff.
fig = plt.figure(figsize=(11, 8.5)) # US Letter

plt.subplot(211)
deg = np.arange(0,n_harmonic_degree)
plt.plot(deg, amp_mean, 'k-o', markersize=10)
for i in range(n_harmonic_degree):
  plt.plot(deg[i]*np.ones(bootstrap_num), amp[:,i], 'k.', markersize=5)
  plt.text(deg[i]+0.1, amp_mean[i], "%.1f"%(amp_mean[i]), va='center',ha='left',fontsize=18)
plt.xlim([0, n_harmonic_degree-1])
plt.ylabel("Amplitude (deg)")
plt.tick_params(axis='both', which='major', labelsize=20)
plt.title("%s Fourier series of Rayleigh wave polarization deviation anlge: %s"%(figure_title,stnm_rec))

plt.subplot(212)
plt.plot(deg, phi_mean, 'k-o', markersize=10)
for i in range(n_harmonic_degree):
  plt.plot(deg[i]*np.ones(bootstrap_num), phi[:,i], 'k.', markersize=5)
  plt.text(deg[i]+0.1, phi_mean[i], "%.1f"%(phi_mean[i]), va='bottom',ha='left',fontsize=18)
plt.xlim([0, n_harmonic_degree-1])
plt.ylim([0, 360])
plt.ylabel("Angle (deg)")
plt.xlabel("Harmonic degree")
plt.tick_params(axis='both', which='major', labelsize=20)

#plt.show()
outfig = "%s/%s_fourier.pdf"%(out_dir,stnm_rec)
fig.savefig(outfig, format="pdf")

#======  plot deviation angle
fig = plt.figure(figsize=(11, 8.5)) # US Letter

plt.plot(az, polar_angle, 'ko', markersize=10)
for i in range(nstn):
  plt.text(az[i]+1, polar_angle[i], " %s(%.2f)"%(stnm_src[i], weight[i]), 
           va='center',ha='left', fontsize=18)

out_file = "%s/%s_polarization.txt"%(out_dir,stnm_rec)
with open(out_file, "w") as f:
  f.write("#receiver_station: %s\n"%(stnm_rec))
  f.write("#az polar_angle ellipticity weight source_station\n")
  for i in range(len(az_list)):
    f.write("%+9.3f  %+9.3f  %+9.3f  %+9.3f  %s\n"%(az[i], polar_angle[i], ellipticity[i], weight[i], stnm_src[i]))

az1 = np.linspace(0,2*np.pi,100)
polar_angle1 = np.zeros(az1.shape)
for i in range(3):
  #plt.plot(np.rad2deg(az1), np.real(c[i]*np.exp(1j*i*az1)), '-')
  plt.plot(np.rad2deg(az1), amp_mean[i]*np.cos(i*(az1 - np.deg2rad(phi_mean[i]))), '-')
  #polar_angle1 += np.real(c[i]*np.exp(1j*i*az1))

#plt.plot(np.rad2deg(az1), polar_angle1, 'k-')
plt.tick_params(axis='both', which='major', labelsize=20)

plt.xlabel('Radial direction (deg)')
plt.ylabel('Deviation angle of Rayleigh polarization (clockwise from radial, deg)')
plt.title('%s receiver station: %s'%(figure_title, stnm_rec))
plt.xlim([0,360])
#plt.show()

#if not outfig:
outfig = "%s/%s_polarization.pdf"%(out_dir,stnm_rec)
fig.savefig(outfig, format="pdf")
