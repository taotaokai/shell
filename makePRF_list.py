#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
#
import numpy as np
#
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
#
from obspy.taup import TauPyModel
from obspy import UTCDateTime, read, Trace
import pyproj
# 
from utils import *
from taper import cosine_taper
from lanczos_interp1 import lanczos_interp1

#====== parameters
sac_list = sys.argv[1]

tau = 0.2
wl_decon = 0.01

# time window for polarization analysis and deconvolution
# zero time denotes the direct arrival of the main phase
twin_polar = [-15.0, -10.0, 20.0, 25.0]
twin_decon = [-50.0, -40.0, 100.0, 110.0]

cmpnm_list = ['.BHE', '.BHN', '.BHZ']

# main phase name
main_phase_names =  ['P']

# bandpass filter
freqmin = 0.008
freqmax = 2
dt = 0.1
nyqfreq = 0.5/dt

# distance range
min_dist = 30
max_dist = 90

#-- model
#model_file = 'best_fit.mod'
#model_file = 'iasp91_tibet_crust.vel'
model_file = 'iasp91.tvel'

#-- depth migration
mig_dep = np.arange(0, 350, 0.1)

#-- output figure name
#outfig = "PRF_a%.2f.pdf"%(tau)

#====== read sediment model
with open(model_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if not l.startswith('#') ]

# for depth registrated model
lyr_dep = np.array([ float(x[0]) for x in lines ])
lyr_vp = np.array([ float(x[1]) for x in lines ])
lyr_vs = np.array([ float(x[2]) for x in lines ])

nlyr = len(lyr_dep)
lyr_z = np.diff(lyr_dep)
lyr_vp = (lyr_vp[1::] + lyr_vp[0:nlyr-1])/2.0
lyr_vs = (lyr_vs[1::] + lyr_vs[0:nlyr-1])/2.0

# for layered model
#lyr_z = np.array([ float(x[0]) for x in lines ])
#lyr_vp = np.array([ float(x[1]) for x in lines ])
#lyr_vs = np.array([ float(x[2]) for x in lines ])

nlyr = len(lyr_z)

#======
nt = np.ceil((twin_decon[-1]-twin_decon[0])/dt) * 2
nt = int(optimal_fft_size(nt, 2))
times = twin_decon[0] + np.arange(nt)*dt

geod = pyproj.Geod(ellps="WGS84")
taup_model = TauPyModel(model="ak135")
 
rfn_list = []
rayp_list = []
 
with open(sac_list, 'r') as f:
  sacFile_list = [ x.split()[0] for x in f.readlines() if not (x.startswith('#')) ]

for sacFile in sacFile_list:

  print("====== ", sacFile)

  #------ read sac files
  for i in range(len(cmpnm_list)):
    if i == 0:
      st = read(sacFile + cmpnm_list[i])
    else:
      st += read(sacFile + cmpnm_list[i])
  
  #------ filter data
  st.detrend(type='linear')
  st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=10, zerophase=True)
  
  #------ calculate main phase arrival time
  stla = st[0].stats.sac['stla']
  stlo = st[0].stats.sac['stlo']
  evla = st[0].stats.sac['evla']
  evlo = st[0].stats.sac['evlo']
  evdp = st[0].stats.sac['evdp']
  
  az, baz, dist = geod.inv(evlo, evla, stlo, stla)
  
  dist_degree = np.rad2deg(dist/6371000.0)
  
  if (dist_degree < min_dist) or (dist_degree > max_dist): 
    print('[WARN] out of distance range: ', dist_degree)
    continue
  
  arrivals = taup_model.get_travel_times(
      source_depth_in_km=evdp,
      distance_in_degree=dist_degree,
      phase_list=main_phase_names,
      )
  if not arrivals:
    print('[WARN] Failed to calculate travletimes')
    print(evdp, dist_degree, main_phase_names)
    continue
  
  origin_time = st[0].stats.starttime + (st[0].stats.sac['o'] - st[0].stats.sac['b'])
  arrtimes = [ x.time for x in arrivals ]
  idx = np.argmin(arrtimes)
  min_arrtime = arrtimes[idx]
  rayp = arrivals[idx].ray_param_sec_degree
  ta = origin_time + min_arrtime # main phase arrival time

  #------ check if record covers the decon time window
  max_starttime = max([ x.stats.starttime  for x in st])
  min_endtime = min([ x.stats.endtime  for x in st])
  if max_starttime > (ta + min(twin_decon)) or min_endtime < (ta + max(twin_decon)):
    print("[WARN] record does not cover the decon time window.")
    continue
  
  #------ resample data
  ENZ = np.zeros((3,nt))
  for i in range(len(cmpnm_list)):
    tr = st[i]
    tb = tr.stats.starttime
    taper_decon = cosine_taper(tr.times()+(tb-ta), twin_decon)
    ENZ[i,:] = lanczos_interp1(tr.data*taper_decon, tr.stats.delta, times+(ta-tb), na=40)
  
  #------ rotation N,E to R
  RZ = np.zeros((2,nt))
  RZ[0,:] = np.sin(np.deg2rad(baz)-np.pi)*ENZ[0,:] + np.cos(np.deg2rad(baz)-np.pi)*ENZ[1,:]
  RZ[1,:] = ENZ[2,:]
  
  #------ polarizatiion analysis
  taper_polar = cosine_taper(times, twin_polar)
  covMat = np.dot(RZ*taper_polar, np.transpose(RZ*taper_polar))
  w, v = np.linalg.eig(covMat)
  
  #------ decompose R-Z into L-Q coordinate
  if w[1] > w[0]:
    u0 = np.dot(v[:,1]*np.sign(v[1,1]), RZ) # main phase
    u1 = np.dot(v[:,0]*np.sign(v[0,0]), RZ) # converted phase
  else:
    u0 = np.dot(v[:,0]*np.sign(v[1,0]), RZ) # main phase
    u1 = np.dot(v[:,1]*np.sign(v[0,1]), RZ) # converted phase 

  #u0 = RZ[1,:]
  #u1 = RZ[0,:]

  #------ deconvolution
  Fu0 = np.fft.rfft(u0)
  Fu1 = np.fft.rfft(u1)
  
  freq = np.fft.rfftfreq(nt, d=dt)
  
  Fu0_sq = np.abs(Fu0)**2
  wl = wl_decon*np.max(Fu0_sq)
  Fu0_sq[Fu0_sq < wl] = wl
  
  Frfn0 = Fu0*np.conj(Fu0)/Fu0_sq * gauss_spectrum(freq*2.0*np.pi, tau)
  Frfn1 = Fu1*np.conj(Fu0)/Fu0_sq * gauss_spectrum(freq*2.0*np.pi, tau)
  
  rfn0 = np.fft.irfft(Frfn0, nt)
  rfn1 = np.fft.irfft(Frfn1, nt)

  #------
  rayp_list.append(rayp)
  rfn_list.append(rfn1)
  
# #------ ouptut sac
# tr = st[0].copy()
# tr.npts = nt
# tr.stattime = origin_time

# tr.stats.sac['kcmpnm'] = 'SRC'
# tr.data = rfn0
# tr.write(sacFile+'.src', format='SAC')
# 
# tr.stats.sac['kcmpnm'] = 'RFN'
# tr.data = rfn1
# tr.write(sacFile+'.rfn', format='SAC')
  
#====== depth migration
nrfn = len(rfn_list)

# migration depth samples

nmig = len(mig_dep)

rfn_mig = np.zeros((nmig, nrfn))

lyr_dep = np.zeros(nlyr+1)
lyr_dep[1::] = np.cumsum(lyr_z)

for i in range(nrfn):
  rayp_s_km = rayp_list[i]/(np.pi*6371./180.)

  # mark layers with post critical S converted P wave
  qp2 = lyr_vp**(-2) - rayp_s_km**2
  idx = qp2 < 0
  qp2[idx] = np.nan

  vttp = lyr_z * qp2**0.5
  vtts = lyr_z * (lyr_vs**(-2) - rayp_s_km**2)**0.5

  lyr_t0p1s = np.zeros(nlyr+1)
  lyr_t0p1s[1::] = np.cumsum(vtts-vttp)

  mig_time = np.interp(mig_dep, lyr_dep,lyr_t0p1s)
  idx = np.isnan(mig_time)
  mig_time[idx] = 0

  rfn_mig[:,i] = lanczos_interp1(rfn_list[i], dt, mig_time, na=40)
  rfn_mig[idx,i] = np.nan

# stack depth migrated RFN
stk_rfn_mig = np.mean(rfn_mig, axis=-1)

#======  plot rfn
fig = plt.figure(figsize=(11, 8.5)) # US Letter

ax1 = fig.add_axes([0.1, 0.1, 0.35, 0.8])
ax2 = fig.add_axes([0.55, 0.1, 0.35, 0.8])

# plot staked depth-migrated rfn
ax2.plot(stk_rfn_mig, mig_dep, 'k-' )
ax2.set_xlabel('RF')
ax2.set_ylabel('Depth (km)')
ax2.invert_yaxis()
stnm = tr.stats.sac['kstnm'].strip()
ax2.set_title("%s: depth migrate"%(stnm))

# plot rfn
yrange = max(rayp_list) - min(rayp_list)
yspacing = 1.5*yrange/nrfn

rfn_mean_amp = np.median(np.array([ np.max(np.abs(x)) for x in rfn_list ]))
yscale_rfn = yspacing/rfn_mean_amp

for i in range(nrfn):
  ax1.plot(times-twin_decon[0], rayp_list[i] + yscale_rfn*rfn_list[i], 'k-')
  ## Ps time
  #rayp_s_km = rayp_list[i]/(np.pi*6371./180.)
  #vttp = lyr_z * (lyr_vp**(-2) - rayp_s_km**2)**0.5
  #vtts = lyr_z * (lyr_vs**(-2) - rayp_s_km**2)**0.5
  ## surface RF
  #t0p1s_sedi_surf = np.sum(vtts[0:ilyr_sedi] - vttp[0:ilyr_sedi])
  #t0p1s_moho_surf = np.sum(vtts[0:ilyr_moho] - vttp[0:ilyr_moho])
  ## mark Ps times
  #ax1.plot(t0p1s_sedi_surf, rayp_list[i], 'r.')
  #ax1.plot(t0p1s_moho_surf, rayp_list[i], 'b.')

ax1.set_xlim([0, 30])
ax1.set_ylim([min(rayp_list)-0.5*yspacing, max(rayp_list)+0.5*yspacing])
ax1.set_xlabel("Time (sec)")
ax1.set_title("%s: P-RF"%(stnm))
ax1.set_ylabel("Ray parameter (s/deg)")

outfig = "%s_PRF_a%.2f.pdf"%(stnm,tau)
fig.savefig(outfig, format="pdf")
