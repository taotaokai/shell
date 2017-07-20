#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
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
import Haskell_PSV_ISO as haskell

#====== parameters
parser = argparse.ArgumentParser(description='RF downward continuation')

parser.add_argument("-d", "--sac-dir", type=str, default="./",
    help="diretory of sac files [./]")

parser.add_argument("-l", "--sac-list", type=argparse.FileType('r'), default="sac.lst",
    help="text file of a list of sac files [sac.lst]")

parser.add_argument("-m", "--model-file", type=argparse.FileType('r'), default='model.lst',
    help="layered model file of lines in: thick(km) vp(km/s) vs(km/s) rho(g/cm^3). [mode.lst] ")

parser.add_argument("-o", "--out-dir", type=str, default='./',
    help="output directory [./]")

parser.add_argument("-a", "--guassian-width", type=float, default=0.1,
    help="gaussian width used for deconvolution [0.1]")

parser.add_argument("-w", "--water-level", type=float, default=0.01,
    help="water level for deconvolution [0.01]")

parser.add_argument("-f", "--filter", metavar="FREQ", type=float, nargs=2, default=[0.01, 1],
    help="bandpass filter frequency range in Hz [0.01 1]")

parser.add_argument("-p", "--phase_name", metavar="PH", type=str, nargs="+", default=["P",],
    help="phase name of the main arrival ['P']")

parser.add_argument("--sac-header", type=str, default=None,
    help="sac header which marks the onset time of the main arrival. If not set, it will be calculated, provided o header is set correctly (not checked by the program)  [None]")

parser.add_argument("-t", "--polarization-twin", metavar='TIME', type=float, nargs=4, default=[-15.0, -10.0, 20.0, 25.0],
    help="time window for polarization decomposition in seconds [-15.0  -10.0  20.0  25.0]")

parser.add_argument("-c", "--decon-twin", metavar='TIME', type=float, nargs=4, default=[-50.0, -40.0, 100.0, 110.0],
    help="time window for deconvolution and downward continuation in seconds [-50.0  -40.0  100.0  110.0]")

parser.add_argument("-n", "--component-name", metavar="COMP", type=str, nargs=3, default=['BHE','BHN','BHZ'],
    help="three component names to be appended to the sac files in the sac_list ['BHE' 'BHN' 'BHZ']")

parser.add_argument("-r", "--dist-range", metavar="DIST", type=float, nargs=2, default=[30, 90],
    help="epi-distance range to be used in degree [30 90]")

parser.add_argument("-s", "--continuation-depth", type=float, default=5,
    help="downward continuation depth in km [5]")

parser.add_argument("--plot-twin", metavar="TIME", type=float, nargs=2, default=[0, 30],
    help="plot time range of RF [0 30]")

parser.add_argument("--outfig", type=str, default="RF_dc.pdf",
    help="output figure name [RF_dc.pdf]")

parser.add_argument("--sampling-interval", type=float, default=0.1,
    help="time interval for resampling data in second [0.1]")

args = parser.parse_args()

sac_dir = args.sac_dir
sac_list = args.sac_list
cmpnm_list = args.component_name

sediment_model_file = args.model_file
dc_depth = args.continuation_depth

dt = args.sampling_interval

tau = args.guassian_width
wl_decon = args.water_level

freqmin = args.filter[0]
freqmax = args.filter[1]

main_phase_names = args.phase_name
sac_header = args.sac_header
min_dist = args.dist_range[0]
max_dist = args.dist_range[1]
twin_polar = args.polarization_twin
twin_decon = args.decon_twin

plot_twin = args.plot_twin
out_dir = args.out_dir
outfig = args.outfig

#DEBUG
#print(args)
#sys.exit()

nyqfreq = 0.5/dt
if freqmax > 0.9*nyqfreq:
  raise Exception("sampling interval too large!")

#====== read sediment model
#with open(sediment_model_file, 'r') as f:
#  lines = [ l.split() for l in f.readlines() if not l.startswith('#') ]
lines = [ l.split() for l in sediment_model_file.readlines() if not l.startswith('#') ]

lyr_z = np.array([ float(x[0]) for x in lines ])
lyr_vp = np.array([ float(x[1]) for x in lines ])
lyr_vs = np.array([ float(x[2]) for x in lines ])
lyr_rho = np.array([ float(x[3]) for x in lines ])

# get the model above the dc_depth
lyr_dep = np.cumsum(lyr_z)
ilyr = next(i for i,v in enumerate(lyr_dep) if v > dc_depth)

dc_lyr_z = lyr_z[0:ilyr+1]
dc_lyr_z[-1] = dc_depth - lyr_dep[ilyr-1]
dc_lyr_vp = lyr_vp[0:ilyr+1]
dc_lyr_vs = lyr_vs[0:ilyr+1]
dc_lyr_rho = lyr_rho[0:ilyr+1]

#DEBUG
#print(dc_lyr_z, dc_lyr_vp, dc_lyr_vs, dc_lyr_rho)
#print(np.cumsum(dc_lyr_z))
#sys.exit()

#======
nt = np.ceil((twin_decon[-1]-twin_decon[0])/dt) * 2
nt = int(optimal_fft_size(nt, 2))
times = twin_decon[0] + np.arange(nt)*dt

geod = pyproj.Geod(ellps="WGS84")
taup_model = TauPyModel(model="ak135")
 
rayp_list = []
rfn_surface_list = []
rfn_dc_list = []

omega = 2.0*np.pi*np.fft.rfftfreq(nt, d=dt)
nomega = omega.size
 
sacfile_list = [ l.split()[0] for l in sac_list.readlines() if not (l.startswith('#')) ]

#print(sacfile_list)
#sys.exit()

for sacfile in sacfile_list:

  print("------ ", sacfile)

  #------ read sac files
  for i in range(len(cmpnm_list)):
    sacfn = "%s/%s.%s"%(sac_dir, sacfile, cmpnm_list[i])
    if i == 0:
      st = read(sacfn)
    else:
      st += read(sacfn)
  
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
    print('[WARN] SKIP out of distance range ', dist_degree)
    continue
  
  arrivals = taup_model.get_travel_times(
      source_depth_in_km=evdp,
      distance_in_degree=dist_degree,
      phase_list=main_phase_names,
      )
  if not arrivals:
    print('[WARN] SKIP failed to calculate travletimes')
    print(evdp, dist_degree, main_phase_names)
    continue
  
  arrtimes = [ x.time for x in arrivals ]
  idx = np.argmin(arrtimes)
  min_arrtime = arrtimes[idx]
  rayp_s_deg = arrivals[idx].ray_param_sec_degree
  rayp_s_km = rayp_s_deg / (6371.0*np.pi/180)

  # main phase arrival time
  if sac_header:
    print("use arrival time in sac header %s"%(sac_header))
    ta = st[0].stats.starttime + (st[0].stats.sac[sac_header] - st[0].stats.sac['b'])
  else:
    print("use calculated arrival time")
    ta = st[0].stats.starttime + (st[0].stats.sac['o'] + min_arrtime - st[0].stats.sac['b'])

  #------ check if record covers the decon time window
  max_starttime = max([ x.stats.starttime  for x in st])
  min_endtime = min([ x.stats.endtime  for x in st])
  if max_starttime > (ta + min(twin_decon)) or min_endtime < (ta + max(twin_decon)):
    #print(ta+min(twin_decon), max_starttime)
    #print(ta+max(twin_decon), min_endtime)
    #sys.exit()
    print("[WARN] SKIP record does not cover the decon time window.")
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

  #------ deconvolution
  Fu0 = np.fft.rfft(u0)
  Fu1 = np.fft.rfft(u1)
  
  Fu0_sq = np.abs(Fu0)**2
  wl = wl_decon*np.max(Fu0_sq)
  Fu0_sq[Fu0_sq < wl] = wl
  
  Frfn1 = Fu1*np.conj(Fu0)/Fu0_sq * gauss_spectrum(omega, tau)
  
  rfn_surface = np.fft.irfft(Frfn1, nt)

  #------ downward continuatin
  V0 = np.array(np.zeros((4,nomega,1)), dtype='complex')
  V0[0,:,0] = np.transpose(np.fft.rfft(RZ[0,:], axis=-1))
  V0[1,:,0] = -1.0*np.transpose(np.fft.rfft(RZ[1,:], axis=-1))

  # beneath sediment layer
  V1 = haskell.dc_psv_iso(dc_lyr_z,dc_lyr_vp,dc_lyr_vs,dc_lyr_rho,omega,rayp_s_km,V0)

  # inverse FFT
  Pu = np.fft.irfft(np.transpose(V1[1,:,0]), nt, axis=-1)
  Su = np.fft.irfft(np.transpose(V1[3,:,0]), nt, axis=-1)

  #------ deconvolution
  Fu0 = np.fft.rfft(Pu)
  Fu1 = np.fft.rfft(Su)
  
  Fu0_sq = np.abs(Fu0)**2
  wl = wl_decon*np.max(Fu0_sq)
  Fu0_sq[Fu0_sq < wl] = wl
  
  #Frfn0 = Fu0*np.conj(Fu0)/Fu0_sq * gauss_spectrum(omega, tau)
  Frfn1 = Fu1*np.conj(Fu0)/Fu0_sq * gauss_spectrum(omega, tau)
  
  #rfn0 = np.fft.irfft(Frfn0, nt)
  rfn_dc = np.fft.irfft(Frfn1, nt)

  #------ store results
  rayp_list.append(rayp_s_deg)
  rfn_surface_list.append(rfn_surface)
  rfn_dc_list.append(rfn_dc)
  
  #------ ouptut sac
  tr = st[0].copy()
  tr.stats.npts = nt
  tr.stats.delta = dt
  tr.stats.sac['delta'] = dt

  t0 = UTCDateTime("2000-01-01:00:00:00")
  tr.stats.starttime = t0
  tr.stats.sac['nzyear'] = 2000
  tr.stats.sac['nzjday'] = 1
  tr.stats.sac['nzhour'] = 0
  tr.stats.sac['nzmin'] = 0
  tr.stats.sac['nzsec'] = 0
  tr.stats.sac['nzmsec'] = 0
  tr.stats.sac['b'] = 0

  tr.stats.sac['kcmpnm'] = 'rf_surf'
  tr.data = rfn_surface
  sacfn = "%s/%s.rfn_surf"%(out_dir,sacfile)
  tr.write(sacfn, format='SAC')
  
  tr.stats.sac['kcmpnm'] = 'rf_dc'
  tr.data = rfn_dc
  sacfn = "%s/%s.rfn_dc"%(out_dir,sacfile)
  tr.write(sacfn, format='SAC')

  tr.stats.starttime = ta + twin_decon[0]
  tr.stats.sac['nzyear'] = ta.year
  tr.stats.sac['nzjday'] = ta.julday
  tr.stats.sac['nzhour'] = ta.hour
  tr.stats.sac['nzmin'] = ta.minute
  tr.stats.sac['nzsec'] = ta.second
  tr.stats.sac['nzmsec'] = int(ta.microsecond/1000)
  tr.stats.sac['b'] = twin_decon[0]

  tr.stats.sac['kcmpnm'] = 'up-P'
  tr.data = Pu
  sacfn = "%s/%s.pu"%(out_dir,sacfile)
  tr.write(sacfn, format='SAC')

  tr.stats.sac['kcmpnm'] = 'up-SV'
  tr.data = Su
  sacfn = "%s/%s.su"%(out_dir,sacfile)
  tr.write(sacfn, format='SAC')
  
# #------ plot
# plt.figure()
# plt.plot(RZ[0,:]*taper_polar, RZ[1,:]*taper_polar, 'k-')
# plt.xlabel('Radial')
# plt.ylabel('Up')
# plt.title('eigen-value: %.1f'%(max(w)/min(w)))
# plt.show()
# 
# # plot
# plt.figure()
# plt.subplot(311)
# plt.plot(times, RZ[0,:], 'r-')
# plt.plot(times, RZ[1,:], 'k-')
# plt.xlim([twin_decon[0], twin_decon[-1]])
# plt.xlabel('time (sec)')
# plt.ylabel('R-Z')
# 
# plt.subplot(312)
# plt.plot(times, u0, 'k-')
# plt.plot(times, u1, 'r-')
# plt.xlim([twin_decon[0], twin_decon[-1]])
# plt.xlabel('time (sec)')
# plt.ylabel('L-Q')
# 
# plt.subplot(313)
# plt.plot(times-twin_decon[0], rfn0, 'r-')
# plt.plot(times-twin_decon[0], rfn1, 'k-')
# plt.xlim([0, 80])
# plt.xlabel('time (sec)')
# plt.ylabel('RF')
# 
# plt.show()

##====== depth migration
#nrfn = len(rfn_surface_list)
#
## migration depth samples
#mig_dep = np.arange(0, 80, 0.1)
#nmig = len(mig_dep)
#
#rfn_mig_surface = np.zeros((nmig, nrfn))
#rfn_mig_dc = np.zeros((nmig, nrfn))
#
#for i in range(nrfn):
#  rayp_s_km = rayp_list[i]/(np.pi*6371./180.)
#  vttp = lyr_z * (lyr_vp**(-2) - rayp_s_km**2)**0.5
#  vtts = lyr_z * (lyr_vs**(-2) - rayp_s_km**2)**0.5
#  lyr_dep = np.zeros(nlyr+1)
#  lyr_t0p1s = np.zeros(nlyr+1)
#  lyr_dep[1::] = np.cumsum(lyr_z)
#  lyr_t0p1s[1::] = np.cumsum(vtts-vttp)
#
#  mig_time = np.interp(mig_dep, lyr_dep,lyr_t0p1s)
#  rfn_mig_surface[:,i] = lanczos_interp1(rfn_surface_list[i], dt, mig_time, na=40)
#
#  t0p1s_dc = lyr_t0p1s[ilyr_dc]
#  mig_time = mig_time - t0p1s_dc
#  rfn_mig_dc[:,i] = lanczos_interp1(rfn_dc_list[i], dt, mig_time, na=40)
#  idx = mig_time <= 0.0
#  rfn_mig_dc[idx,i] = 0.0
#
## stack depth migrated RFN
#stk_rfn_mig_surface = np.mean(rfn_mig_surface, axis=-1)
#stk_rfn_mig_dc = np.mean(rfn_mig_dc, axis=-1)
#
## plot rfn
#fig = plt.figure(figsize=(11, 8.5)) # US Letter
#plt.plot(stk_rfn_mig_surface, mig_dep, 'k-' )
#plt.plot(stk_rfn_mig_dc, mig_dep, 'r-' )
#plt.xlabel('RF')
#plt.ylabel('Depth (km)')
#plt.legend(['surface RF', 'sub-sediment RF'])
#plt.gca().invert_yaxis()
#plt.title(tr.stats.sac['kstnm'])
#plt.savefig(outfig_mig, format='pdf')

#======  plot rfn
fig = plt.figure(figsize=(11, 8.5)) # US Letter
ax1 = fig.add_axes([0.1, 0.1, 0.35, 0.8])
ax2 = fig.add_axes([0.55, 0.1, 0.35, 0.8])

yrange = max(rayp_list) - min(rayp_list)
nrfn = len(rfn_surface_list)
yspacing = 1.5*yrange/nrfn

rfn_surface_stk = rfn_surface_list[0]
rfn_dc_stk = rfn_dc_list[0]

rfn_mean_amp_surface = np.median(np.array([ np.max(np.abs(x)) for x in rfn_surface_list ]))
yscale_rfn_surface = yspacing/rfn_mean_amp_surface

rfn_mean_amp_dc = np.median(np.array([ np.max(np.abs(x)) for x in rfn_dc_list ]))
yscale_rfn_dc = yspacing/rfn_mean_amp_dc

for i in range(nrfn):
  ax1.plot(times-twin_decon[0], rayp_list[i] + yscale_rfn_surface*rfn_surface_list[i], 'k', lw=0.5)

  # Ps time
  #rayp_s_km = rayp_list[i]/(np.pi*6371./180.)
  #vttp = dc_lyr_z * (dc_lyr_vp**(-2) - rayp_s_km**2)**0.5
  #vtts = dc_lyr_z * (dc_lyr_vs**(-2) - rayp_s_km**2)**0.5
  #t0p1s_surface_dc = np.sum(vtts - vttp)

  ax2.plot(times-twin_decon[0], rayp_list[i] + yscale_rfn_dc*rfn_dc_list[i], 'k', lw=0.5)

  ## surface RF
  #t0p1s_sedi_surf = np.sum(vtts[0:ilyr_sedi] - vttp[0:ilyr_sedi])
  #t0p1s_moho_surf = np.sum(vtts[0:ilyr_moho] - vttp[0:ilyr_moho])
  ## sub-sedi RF
  #t0p1s_moho_dc = np.sum(vtts[ilyr_dc:ilyr_moho] - vttp[ilyr_dc:ilyr_moho])
  ## mark Ps times
  #ax1.plot(t0p1s_sedi_surf, rayp_list[i], 'r.')
  #ax1.plot(t0p1s_moho_surf, rayp_list[i], 'b.')
  #ax2.plot(t0p1s_moho_dc, rayp_list[i], 'b.')

  #if i > 0: 
  #   rfn_surface_stk += rfn_surface_list[i]
  #   rfn_dc_stk += rfn_dc_list[i]

#ax1.plot(times-twin_decon[0], np.mean(rayp_list) + yscale*rfn_surface_stk/nrfn, 'r-', lw=2 )
#ax2.plot(times-twin_decon[0], np.mean(rayp_list) + yscale*rfn_dc_stk/nrfn, 'r-', lw=2 )

ax1.grid(color='b', linestyle='-', linewidth=0.2)
ax2.grid(color='b', linestyle='-', linewidth=0.2)

ax1.set_xlim(plot_twin)
ax2.set_xlim(plot_twin)

ax1.set_ylim([min(rayp_list)-1.5*yspacing, max(rayp_list)+1.5*yspacing])
ax2.set_ylim([min(rayp_list)-1.5*yspacing, max(rayp_list)+1.5*yspacing])

ax1.set_xlabel("Time (sec)")
ax2.set_xlabel("Time (sec)")

ax1.set_title("surface P-RF")
ax2.set_title("sub-sediment P-RF (%.1f km)"%(dc_depth))

ax1.set_ylabel("Ray parameter (s/deg)")

fig.savefig(outfig, format="pdf")
