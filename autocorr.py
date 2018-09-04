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
import scipy.signal as signal
#
#import matplotlib
##matplotlib.use("pdf")
#import matplotlib.pyplot as plt
#
from obspy import read, Trace, Stream, UTCDateTime
from obspy.io.sac import SACTrace
#import pyproj
# 
#from taper import cosine_taper
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
  description='moving window cross-correlation using power-normalized coherence method',
  #epilog="It is assumed that the zero time in sac record corresponds to the cross-correlation zero lag time.",
  )

parser.convert_arg_line_to_args = convert_arg_line_to_args

parser.add_argument("--data_list", type=argparse.FileType('r'), default=None,
  help="List of data files. Wildcard could be used, e.g. Data/Sta/*.BH?")
 
parser.add_argument("--start_time", type=str, default="1990-01-01T00:00:00", 
  help="starttime of the data to process.")

parser.add_argument("--end_time", type=str, default="1990-01-02T00:00:00", 
  help="endtime of the data to process.")

parser.add_argument("--window_length", type=float, default=600,
  help="cross-correlation time window length in seconds")

parser.add_argument("--overlap_ratio", type=float, default=0.9,
  help="cc window overlap percentages (ranging from 0 to 1).")

parser.add_argument("--min_freq", type=float, default=0.01, 
  help="minimum frequency in Hz")

parser.add_argument("--max_freq", type=float, default=1, 
  help="maximum frequency in Hz")

parser.add_argument("--sampling_interval", type=float, default=0.1, 
  help="sampling interval in second")

parser.add_argument("--taper_percentage", type=float, default=0.1,
    help="Decimal percentage of taper at both ends (ranging from 0. to 0.5)")

parser.add_argument("--npts_freq_smooth", type=float, default=0,
    help="n points frequency running average")

parser.add_argument("--out_dir", type=str, default=None,
    help="output directory for stacking results")

args = parser.parse_args()

#--- validity check
nyqfreq = 0.5/args.sampling_interval
if args.max_freq > 0.9*nyqfreq:
  raise Exception("sampling interval too large!")



#====== read data
print("=======================================")
print("Read in data")
print("=======================================\n")

#st = read(args.data_file, starttime=UTCDateTime(args.start_time), endtime=UTCDateTime(args.end_time))
data_starttime = UTCDateTime(args.start_time)
data_endtime = UTCDateTime(args.end_time)

st = Stream()
lines = [ l.split()[0] for l in args.data_list if not l.startswith('#') ]
for fname in lines:
  st += read(fname, starttime=data_starttime, endtime=data_endtime)

st.merge(method=1, interpolation_samples=-1, fill_value=None)

station_list = set([tr.stats.station for tr in st])

print(st)

#====== compute power normalized cross-coherence matrix
print("=======================================")
print("Process each correlation time window")
print("=======================================\n")

#--- determine cc windows
cc_window_num = int((data_endtime - data_starttime)/args.window_length/(1-args.overlap_ratio)-1)
times = np.arange(0, args.window_length, args.sampling_interval)
nt = len(times)

#--- 
cc_times = np.arange(-(nt-1), nt)*args.sampling_interval
idx = np.arange(-(nt-1),nt)
idx[idx<0] += (2*nt-1)

# n-point frequency running average
if args.npts_freq_smooth > 0:
  smth_win = np.ones(args.npts_freq_smooth)

cc_window = {}
for iwin in range(cc_window_num):
  #--- cut window
  cc_window_starttime = data_starttime + iwin*args.window_length*(1-args.overlap_ratio)
  cc_window_endtime = cc_window_starttime + args.window_length 

  print("#####", iwin, ':', cc_window_starttime, '--->', cc_window_endtime)
  print("\n")

  time_margin = 5*np.max([tr.stats.delta for tr in st])
  st_win = st.slice(starttime=cc_window_starttime-time_margin, endtime=cc_window_endtime+time_margin)

  #-- remove traces with missing data
  for tr in st_win:
    if (tr.stats.starttime > cc_window_starttime) \
      or (tr.stats.endtime < cc_window_endtime) \
      or np.ma.is_masked(tr.data):
      #warnings.warn("%s: ignore this trace with data gaps!"%(tr.id))
      #print(tr.stats.starttime, cc_window_starttime)
      #print(tr.stats.endtime, cc_window_endtime)
      #print(np.ma.is_masked(tr.data))
      print("[WARN] %s is removed for missing data!"%(tr.id))
      st_win.remove(tr)

  if not st_win:
    print("No traces in this time window.")
    continue

  #--- filter data
  st_win.detrend()
  st_win.taper(args.taper_percentage)
  st_win.filter('bandpass', freqmin=args.min_freq, freqmax=args.max_freq, corners=10, zerophase=True)

  #-- loop each station
  for stnm in station_list:
    print('Station:', stnm)

    if stnm not in cc_window:
      cc_window[stnm] = {}
    cc_station = cc_window[stnm]

    st_select = st_win.select(station=stnm)
    st_select.sort()

    comp_list = [tr.stats.channel for tr in st_select]
    ncomp = len(comp_list)
    data = np.zeros((ncomp,nt))
    #st_select.plot()
    for icomp in range(ncomp):
      tr = st_select[icomp]
      t0 = cc_window_starttime - tr.stats.starttime
      data[icomp,:] = lanczos_interp1(tr.data, tr.stats.delta, t0 + times, na=40)

    #plt.subplot(311); plt.plot(times, data[0,:])
    #plt.subplot(312); plt.plot(times, data[1,:])
    #plt.subplot(313); plt.plot(times, data[2,:])
    #plt.show()

    #- cross-coherence
    #taper_corner = [0, args.taper_width, args.window_length-args.taper_width, args.window_length]
    #taper_func = cosine_taper(times, taper_corner)
    data_fft = np.fft.rfft(data, n=2*nt-1, axis=1)
    data_fft_power = np.sum(np.abs(data_fft)**2, axis=0)
    if args.npts_freq_smooth > 0:
      data_fft_power = signal.convolve(data_fft_power, smth_win, mode='same') / sum(smth_win)

    for icomp in range(ncomp):
      for jcomp in range(icomp,ncomp):
        cc_ij = np.fft.irfft(data_fft[icomp,:]*np.conj(data_fft[jcomp,:])/data_fft_power, n=2*nt-1)
        cc_ij = cc_ij[idx]
        key_comp_ij = (comp_list[icomp], comp_list[jcomp])
        if key_comp_ij not in cc_window[stnm]:
          cc_station[key_comp_ij] = []
        #cc_station[key_comp_ij].append({'window_idx':iwin, 'cc':cc_ij})
        cc_station[key_comp_ij].append(cc_ij)
        #DEBUG
        #print(key_comp_ij)
        #plt.plot(cc_times, cc_ij)
        #plt.show()
        #sys.exit(-1)

    #ax=plt.subplot(331); plt.plot(cc_times, cc_station[(comp_list[0], comp_list[0])][iwin]); ax.set_xlim([-5,5]); plt.title(comp_list[0]); plt.ylabel(comp_list[0])
    #ax=plt.subplot(332); plt.plot(cc_times, cc_station[(comp_list[0], comp_list[1])][iwin]); ax.set_xlim([-5,5]); plt.title(comp_list[1]) 
    #ax=plt.subplot(333); plt.plot(cc_times, cc_station[(comp_list[0], comp_list[2])][iwin]); ax.set_xlim([-5,5]); plt.title(comp_list[2])
    #ax=plt.subplot(335); plt.plot(cc_times, cc_station[(comp_list[1], comp_list[1])][iwin]); ax.set_xlim([-5,5]); plt.ylabel(comp_list[1])
    #ax=plt.subplot(336); plt.plot(cc_times, cc_station[(comp_list[1], comp_list[2])][iwin]); ax.set_xlim([-5,5]); 
    #ax=plt.subplot(339); plt.plot(cc_times, cc_station[(comp_list[2], comp_list[2])][iwin]); ax.set_xlim([-5,5]); plt.ylabel(comp_list[2])
    #plt.show()

#====== stack
print("=======================================")
print("Stack all correlation  windows")
print("=======================================")
print("# comp_ij  nwin")
for stnm in station_list:
  cc_station = cc_window[stnm]
  for comp_ij in cc_station:
    #print(comp_ij)
    cc_comp_ij = cc_station[comp_ij]
    nwin = len(cc_comp_ij)
    print(comp_ij, nwin)
    cc_stack = np.zeros(2*nt-1)
    for cc_ij in cc_comp_ij:
      cc_stack += cc_ij
    cc_stack /= nwin
    #plt.plot(cc_times, cc_stack[comp_ij])
    #plt.show()
    header = {'kstnm':stnm, 'kcmpnm':"%s.%s"%(comp_ij[0],comp_ij[1]), 'delta':args.sampling_interval}
    sac = SACTrace(data=cc_stack, **header)
    sac.reftime = data_starttime
    sac.o = 0.0
    sac.b = -(nt-1)*args.sampling_interval
    sac.write("%s/%s.%s.%s.sac"%(args.out_dir,stnm,comp_ij[0],comp_ij[1]))

  #ax=plt.subplot(331); plt.plot(cc_times, cc_stack[(comp_list[0], comp_list[0])]); ax.set_xlim([-5,5]); plt.title(comp_list[0]); plt.ylabel(comp_list[0])
  #ax=plt.subplot(332); plt.plot(cc_times, cc_stack[(comp_list[0], comp_list[1])]); ax.set_xlim([-5,5]); plt.title(comp_list[1]) 
  #ax=plt.subplot(333); plt.plot(cc_times, cc_stack[(comp_list[0], comp_list[2])]); ax.set_xlim([-5,5]); plt.title(comp_list[2])
  #ax=plt.subplot(335); plt.plot(cc_times, cc_stack[(comp_list[1], comp_list[1])]); ax.set_xlim([-5,5]); plt.ylabel(comp_list[1])
  #ax=plt.subplot(336); plt.plot(cc_times, cc_stack[(comp_list[1], comp_list[2])]); ax.set_xlim([-5,5]); 
  #ax=plt.subplot(339); plt.plot(cc_times, cc_stack[(comp_list[2], comp_list[2])]); ax.set_xlim([-5,5]); plt.ylabel(comp_list[2])
  #plt.show()


##if not outfig:
#outfig = "%s/%s_polarization.pdf"%(out_dir,stnm_rec)
#fig.savefig(outfig, format="pdf")
