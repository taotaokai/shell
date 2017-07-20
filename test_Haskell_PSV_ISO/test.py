#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
#
import numpy as np
#
import matplotlib
import matplotlib.pyplot as plt
#
import Haskell_PSV_ISO as haskell
#
from obspy import UTCDateTime, read, Trace
# 

def prime_factor(n):
  """
  Prime factorization
  """
  factor = []
  d = 2
  while True:
    while n%d == 0:
      factor.append(d)
      n //= d
    if n == 1:
      return factor
    d += 1
    if d**2 > n:
      factor.append(n)
      return factor

def optimal_fft_size(n, d):
  """
  Return closest upper integer of n whose largest-prime-factor is smaller than d
  """
  if d < 2:
    raise ValueError("d should >= 2, while %d is used" % (d))
  while True:
    factor = prime_factor(n)
    if max(factor) > d:
      n += 1
    else:
      return n

def stf_gauss_spectrum(omega, tau):
  """ 
  spectrum of the Gaussian STF of unit area: 
    stf(t;t0,tau) = 1/sqrt(PI)/tau * exp(-((t-t0)/tau)^2)
    here, we take t0 = 0, then
    F_stf = exp(- pi^2 * f^2 * tau^2)
    F_stf = exp(- w^2 * tau^2 / 4)
  """
  #F_src = np.exp(-np.pi**2 * f**2 * tau**2)
  F_src = np.exp(-0.25 * omega**2 * tau**2)
  return F_src

# parameters
tau = 0.1
dt = 0.01
nt = optimal_fft_size(5000, 10)
alpha = 1 # suppress amplitude by a factor of alpha at the synthetic end time
inc_wave_type = 0 # P=0, SV=1
time_shift = 10.0 # time shift of direct P from the synthetic begin time 

# layered model
z   = np.array( [1.00, 40.00, 0.00] )
vp  = np.array( [4.00,  6.25, 8.27] )
vs  = np.array( [2.00,  3.65, 4.70] )
rho = np.array( [2.45,  2.77, 3.34] )

#z   = np.array( [10.00] )
#vp  = np.array( [8.27]  ) 
#vs  = np.array( [4.70]  )
#rho = np.array( [3.34]  )

#z   = np.array( [30.00, 0.00] )
#vp  = np.array( [ 6.25, 8.27] )
#vs  = np.array( [ 3.65, 4.70] )
#rho = np.array( [ 2.77, 3.34] )

# 
anginc_deg = np.arange(0,30,1)

#====== receiver response for P-wave incidence
beta = -np.log(alpha)/dt/nt

rayp = np.sin(np.deg2rad(anginc_deg))/vp[-1]
nrayp = rayp.size

omega = 2.0*np.pi*np.fft.rfftfreq(nt, d=dt)
nomega = len(omega)

# time domain anti-aliasing
omega1 = omega - 1j*beta

# frequency domain response
F_Vr, F_Vu = haskell.rf_psv_iso(z,vp,vs,rho,omega1,rayp,inc_wave_type)

# phase shift
nz = len(z)
nrayp = rayp.size
Tp = np.zeros((nz, nrayp))
#Ts = np.zeros((nz, nrayp))
for ip in range(nrayp):
  Tp[:,ip] = (vp**(-2) - rayp[ip]**2)**0.5 * z
  #Ts[:,ip] = (vs**(-2) - rayp[ip]**2)**0.5 * z
  phase_factor = np.exp(1j*omega1*(np.sum(Tp[:,ip]) - time_shift))
  F_Vr[:,ip] = F_Vr[:,ip] * phase_factor
  F_Vu[:,ip] = F_Vu[:,ip] * phase_factor

# source spectrum
F_src = stf_gauss_spectrum(omega1, tau)

# inverse FFT
Vr = np.fft.irfft(np.transpose(F_Vr)*F_src, nt, axis=-1)
Vu = np.fft.irfft(np.transpose(F_Vu)*F_src, nt, axis=-1)

#====== plot
times = dt*np.arange(nt) - time_shift

fig = plt.figure()
for ip in range(nrayp):
  amp_factor = 0.02/np.max(np.abs(Vu[ip,:]))
  line_Vr, = plt.plot(times, amp_factor*Vr[ip,:]+rayp[ip], 'r')
  line_Vu, = plt.plot(times, amp_factor*Vu[ip,:]+rayp[ip], 'k')
  #plt.plot(np.sum(Tp, 0, 'r^')
  #plt.plot(Ts, 0, 'b^')
  #plt.plot(3*Tp, 0, 'r^')
  #plt.plot(2*Tp+Ts, 0, 'c^')
  #plt.plot(Tp+2*Ts, 0, 'gv')

plt.title("Receiver site response")
plt.legend([line_Vr, line_Vu], ['Radial', 'Up'])
plt.xlabel("Time after direct P (sec)")
plt.ylabel("Ray parameter (s/km)")
#plt.show()
plt.savefig("receiver_site_response.pdf", format='pdf')

#====== downward continuation
# FFT
F_Vr = np.fft.rfft(Vr, axis=-1)
F_Vu = np.fft.rfft(Vu, axis=-1)

V0 = np.array( np.zeros((4,nomega,nrayp)), dtype='complex' )

V0[0,:,:] = np.transpose(F_Vr)
V0[1,:,:] = -1.0*np.transpose(F_Vu)

# To crust
z1 = np.copy(z[0:2])
z1[-1] = 0.0
V1 = haskell.dc_psv_iso(z1,vp[0:2],vs[0:2],rho[0:2],omega,rayp,V0)

# inverse FFT
Pu = np.fft.irfft(np.transpose(V1[1,:,:]), nt, axis=-1)
Su = np.fft.irfft(np.transpose(V1[3,:,:]), nt, axis=-1)

fig = plt.figure()
for ip in range(nrayp):
  amp_factor = 0.02/np.max(np.abs(Pu[ip,:]))
  line_Pu, = plt.plot(times, amp_factor*Pu[ip,:]+rayp[ip], 'k')
  line_Su, = plt.plot(times, amp_factor*Su[ip,:]+rayp[ip], 'r')

plt.title("Downward continuation to crust")
plt.legend([line_Pu, line_Su], ['up-P', 'up-S'])
plt.xlabel("Time after direct P (sec)")
plt.ylabel("Ray parameter (s/km)")
#plt.show()
plt.savefig("downward_continuation_to_crust.pdf", format='pdf')

# To mantle 
V1 = haskell.dc_psv_iso(z,vp,vs,rho,omega,rayp,V0)

# inverse FFT
Pu = np.fft.irfft(np.transpose(V1[1,:,:]), nt, axis=-1)
Su = np.fft.irfft(np.transpose(V1[3,:,:]), nt, axis=-1)

fig = plt.figure()
for ip in range(nrayp):
  amp_factor = 0.02/np.max(np.abs(Pu[ip,:]))
  line_Pu, = plt.plot(times, amp_factor*Pu[ip,:]+rayp[ip], 'k')
  line_Su, = plt.plot(times, amp_factor*Su[ip,:]+rayp[ip], 'r')

plt.title("Downward continuation to mantle")
plt.legend([line_Pu, line_Su], ['up-P', 'up-S'])
plt.xlabel("Time after direct P (sec)")
plt.ylabel("Ray parameter (s/km)")
#plt.show()
plt.savefig("downward_continuation_to_mantle.pdf", format='pdf')
