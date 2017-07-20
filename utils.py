#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

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

def gauss_spectrum(freq, tau):
  """ 
  spectrum of the Gaussian STF of unit area: 
    stf(t;t0,tau) = 1/sqrt(PI)/tau * exp(-((t-t0)/tau)^2)
    here, we take t0 = 0, then
    F_stf = exp(- pi^2 * f^2 * tau^2)
    F_stf = exp(- w^2 * tau^2 / 4)
  """
  F_src = np.exp(-np.pi**2 * freq**2 * tau**2)
  #F_src = np.exp(-0.25 * omega**2 * tau**2)
  return F_src
