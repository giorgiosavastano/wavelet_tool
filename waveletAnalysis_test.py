import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

__author__ = 'Evgeniya Predybaylo'


# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset
#
# See "http://paos.colorado.edu/research/wavelets/"
# The Matlab code written January 1998 by C. Torrence is modified to Python by Evgeniya Predybaylo, December 2014
#
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
#   changed all "log" to "log2", changed logarithmic axis on GWS to
#   a normal axis.
# ------------------------------------------------------------------------------------------------------------------

# READ THE DATA
sst1 = np.genfromtxt('ainp_4_test.txt', skip_header=1, usecols=1)  # input SST time series)  # input SST time series
tsunami   = 8*60*60+30*60 

#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
variance = np.std(sst1, ddof=1) ** 2
sst = (sst1 - np.mean(sst1)) / np.std(sst1, ddof=1)
n = len(sst)
dt = 30.0
time = np.arange(len(sst)) * dt  + 28800.0 # construct time array 
xlim = ([28800.0 , 36600.0])  # plotting range
pad = 1  # pad the time series with zeroes (recommended)
dj = 0.25  # this will do 4 sub-octaves per octave
s0 = 2 * dt  # this says start at a scale of 6 months
j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72  # lag-1 autocorrelation for red noise background
mother = 'PAUL'  # 'PAUL'

# Wavelet transform:
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum

# Significance levels: (variance=1 for the normalized SST)
signif = wave_signif(([1.0]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
global_ws = variance * (np.sum(power, axis=1) / n)  # time-average over all times
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother)

# Scale-average between El Nino periods of 2--8 years
avg = np.logical_and(scale >= 2, scale < 8)
Cdelta = 0.776  # this is for the MORLET wavelet
scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
scale_avg = power / scale_avg  # [Eqn(24)]
scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother)

#------------------------------------------------------ Plotting

#--- Plot time series
plt.figure(figsize=(18, 9))
plt.subplot(211)
plt.plot(time, sst1)
plt.vlines(tsunami, -0.8 ,0.8, lw=2 )
plt.xlim(xlim[:])
plt.ylim(-0.35,0.35)
plt.xlabel('Time (seconds)')
plt.ylabel('TEC (TECU)')
plt.title('a) AINP TEC values, PRN=4')
plt.grid()
plt.hold(False)

# --- Plot 2--8 yr scale-average time series
#plt.subplot(222)
#plt.plot(time, scale_avg)
#plt.xlim(xlim[:])
#plt.xlabel('Time (year)')
#plt.ylabel('Avg variance (degC^2)')
#plt.title('d) 2-8 yr Scale-average Time Series')
#plt.hold(True)
#plt.plot(xlim, scaleavg_signif + [0, 0], '--')
#plt.hold(False)

#--- Contour plot wavelet power spectrum
plt3 = plt.subplot(212)
# levels = [1.0 / 65536.0, 1.0 / 32768.0, 1.0 / 16384.0, 1.0 / 8192.0, 1.0 / 4096.0, 1.0 / 2048.0, 1.0 / 1024.0, 1.0 / 512.0, \
# 			 1.0 / 256.0, 1.0 / 128.0, 1.0 / 64.0, 1.0 / 32.0, \
# 				1.0 / 16.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0 , \
# 					1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
levels = [1.0 / 64.0, 1.0 / 32.0, \
				1.0 / 16.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0 , \
					1, 2, 4, 8, 16, 32, 64]
CS = plt.contourf(time, period, np.log2(power), len(levels))  #*** or use 'contour'
im = plt.contourf(CS, levels=np.log2(levels))
plt.xlabel('Time (seconds)')
plt.ylabel('Period (seconds)')
plt.title('b) TEC Wavelet Power Spectrum')
plt.xlim(xlim[:])
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
plt.hold(True)
plt.contour(time, period, sig95, [-99, 1], colors='k')
# cone-of-influence, anything "below" is dubious
plt.plot(time, coi, 'k')
plt.hold(False)
# format y-scale
plt3.set_yscale('log', basey=2, subsy=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt3.ticklabel_format(axis='y', style='plain')
plt3.invert_yaxis()
# set up the size and location of the colorbar
divider = make_axes_locatable(plt3)
cax = divider.append_axes("bottom", size="5%", pad=0.5)
plt.colorbar(im, cax=cax, orientation='horizontal')

##--- Plot global wavelet spectrum
#plt4 = plt.subplot(224)
#plt.plot(global_ws, period)
#plt.hold(True)
#plt.plot(global_signif, period, '--')
#plt.hold(False)
#plt.xlabel('Power (degC^2)')
#plt.ylabel('Period (seconds)')
#plt.title('c) Global Wavelet Spectrum')
#plt.xlim([0, 1.25 * np.max(global_ws)])
## format y-scale
#plt4.set_yscale('log', basey=2, subsy=None)
#plt.ylim([np.min(period), np.max(period)])
#ax = plt.gca().yaxis
#ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#plt4.ticklabel_format(axis='y', style='plain')
#plt4.invert_yaxis()

plt.tight_layout()

plt.show()

# end of code
