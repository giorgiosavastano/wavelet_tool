import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

# ------------------------------------------------------------------------------------------------------------------

lista = glob.glob( '*_test.txt' )
lista.sort()
dt_list = []
periods_general = []
maxs_general = []

for k in xrange(len(lista)):

	# READ THE DATA
	time = np.genfromtxt( lista[k], skip_header=1, usecols=0 )  # input SST time series)  # input SST time series

	sst1 = np.genfromtxt( lista[k], skip_header=1, usecols=1 )  # input SST time series)  # input SST time series

	#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

	# normalize by standard deviation (not necessary, but makes it easier
	# to compare with plot on Interactive Wavelet page, at
	# "http://paos.colorado.edu/research/wavelets/plot/"
	variance = np.std(sst1, ddof=1) ** 2
	sst = (sst1 - np.mean(sst1)) / np.std(sst1, ddof=1)
	n = len(sst)
	dt = 60.0
	i = 0
	while dt > 30:
	   dt = time[i+1]-time[i]
	   i += 1
	dt_list.append(dt)  
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

	max_list_power = []
	
	for i in xrange( len( power ) ):
	        max_list_power.append( np.max( power[i] ) )

        inde = np.argmax( max_list_power )
        periods_general.append( period[inde] )
        
        
	
