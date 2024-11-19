#####################################################################
## This code is translated from Roberto Galvan-Madrid's Mathematica
## code.
##
## To do list:
## 		- Add the second velocity component.
## 		- Include models of the other J transitions.
##		- Add the function of loading a FITS cube and fitting all pixels.
##		- Incorporate LVG models (e.g. RADEX).
## 
## -- xlu, Tue Apr  7 17:23:21 EDT 2015
#####################################################################

import os, time, sys
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants
from lmfit import minimize, Parameters, report_fit

start = time.clock()
print 'Start the timer...'

# Define some useful constants first:
c = constants.c.cgs.value # Speed of light (cm/s)
k_B = constants.k_B.cgs.value # Boltzmann coefficient (erg/K)
h = constants.h.cgs.value # Planck constant (erg*s)

def __readascii__(infile):
	"""
	Read in an ASCII file, first column is velocity/frequency axis,
	second column is intensity/brightness temperature.
	Return two numpy arrays.
	"""
	temp = open(infile, 'r')
	text = temp.readlines()
	spec = np.empty(len(text))
	vaxis = np.empty(len(text))
	for line in range(len(text)):
		vaxis[line] = np.float(text[line].split()[0])
		spec[line] = np.float(text[line].split()[1])
	temp.close()
	del temp

	return spec, vaxis

def __ch3cn_init__():
	"""
	Initialize the parameters. No input keywords.
	"""
	ch3cn_info = {}
	ch3cn_info['E'] = [68.866, 76.010, 97.442, 133.157, 183.147, 247.399, 325.899, 418.630, \
	525.571, 646.696] # Upper level energy (K)
	ch3cn_info['frest'] = [220.7472, 220.7430, 220.7302, 220.7090, 220.6792, 220.6410, 220.5944, \
	220.5393, 220.4758, 220.4039] # Rest frequencies (GHz)
	ch3cn_info['gk'] = [4.0, 4.0, 4.0, 8.0, 4.0, 4.0, 8.0, 4.0, 4.0, 8.0] # Degeneracy = g_K*g_nuclear
	ch3cn_info['mu'] = 3.922e-18 # Permanent dipole moment (esu*cm)
	ch3cn_info['A'] = 158.099 # Rotational constant (GHz)
	ch3cn_info['B'] = 9.199 # Another rotational constant (GHz)

	return ch3cn_info

def __gauss_tau__(axis,p):
	"""
	Genenerate a Gaussian model given an axis and a set of parameters.
	For detail see Araya et al. 2005, ApJS, 157, 279.
	p: [T, Ntot, fsky, sigma, K]
	"""
	T= p[0]; Ntot = p[1]; fsky = p[2]; sigma = p[3]; K = p[4]

	phijk = 1/sqrt(2 * pi) / (sigma * 1e9) * np.exp(-0.5 * (axis - fsky)**2 / sigma**2)
	Ajk = (64 * pi**4 * (ch3cn_info['frest'][K] * 1e9)**3 * ch3cn_info['mu']**2 / 3 / h / c**3) * (J**2 - K**2) / (J * (2*J + 1))
	gjk = (2*J + 1) * ch3cn_info['gk'][K]
	Q = 3.89 * T**1.5 / (-1.0 * expm1(-524.8 / T))**2
	Njk = Ntot * (gjk / Q) * exp(-1.0 * ch3cn_info['E'][K] / T)

	tau = (h * c**2 * Njk * Ajk) / (8 * pi * ch3cn_info['frest'][K] * 1e9 * k_B * T) * phijk	

	return tau

def __tau__(Ntot, sigma, T, K):
	"""
	Calculate opacity at the line center
	"""
	phijk = 1/sqrt(2 * pi) / (sigma * 1e9)
	Ajk = (64 * pi**4 * (ch3cn_info['frest'][K] * 1e9)**3 * ch3cn_info['mu']**2 / 3 / h / c**3) * (J**2 - K**2) / (J * (2*J + 1))
	gjk = (2*J + 1) * ch3cn_info['gk'][K]
	Q = 3.89 * T**1.5 / (-1.0 * expm1(-524.8 / T))**2
	Njk = Ntot * (gjk / Q) * exp(-1.0 * ch3cn_info['E'][K] / T)

	tau = (h * c**2 * Njk * Ajk) / (8 * pi * ch3cn_info['frest'][K] * 1e9 * k_B * T) * phijk

	return tau

def __model_11__(params, faxis, spec):
	"""
	Model hyperfine components of CH3CN.
	Then subtract data.
	params: [T, Ntot, fsky, sigma]
	"""
	T = params['T'].value
	Ntot = params['Ntot'].value
	fsky = params['fsky'].value
	sigma = params['sigma'].value

	fsky_k = ch3cn_info['frest'] + (fsky - ch3cn_info['frest'][0])		

	tau = np.zeros(len(faxis))
	for k in arange(Kladder):
		tau += __gauss_tau__(faxis, [T, Ntot, fsky_k[k], sigma, k])
        model = T * (1 - np.exp(-1.0 * tau))
	return model - spec

clickvalue = []
def onclick(event):
	print 'The frequency you select: %f' % event.xdata
	clickvalue.append(event.xdata)

def fit_spec(spec, faxis, Jupp=12, K_fit=7, cutoff=0.009, varyf=2, interactive=True, mode='single'):
	"""
	Fit the hyperfine lines of CH3CN, derive best-fitted Trot and Ntot.
	Input:
		spec: the spectra
		faxis: the frequency axis
		Jupp: J of the upper level
		K_fit: the number of K transitions to fit (e.g. K=0-6, then K_fit=7)
		cutoff: not used...
		varyf: number of channels to vary after you select the Vlsr
		interactive: true or false
		mode: single or double, components along the line of sight (to be done...)
	"""
	# Define the J, K numbers as global variables:
	# (Not recommended by Python experts...)
	global J
	global Kladder
	J = Jupp
	Kladder = K_fit

	if interactive:
		plt.ion()
		f = plt.figure(figsize=(14,8))
		ax = f.add_subplot(111)

	unsatisfied = True
	while unsatisfied:
		if interactive:
			f.clear()
			plt.ion()
			plt.plot(faxis, spec, 'k-', label='Spectrum')
			cutoff_line = [cutoff] * len(faxis)
			cutoff_line_minus = [-1.0*cutoff] * len(faxis)
			plt.plot(faxis, cutoff_line, 'r-')
			plt.plot(faxis, cutoff_line_minus, 'r-')
			plt.xlabel(r'Sky Frequency (GHz)', fontsize=20, labelpad=20)
			plt.ylabel(r'$T_{\nu}$ (K)', fontsize=20)
			plt.text(0.02, 0.92, sourcename, transform=ax.transAxes, color='r', fontsize=15)
			#plt.ylim([-10,60])
			#clickvalue = []
			if mode == 'single':
				cid = f.canvas.mpl_connect('button_press_event', onclick)
				raw_input('Click on the plot to select a Vlsr...')
				#print clickvalue
				if len(clickvalue) >= 1:
					print 'Please select at least one velocity! The last one will be used.'
					vlsr1 = clickvalue[-1]
					vlsr1 = c * (1 - vlsr1/ch3cn_info['frest'][0]) / 1e5
				elif len(clickvalue) == 0:
					vlsr1 = 0.0
				print 'Or input one velocity manually:'
				manualv = raw_input()
				manualv = manualv.split()
				if len(manualv) == 1:
					vlsr1 = np.float_(manualv)
				else:
					print 'Invalid input...'
				print 'The Vlsr is %0.2f km/s' % vlsr1
				raw_input('Press any key to start fitting...')
				f.canvas.mpl_disconnect(cid)
				vlsr2 = 0.0
			elif mode == 'double':
				cid = f.canvas.mpl_connect('button_press_event', onclick)
				raw_input('Click on the plot to select Vlsrs...')
				print clickvalue
				if len(clickvalue) >= 2:
					print 'Please select at least two velocities! The last two will be used.'
					vlsr1,vlsr2 = clickvalue[-2],clickvalue[-1]
				elif len(clickvalue) == 1:
					vlsr1 = clickvalue[-1]
					vlsr2 = 0.0
				elif len(clickvalue) == 0:
					vlsr1,vlsr2 = 0.0,0.0
				print 'Or input two velocities manually:'
				manualv = raw_input()
				manualv = manualv.split()
				if len(manualv) == 2:
					vlsr1,vlsr2 = np.float_(manualv)
				else:
					print 'Invalid input...'
				print 'The two Vlsrs are %0.2f km/s and %0.2f km/s.' % (vlsr1,vlsr2)
				raw_input('Press any key to start fitting...')
				f.canvas.mpl_disconnect(cid)
			else:
				vlsr1,vlsr2 = 0.0,0.0
		else:
			if mode == 'single':
				if spec_low.max() >= cutoff:
					print 'Reserved space...'
				else:
					vlsr1 = 0.0
				vlsr2 = 0.0
			elif mode == 'double':
				vlsr1,vlsr2 = 86.0,88.0
			else:
				vlsr1,vlsr2 = 0.0,0.0

		plt.text(0.02, 0.85, r'V$_{lsr}$=%.1f km/s' % vlsr1, transform=ax.transAxes, color='r', fontsize=15)
		fsky_init = ch3cn_info['frest'][0] * (1 - vlsr1 * 1e5 / c)

        # Add 4 parameters:
		params = Parameters()
		if vlsr1 != 0:
			params.add('Ntot', value=1e15, min=0, max=1e25)
			params.add('T', value=100, min=10)
			#params.add('sigma', value=0.0035, vary=False)
			params.add('sigma', value=0.0027, min=0, max=0.04)
			if varyf > 0:
				params.add('fsky', value=fsky_init, min=fsky_init-varyf*chanwidth, \
				max=fsky0_init+varyf*chanwidth)
			elif varyf == 0:
				params.add('fsky', value=fsky_init, vary=False)
		if vlsr2 != 0:
			print 'Reserved for two-component fitting.'
		
		# Run the non-linear minimization:
		if vlsr1 != 0 and vlsr2 != 0:
			result = minimize(__model_11_2c__, params, args=(faxis, spec))
		elif vlsr1 != 0 or vlsr2 != 0:
			result = minimize(__model_11__, params, args=(faxis, spec))
		else:
			unsatisfied = False
			continue

		final = spec + result.residual
		#report_fit(params)

		if interactive:
			plt.plot(faxis, final, 'r', label='Best-fitted model')
			if vlsr1 != 0 and vlsr2 != 0:
				print 'Reserved for two-component fitting.'
			elif vlsr1 != 0 or vlsr2 != 0:
				plt.text(0.02, 0.80, r'T$_{rot}$=%.1f($\pm$%.1f) K' % (result.params['T'].value,result.params['T'].stderr), transform=ax.transAxes, color='r', fontsize=15)
				plt.text(0.02, 0.75, r'N$_{tot}$=%.2e($\pm$%.2e) cm$^{-2}$' % (result.params['Ntot'].value,result.params['Ntot'].stderr), transform=ax.transAxes, color='r', fontsize=15)
				plt.text(0.02, 0.70, r'FWHM=%.2f($\pm$%.2f) km/s' % (c*result.params['sigma'].value/ch3cn_info['frest'][0]/1e5*2.355,c*result.params['sigma'].stderr/ch3cn_info['frest'][0]/1e5*2.355), transform=ax.transAxes, color='r', fontsize=15)
				plt.text(0.02, 0.65, r'$\tau_\mathrm{K=0}$=%.1e' % (__tau__(result.params['Ntot'].value,result.params['sigma'].value,result.params['T'].value,0)), transform=ax.transAxes, color='r', fontsize=15)
			plt.legend()
			plt.draw()
			print 'Is the fitting ok? y/n'
			yn = raw_input()
			if yn == 'y':
				unsatisfied = False
				currentT = time.strftime("%Y-%m-%d_%H:%M:%S")
				plt.savefig('CH3CN_fitting_'+currentT+'.png')
			else:
				unsatisfied = True
			#raw_input('Press any key to continue...')
			f.clear()
		else:
			unsatisfied = False

###############################################################################

ch3cn_info = __ch3cn_init__()

# Read the ASCII file.
# The frequency axis is assumed to be in GHz, and the y axis is T_B in K.
#spec, faxis = __readascii__('ch3cn_sgrb2.txt')
spec, faxis = __readascii__('dustridge.txt')
sourcename = 'Dust ridge core X'
chanwidth = abs(faxis[0] - faxis[-1]) / len(faxis)
print 'Channel width is %.4f GHz' % chanwidth
print 'Channel number is %d' % len(faxis)

# For SgrB2: flag negtiva channels
#spec[np.where(spec<0)] = 0
#spec = spec[50:390]
#faxis = faxis[50:390]

# This file is in Jy/beam so we convert it first:
#spec = 1.224e6 * spec / (ch3cn_info['frest'][0])**2 / (3.3 * 3.2)

# Reverse the axis:
#spec = spec[::-1]
#faxis = faxis[::-1]

spec +=0.1

# Run the fitting:
fit_spec(spec, faxis, Jupp=12, K_fit=8, cutoff=0.5, varyf=0, interactive=True, mode='single')
# Use a VLSR of 38.2 km/s in the fitting.

elapsed = (time.clock() - start)
print 'Stop the timer...'
print 'Time used: %0.0f seconds, or %0.1f minutes.' % (elapsed, elapsed/60.)

