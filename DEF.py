# -*- coding: utf-8 -*-
"""
Debye-Einstein-Fit:
Fitting temperature dependent pyroelectric coefficients with a Debye and up to
six Einstein functions
Author: Sven Jachalke
Date: 2017-11-04
"""

from pylab import *
from tubafcdpy import *
import pandas as pd
import lmfit as fit
from mpmath import polylog
from numpy import log10
import warnings
warnings.filterwarnings("ignore")

# Debye-Einstein-functions ----------------------------------------------------
# non-integral Debye function (Dubinov et al, 2008)
def debye_non_integral(x):
	'''
	Non-Integral function for the Debye Integral
	(has only to be multiplied with factor A!)
	x = T_D/T (Debye temperature and absolute temperature)
	'''
	
	# individual summands of integral replacement
	s1 = 4./5. * pi**4 / x**3
	s2 = 3*x*exp(-x) / (exp(-x)-1)
	s3 = 12*log(1-exp(-x))
	s4 = -36/x * polylog(2, exp(-x))
	s5 = -72/x**2 * polylog(3, exp(-x))
	s6 = -72/x**3 * polylog(4, exp(-x))
	
	DC3 = (s1+s2+s3+s4+s5+s6)
						
	return DC3	

def Debye(x):
	'''
	vectorized function which takes a nested sequence of objects or numpy 
	arrays as inputs and returns an single or tuple of numpy array as output.
	'''
	f =  vectorize(debye_non_integral, otypes=[complex])
	return f(x).real	

# Einstein function
def Einstein(T,T_E):
	'''
	Einstein function
	T ... array of temperature range
	T_E ... Einstein temperature
	'''
	return (T_E/T)**2 * ((exp(T_E/T))/(exp(T_E/T)-1)**2)

# Combination of Debye and several Einstein functions with factors
def Debye_Einstein_Modell(params, T, data=None, n_Einstein=1):
	'''
	Combination of one Debye and several Einstein functions
	parameteters:
		params - lmfit ParamsDict with factors and Debye/Einstein temperatures
		T - numpy array with temperatures
		data - numpy array with pyroelectric coefficients
		n_Einstein - number of needed Einstein contributions 
	'''
	
	# DebyePart
	T_D = params['T_D'].value
	A = params['A'].value
	
	DEM = A*Debye(T_D/T)
	
	# add EinsteinPart if needed
	if n_Einstein != 0:
		#EinsteinPart
		for i in range(n_Einstein):
			factor = 'B'+str(i+1)
			temperature = 'T_E'+str(i+1)
			
			B = params[factor].value
			T_E = params[temperature].value
			
			DEM = DEM + B*Einstein(T,T_E)
		
	if data is None:
		return DEM
	else:
		return DEM-data

# Conversion functions --------------------------------------------------------
# Temperature to cm-1 conversion
def Temp2cm1(temperature):
	"""
	converts temperature (K) to wave number (cm^-1)
	"""
	return temperature*0.695

# cm-1 to Temperature conversion
def cm2Temp(cm):
	"""
	converters wave number (cm^-1) to temperature (K)
	"""
	return cm/0.695

# Temperature to eV conversion
def Temp2eV(temperature):
	"""
	converts temperature (K) to energy (eV)
	"""
	return temperature*8.621738e-5	

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
style.use('science')
fit_flag = True 				#True - perform fit, False - use defined values
name = 'Debye-Einstein-Fit' 	#adjust name for plots and log files

T_start = 0					#lower limit of temperature range
T_end = 600					#upper limit of temperature range
T_step = 150					#temperature steps between low and high

# start parameters
# column1 = Debye Parameters (factor and Debye temperature)
# column2 to column 7 = Einstein Parameters (factor and Einstein temperature)
# vary_factor: True = factor will be fitted, False = factor remains untouched
# vary_temperatures: True = temperatur will be fitted, False = temperature remains untouched

#Lines:1977 params (only on low T data!) -- use of 2x Einstein without Debye
#factors = 		[	1e-10,		2.556e-5,		1.4337e-4,	0,			0,			0,		0]
#temperatures = 	[	0.1,			116.54,		313.669,		0,			0,			0,		0]
#vary_factors =	[	False,		True,		True,		False,		False,		False,	False]
#vary_temperatures =[	False,		False,		False,		False,		False,		False,	False]
#redX2: 5.061385e-11

##1x Einstein 
factors = 		[	1.4755e-4,	7.42261e-4,	0,			0,			0,			0,		0]
temperatures = 	[	339.061,		2030,71,		0,			0,			0,			0,		0]
vary_factors =	[	True,		True,		False,		False,		False,		False,	False]
vary_temperatures =[	True,		True,		False,		False,		False,		False,	False]
##redX2: 4.481420e-12 --> T_D näher an CV Daten!

##1x Einstein mit ThomasDaten
#factors = 		[	1.563e-4,		8.90189e-4,	0,			0,			0,			0,		0]
#temperatures = 	[	358.969,		2189.93,		0,			0,			0,			0,		0]
#vary_factors =	[	True,		True,		False,		False,		False,		False,	False]
#vary_temperatures =[	True,		False,		False,		False,		False,		False,	False]
##redX2: 7.405409e-12 --> etwas schlechter

#1x Einstein mit 1522cm und mittlerer Debye-Temperatur aus Literatur (451K)
#factors = 		[	1.563e-4,		8.36e-4,		0,			0,			0,			0,		0]
#temperatures = 	[	451,			2189.93,		0,			0,			0,			0,		0]
#vary_factors =	[	True,		True,		True,		False,		False,		False,	False]
#vary_temperatures =[	False,		False,		False,		False,		False,		False,	False]
#redX2: 7.405409e-12 --> etwas schlechter

#2x Einstein 
#factors = 		[	3.996e-4,		7.396e-4,		1.071e-4,			0,			0,			0,		0]
#temperatures = 	[	214.18,		2026,33,		286.31,			0,			0,			0,		0]
#vary_factors =	[	True,		True,		True,		False,		False,		False,	False]
#vary_temperatures =[	True,		True,		True,		False,		False,		False,	False]
#redX2: 4.541127e-12

#2x Einstein mit ThomasDaten -- 1755cm-1 Bande sinnlos!
#factors = 		[	1.563e-4,		8.90189e-4,	5.6997e-10,			0,			0,			0,		0]
#temperatures = 	[	358.969,		2189.93,		2525.17,			0,			0,			0,		0]
#vary_factors =	[	True,		True,		True,		False,		False,		False,	False]
#vary_temperatures =[	False,		False,		False,		False,		False,		False,	False]
#redX2: 7.405409e-12 --> etwas schlechter

#2x Einstein mit ThomasDaten -- 600cm-1 Bande sinnlos!
#factors = 		[	1.563e-4,		1e-4,		8.90189e-4,			0,			0,			0,		0]
#temperatures = 	[	287.76,		863.3,		2189.93,			0,			0,			0,		0]
#vary_factors =	[	True,		True,		True,		False,		False,		False,	False]
#vary_temperatures =[	False,		False,		False,		False,		False,		False,	False]
#redX2: 7.405409e-12 --> etwas schlechter

T = linspace(T_start,T_end,T_step)

n_Einstein = sum(array(factors)!=0)-1	#Count number of Einsteins defined
n_Debye = 1
n = n_Einstein+n_Debye 				#calc total number of modes

# load data -------------------------------------------------------------------
# (has to be adjusted for individual files. This example uses data from the
# literature and own data. In the end a combined pandas DataFrame is created)

# Literature data
Lit_data = pd.read_csv("LT-LandB-77L2-Lines1977.txt",skiprows=1,names=['temp','p'],delimiter=',',decimal='.')
Lit_data['p'] = -abs(Lit_data.p)

# own data
Own_data = pd.read_csv("2017-06-27_15-30_LiTaO3-C-LT-F1_SineWave+LinRamp_PyroData.txt",skiprows=1,usecols=[1,2,8],names=['temp','p','perror'],delimiter='\t')
Own_data['p'] = -abs(Own_data.p) #+ 5e-6

# fit combined data set of Shaldin and own data
data = pd.concat([Lit_data,Own_data])
data = data.sort_values(by=['temp'])

# Initialize ParametersDicts with scaling factors and temperatures ------------
Sim = fit.Parameters()
Sim.add('A', value=factors[0], vary=vary_factors[0])
Sim.add('T_D', value=temperatures[0], vary=vary_temperatures[0],min=0)
if n_Einstein >= 1:
	Sim.add('B1', value=factors[1], vary=vary_factors[1],min=0)
	Sim.add('T_E1', value=temperatures[1], vary=vary_temperatures[1],min=0)
if n_Einstein >= 2:
	Sim.add('B2', value=factors[2], vary=vary_factors[2],min=0)
	Sim.add('T_E2', value=temperatures[2], vary=vary_temperatures[2],min=0)
if n_Einstein >= 3:
	Sim.add('B3', value=factors[3], vary=vary_factors[3],min=0)
	Sim.add('T_E3', value=temperatures[3], vary=vary_temperatures[3],min=0)
if n_Einstein >= 4:
	Sim.add('B4', value=factors[4], vary=vary_factors[4],min=0)
	Sim.add('T_E4', value=temperatures[4], vary=vary_temperatures[4],min=0)
if n_Einstein >= 5:
	Sim.add('B5', value=factors[5], vary=vary_factors[5],min=0)
	Sim.add('T_E5', value=temperatures[5], vary=vary_temperatures[5],min=0)
if n_Einstein >= 6:
	Sim.add('B6', value=factors[6], vary=vary_factors[6],min=0)
	Sim.add('T_E6', value=temperatures[6], vary=vary_temperatures[6],min=0)

# -----------------------------------------------------------------------------
# fitting
if fit_flag == True:
	result = fit.minimize(Debye_Einstein_Modell, Sim, args=(abs(data.temp), abs(data.p), n_Einstein))
	fitted_params = result.params
	print('red X2: %e'%result.redchi)
else:
	fitted_params = Sim

# plotting --------------------------------------------------------------------
f = figure('Shaldin',figsize=(8,6))
ax = f.add_subplot(111)

# experimental data
ax.plot(Lit_data.temp,Lit_data.p*1e6,color=tubafblue(),linestyle='',marker='o',label=r'\textsc{Lines} \textit{et al.}, Phys. Rev. Lett., Vol. 39(21), pp. 1362--1365, 1977')
#ax.plot(LN_Shaldin_air.temp,LN_Shaldin_air.p*1e6,color=tubafred(),linestyle='',marker='o',label='Shaldin2008 - air (1000 K, 7h)')
#ax.plot(LN_init.temp,LN_init.p*1e6,color=tubafgreen(),linestyle='',marker='o',label='own data ($+\SI{5}{\micro\coulomb\per\kelvin\per\square\meter}$)')
ax.plot(Own_data.temp,Own_data.p*1e6,color=tubafgreen(),linestyle='',marker='o',label='own data')

# plot contributions
for i in range(n):
	# calculate simulated p --------------------------------------------------
	p_simulation = Debye_Einstein_Modell(fitted_params,T,data=None,n_Einstein=i)-fitted_params['A'].value*Debye(fitted_params['T_D']/T)
	
	if i==0:
		p_simulation = Debye_Einstein_Modell(fitted_params,T,data=None,n_Einstein=i)
		B = ('%.2e'%fitted_params['A']).split('e-')
		ax.plot(T,-p_simulation*1e6,color=tubafred(),label=r'$T_{\mathrm{D}} = %.0f$ K (%.0f cm$^{-1}$), $A = %se^{-%s}$' % (fitted_params['T_D'].value,Temp2cm1(fitted_params['T_D'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['A'].value,fitted_params['T_D'].value,Temp2cm1(fitted_params['T_D'].value)))
	if i==1:
		B = ('%.2e'%fitted_params['B1']).split('e-')
		ax.plot(T,-p_simulation*1e6,color=tubaforange(),label=r'$T_{\mathrm{E1}} = %.0f$ K (%.0f cm$^{-1}$), $B_1 = %se^{-%s}$' % (fitted_params['T_E1'].value,Temp2cm1(fitted_params['T_E1'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B1'].value,fitted_params['T_E1'].value,Temp2cm1(fitted_params['T_E1'].value)))
	if i==2:
		B = ('%.2e'%fitted_params['B2']).split('e-')
		ax.plot(T,-p_simulation*1e6,color=tubafgreen(),label=r'$T_{\mathrm{E2}} = %.0f$ K (%.0f cm$^{-1}$), $B_2 = %se^{-%s}$' % (fitted_params['T_E2'].value,Temp2cm1(fitted_params['T_E2'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B2'].value,fitted_params['T_E2'].value,Temp2cm1(fitted_params['T_E2'].value)))
	if i==3:
		B = ('%.2e'%fitted_params['B3']).split('e-')
		ax.plot(T,-p_simulation*1e6,color=tubafcyan(),label=r'$T_{\mathrm{E3}} = %.0f$ K (%.0f cm$^{-1}$), $B_3 = %se^{-%s}$' % (fitted_params['T_E3'].value,Temp2cm1(fitted_params['T_E3'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B3'].value,fitted_params['T_E3'].value,Temp2cm1(fitted_params['T_E3'].value)))
	if i==4:
		B = ('%.2e'%fitted_params['B4']).split('e-')
		ax.plot(T,-p_simulation*1e6,color=tubafblue(),label=r'$T_{\mathrm{E4}} = %.0f$ K (%.0f cm$^{-1}$), $B_4 = %se^{-%s}$' % (fitted_params['T_E4'].value,Temp2cm1(fitted_params['T_E4'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B4'].value,fitted_params['T_E4'].value,Temp2cm1(fitted_params['T_E4'].value)))
	if i==5:
		B = ('%.2e'%fitted_params['B5']).split('e-')
		ax.plot(T,-p_simulation*1e6,label=r'$T_{\mathrm{E5}} = %.0f$ K (%.0f cm$^{-1}$), $B_5 = %se^{-%s}$' % (fitted_params['T_E5'].value,Temp2cm1(fitted_params['T_E5'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B5'].value,fitted_params['T_E5'].value,Temp2cm1(fitted_params['T_E5'].value)))
	if i==6:
		B = ('%.2e'%fitted_params['B6']).split('e-')
		ax.plot(T,-p_simulation*1e6,label=r'$T_{\mathrm{E6}} = %.0f$ K (%.0f cm$^{-1}$), $B_6 = %se^{-%s}$' % (fitted_params['T_E6'].value,Temp2cm1(fitted_params['T_E6'].value),B[0],B[1][1]) )
		print('%i: %.2e - %0.f K - %.2f cm-1'%(i,fitted_params['B6'].value,fitted_params['T_E6'].value,Temp2cm1(fitted_params['T_E6'].value)))
	
# plot complete fit
p_simulation = Debye_Einstein_Modell(fitted_params,T,data=None,n_Einstein=n_Einstein)
ax.plot(T,-p_simulation*1e6,label=r'complete fit',color='k')
	
ax.legend(loc=3)
ax.grid()
ax.set_xlabel(r'$T$ (K)')
ax.set_ylabel(r'$p$ (\si{\micro\coulomb\per\kelvin\per\square\meter})')
ax.set_ylim(-450,0)

f.tight_layout()
f.savefig(name+'.pdf')