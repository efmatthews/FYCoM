#Eric Matthews
#October 4, 2019
#Fission Yield Covariance Matrix Generation (FYCoM)
#Fitting P(nu,A) distributions to the England and Rider Evaluation

#Import Statements
#---------------------------------------------------------------------------------------------
import numpy
numpy.random.seed(0)
import scipy
from math import erfc
from math import erf
from math import sqrt
from math import exp
from math import pi
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from decimal import Decimal
#---------------------------------------------------------------------------------------------



#Global Variables
#---------------------------------------------------------------------------------------------
AMIN = 0
AMID = 0
AMAX = 0
ZCN = 0 
ACN = 0
YIELDS = {}
#---------------------------------------------------------------------------------------------



#Functions
#---------------------------------------------------------------------------------------------
def gauss_trunc(x,mu,sigma):
	"""
	Truncated Gaussian function
	x = point along distribution
	mu = centroid of distribution
	sigma = width of distribution 
	"""
	#norm = 1.0 / ( 1.0 - ( erfc( mu * sqrt( 1.0 / sigma**2.0 ) / sqrt(2.0) ) / ( 2.0 * mu * sqrt( 1.0 / sigma**2.0 ) ) ) )
	norm = 1.0 / ( 0.5 + 0.5 * erf( mu / (sqrt(2.0)*sigma) ) )
	return norm * ( 1.0 / (sigma * sqrt(2.0*pi)) ) * numpy.exp( -0.5 * ((x-mu)/sigma)**2.0 )

def gauss_trunc_integ(b,c,mu,sigma):
	"""
	Definite integral of truncated Gaussian function
	x = point along distribution
	mu = centroid of distribution
	sigma = width of distribution 
	"""
	norm = 1.0 / ( 0.5 + 0.5 * erf( mu / (sqrt(2.0)*sigma) ) )
	K = norm * ( 1.0 / (sigma * sqrt(2.0*pi)) )
	return K * sqrt(pi/2.0) * sigma * ( scipy.special.erf( (mu-b)/(sqrt(2.0)*sigma) ) - scipy.special.erf( (mu-c)/(sqrt(2.0)*sigma) ) ) 

def gauss_trunc_int(nu,mu,sigma,mid=0.5):
	"""
	Integer value of truncated Gaussian
	nu = integer point of distribution
	mu = centroid of distribution
	sigma = width of distribution 
	mid = midway point to integrate over (nu-(1-mid),nu+mid)
	"""
	nu = int(nu)

	lower = nu-(1.0-mid)
	if(lower < 0.0):
		lower = 0.0
	upper = nu+mid

	return gauss_trunc_integ(lower,upper,mu,sigma)

def gauss_trunc_bar(mu,sigma,h=0.001,upper=20):
	"""
	Get the expectation value of a truncated Gaussian
	mu = centroid of distribution
	sigma = width of distribution 
	h = integration constant
	upper = max value to integrate up to
	"""
	val = 0.0
	for x in numpy.arange(0,upper,h):
		val += gauss_trunc_integ(x,x+h,mu,sigma) * (x+(h/2.0))

	return val

def chi2A(x,return_yields=False):
	"""
	Function to calculate chi2 for P_nu_A distribution against ER yields
	"""
	mu = x[0]
	sigma = x[1]
	nu_max = 10
	A_upper = ACN - ACUR
	A_lower = ACN - ACUR - nu_max

	#Calculate yields in A chain
	yields_new = {}
	for key in YIELDS.keys():
		if( (key[1] >= A_lower) and (key[1] <= A_upper) ):
			for j in range(0,nu_max):
				if( (ACN-key[1]-j) == ACUR ):
					try:
						yields_new[ZCN-key[0],ACN-key[1]-j] += gauss_trunc_int(j,mu,sigma) * YIELDS[key]
					except KeyError as e:
						yields_new[ZCN-key[0],ACN-key[1]-j] = gauss_trunc_int(j,mu,sigma) * YIELDS[key]

	#Calculate chi2
	chi2 = 0.0
	for key in yields_new.keys():
		try:
			chi2 += ( YIELDS[key] - yields_new[key] )**2.0 / YIELDS[key]
		except KeyError as e:
			pass

	if( return_yields ):
		return chi2, yields_new
	else:
		return chi2
#---------------------------------------------------------------------------------------------



#Import the list of fissioning systems to evaluate covariance matrices for
#---------------------------------------------------------------------------------------------
file = open('yields/systems.txt','r')
lines = file.readlines()
file.close()
systems = []
ZAs = []
for line in lines:
	systems.append( line.split(',')[0].strip() )
	Z = int( line.split(',')[1] )
	A = int( line.split(',')[2] )
	ZAs.append( (Z,A) )
#---------------------------------------------------------------------------------------------



#Iterate over systems, fit P(nu,A) for each
#---------------------------------------------------------------------------------------------
for p in range(0,len(systems)):
	system = systems[p]
	#Read in the yields for the system
	#-----------------------------------------------------------------------------------------
	file = open( 'yields/' + system + '.csv', 'r' )
	lines = file.readlines()
	file.close()
	A_min = 300
	A_max = 0
	yields = {}
	yields_unc = {}
	for line in lines:
		parts = line.split(',')
		Z = int( parts[0] )
		A = int( parts[1] )
		I = int( parts[2] )
		Y = float( parts[3] )
		Y_unc = float( parts[4] )
		if( A < A_min ):
			A_min = A
		if( A > A_max ):
			A_max = A
		if( I == 0 ):
			yields[Z,A] = Y
			yields_unc[Z,A] = Y_unc
	YIELDS = yields
	YIELDS_UNC = yields_unc
	ZCN = ZAs[p][0]
	ACN = ZAs[p][1]
	AMIN = A_min
	AMID = int( (A_max - A_min)/2.0 ) + A_min
	AMAX = A_max
	#-----------------------------------------------------------------------------------------


	#Perform differential evolution to fit the values for each A value and save to file
	#-----------------------------------------------------------------------------------------
	file = open( 'yields/nu_reports/' + system + '_nu_report.csv', 'w' )
	fit_last = [ None, None ]
	nu_bars = []
	fits = []
	for ACUR in range(AMIN,AMAX+1):
		#Fit values
		#-------------------------------------------------------------------------------------
		guess = [ (0.01,5.0), (0.01,3.0) ]
		fit_P_nu_A = differential_evolution( chi2A, guess )
		fit_P_nu_A = fit_P_nu_A.x
		flag = 0
		chi2 = chi2A(fit_P_nu_A)


		#If fit railed on any of the edges, use last A's fit
		if( (fit_P_nu_A[0] == guess[0][0]) or (fit_P_nu_A[0] == guess[0][1]) or (fit_P_nu_A[1] == guess[1][0]) or (fit_P_nu_A[1] == guess[1][1]) or (chi2 == 0.0) ):
			if( fit_last[0] != None ):
				fit_P_nu_A = fit_last
				flag = 1
				chi2 = chi2A(fit_P_nu_A)

		nu_bar = gauss_trunc_bar(*fit_P_nu_A)
		nu_bars.append( nu_bar )
		fits.append( fit_P_nu_A )
		fit_last = fit_P_nu_A
		#-------------------------------------------------------------------------------------



		#Save results to file
		#-------------------------------------------------------------------------------------
		file.write( str(ACUR) + ', ' + str(fit_P_nu_A[0]) + ', ' + str(fit_P_nu_A[1]) + ', ' + str(chi2) + ', ' + str(nu_bar) + ', ' + str(flag) + '\n' )
		#-------------------------------------------------------------------------------------
	file.close()
	#-----------------------------------------------------------------------------------------



	#Plot nu_bar(A)
	#-----------------------------------------------------------------------------------------
	plt.plot( range(AMIN,AMAX+1), nu_bars, 'bo' )
	plt.xlabel( 'Mass Number (A)' )
	plt.ylabel( r'$\bar{\nu}$' )
	plt.savefig( 'yields/nu_reports/' + system + '_nu_bar.png', dpi=500 )
	plt.clf()
	#-----------------------------------------------------------------------------------------



	#Output P_nu distribution to file
	#-----------------------------------------------------------------------------------------
	file = open( 'yields/' + system + '_nu.csv', 'w' )
	for ACUR in range(AMIN,AMAX+1):
		file.write( str(ACUR) + ',' + str( gauss_trunc_int(0,fits[ACUR-AMIN][0],fits[ACUR-AMIN][1]) ) )
		for j in range(1,10):
			file.write( ', ' + str( gauss_trunc_int(j,fits[ACUR-AMIN][0],fits[ACUR-AMIN][1]) ) )
		file.write('\n')
	file.close()
	#-----------------------------------------------------------------------------------------



	#Plot yield residuals for one fissioning system
	#-----------------------------------------------------------------------------------------
	if( (system == 'U235F') and False ):
		for ACUR in range(AMIN,AMAX+1):
			chi2, yields_new = chi2A( fits[ACUR-AMIN], return_yields=True )

			Zs = []
			yields_new_vals = []
			yields_old_vals = []
			yields_old_vals_unc = []
			stdevs = []
			for key in yields_new.keys():
				Zs.append( key[0] )
				yields_new_vals.append( yields_new[key] )
				try:
					yields_old_vals.append( YIELDS[key] )
					yields_old_vals_unc.append( YIELDS_UNC[key] )
					stdev = abs((YIELDS[key] - yields_new[key])/YIELDS_UNC[key])
					stdevs.append( stdev )
				except KeyError as e:
					yields_old_vals.append( 0.0 )
					yields_old_vals_unc.append( 0.0 )
					stdevs.append( float('nan') )

			fig = plt.figure()
			frame1 = fig.add_axes((.1,.3,.8,.6))
			plt.errorbar( Zs, yields_old_vals, yerr=yields_old_vals_unc, fmt='bo', label='Evaluation' )
			plt.errorbar( Zs, yields_new_vals, fmt='r.', label=r'From $P(\nu,A)$' )
			plt.ylabel( 'Yield (%)' )
			plt.title( 'A = ' + str(ACUR) + ', ' + r'$\bar{\nu}$' + ' = ' + str(round(gauss_trunc_bar(*fits[ACUR-AMIN]),2)) + ', ' + r'$\chi^2$' + ' = ' + ('%.2E' % Decimal(str(chi2))) )
			plt.xlim( (min(Zs)-0.5,max(Zs)+0.5) )
			plt.yscale( 'log' )
			plt.legend()
			plt.tight_layout()
			frame2 = fig.add_axes((.1,.1,.8,.2))
			plt.plot( Zs, stdevs, 'k.-' )
			plt.xlabel( 'Atomic Number (Z)' )
			plt.ylabel( r'Deviation ($\sigma$)' )
			plt.xlim( (min(Zs)-0.5,max(Zs)+0.5) )
			plt.yscale( 'log' )
			plt.tight_layout()
			plt.savefig( 'yields/nu_reports/residuals/' + system + '_' + str(ACUR) + '_nu_bar.png', dpi=500 )
			plt.clf()
	#-----------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------------