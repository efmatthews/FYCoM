#Eric Matthews
#October 4, 2019
#Fission Yield Covariance Matrix Generation (FYCoM)
#Fitting P(nu,A) distributions to the England and Rider Evaluation

#Import Statements
#---------------------------------------------------------------------------------------------
import numpy
import scipy
from math import erfc
from math import erf
from math import sqrt
from math import exp
from math import pi
from math import isnan
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'figure.autolayout': True})
from decimal import Decimal
import warnings
warnings.filterwarnings("ignore")
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

def chi2A(x,ACUR,ACN,YIELDS,return_yields=False):
    """
    Function to calculate chi2 for P_nu_A distribution against ER yields
    Single A chain
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
    
def chi2A_fit(ACUR,ACN,YIELDS):
    """
    Function to perform differential evolution of data to chi2A
    """
    guess = [ (0.01,5.0), (0.01,3.0) ]
    fit_P_nu_A = differential_evolution( chi2A, guess, args=(ACUR,ACN,YIELDS) )
    fit_P_nu_A = fit_P_nu_A.x
    chi2 = chi2A(fit_P_nu_A,ACUR,ACN,YIELDS)
    
    return fit_P_nu_A, chi2

def chi2A_mod(P_nu_A,ACUR,ACN,YIELDS,return_yields=False):
    """
    Function to calculate chi2 for P_nu_A distribution against ER yields
    Single A chain
    """
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
                        yields_new[ZCN-key[0],ACN-key[1]-j] += P_nu_A[key[1]][j] * YIELDS[key]
                    except KeyError as e:
                        yields_new[ZCN-key[0],ACN-key[1]-j] = P_nu_A[key[1]][j] * YIELDS[key]

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
	Zp = ZAs[p][0]
	Ap = ZAs[p][1]

	#Set random seed for reproducability
	#-----------------------------------------------------------------------------------------
	numpy.random.seed(0)
	#-----------------------------------------------------------------------------------------


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
	        try:
	            yields[Z,A] += Y
	            yields_unc[Z,A] += Y_unc**2.0
	        except KeyError as e:
	            yields[Z,A] = Y
	            yields_unc[Z,A] = Y_unc**2.0
	YIELDS = yields
	for key in yields_unc:
	    yields_unc[key] = sqrt( yields_unc[key] )
	YIELDS_UNC = yields_unc
	ZCN = Zp
	ACN = Ap
	AMIN = A_min
	AMID = int( (A_max - A_min)/2.0 ) + A_min
	AMAX = A_max
	#-----------------------------------------------------------------------------------------


	#Perform differential evolution to fit the values for each A value 
	#-----------------------------------------------------------------------------------------
	fits = {}
	chi2s_before = []
	flags = []
	fit_last = None
	guess = [ (0.01,5.0), (0.01,3.0) ]

	for ACUR in range(AMIN,AMAX+1):
	    #Fit values
	    fit_P_nu_A, chi2 = chi2A_fit(ACUR,ACN,YIELDS)
	    flag = 0
	    chi2A(fit_P_nu_A,ACUR,ACN,YIELDS)

	    #If fit railed on any of the edges, use last A's fit
	    if( (fit_P_nu_A[0] == guess[0][0]) or (fit_P_nu_A[0] == guess[0][1]) or (fit_P_nu_A[1] == guess[1][0]) or (fit_P_nu_A[1] == guess[1][1]) or (chi2 == 0.0) ):
	        if( type(fit_last) != type(None) ):
	            fit_P_nu_A = fit_last
	            flag = 1
	            chi2 = chi2A(fit_P_nu_A,ACUR,ACN,YIELDS)
	    
	    chi2s_before.append( chi2 )
	    flags.append( flag )
	    fits[ACUR] = fit_P_nu_A
	    fit_last = fit_P_nu_A
	#-----------------------------------------------------------------------------------------


	#Generate P(v,A) table
	#-----------------------------------------------------------------------------------------
	P_nu_A = {}
	for A in range(AMIN,AMAX+1):
	    P_nu = []
	    for j in range(0,10):
	        mu, sigma = fits[A]
	        val = round( gauss_trunc_int(j,mu,sigma), 5 )
	        P_nu.append( val )
	    
	    norm = 1.0/sum(P_nu)
	    for j in range(0,10):
	        P_nu[j] = P_nu[j] * norm
	    
	    P_nu_A[A] = P_nu
	#-----------------------------------------------------------------------------------------


	#Output P_nu distribution to file
	#-----------------------------------------------------------------------------------------
	file = open( 'yields/P_nu_A/' + system + '_nu_data.csv', 'w' )
	for A in range(AMIN,AMAX+1):
	    file.write( str(A) )
	    for j in range(0,10):
	        file.write( ', ' + str( P_nu_A[A][j] ) )
	    file.write( '\n' )
	file.close()

	file = open( 'yields/P_nu_A/nu_reports/' + system + '_fit.csv', 'w' )
	file.write( 'A, mu, sigma, chi2, flag\n' )
	for A in range(AMIN,AMAX+1):
	    file.write( str(A) )
	    for j in range(0,2):
	        file.write( ', ' + str( fits[A][j] ) )
	    file.write( ', ' + str(chi2s_before[A-AMIN]) )
	    file.write( ', ' + str(flags[A-AMIN]) )
	    file.write('\n')
	file.close()
	#-----------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------------