#Eric Matthews
#September 23, 2019
#Fission Yield Covariance Matrix Generation (FYCoM)
#Independent covariances to cumulative covariances



#Constants
#---------------------------------------------------------------------------------------------
TRIALS = 1000
energies = { 'T':'thermal', 'F':'fission', 'H':'DT', 'D':'DD' }
elements = { 90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 98:'Cf', 99:'Es', 100:'Fm' }
#---------------------------------------------------------------------------------------------



#Import Statements
#---------------------------------------------------------------------------------------------
from math import sqrt
import numpy
import matplotlib.pyplot as plt
from pylab import rcParams
import warnings
warnings.filterwarnings("ignore")
import os
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



#Import a list of radioactive decays and branching ratios
#---------------------------------------------------------------------------------------------
file = open( 'FIER/FIER/input_data/decays.csv', 'r' )
lines = file.readlines()
file.close()

decays = {}
for line in lines: 
    Zp = int( line.split(',')[0] )
    Ap = int( line.split(',')[1] )
    Ip = int( line.split(',')[2] )
    BR = float( line.split(',')[5] ) / 100.0
    BR_unc = float( line.split(',')[6] ) / 100.0
    Zd = int( line.split(',')[7] )
    Ad = int( line.split(',')[8] )
    Id = int( line.split(',')[9] )
    try:
        decays[Zp,Ap,Ip][Zd,Ad,Id] = (BR, BR_unc)
    except KeyError as e:
        decays[Zp,Ap,Ip] = {}
        decays[Zp,Ap,Ip][Zd,Ad,Id] = (BR, BR_unc)
#---------------------------------------------------------------------------------------------



#Iterate over systems, covert covariance matrix for each
#---------------------------------------------------------------------------------------------
for p in range(0,len(systems)):
	system = systems[p]
	Zp = ZAs[p][0]
	Ap = ZAs[p][1]


	#Run FIER to generate decay stems for the system
	#-----------------------------------------------------------------------------------------
	file = open( 'FIER/FIER/deck.txt', 'w' )
	file.write( 'MODE:SINGLE\n' )
	file.write( 'ON DECAY PREDICTION\n' )
	file.write( 'input_data/isotopes.csv     ISOTOPES FILE\n' )
	file.write( 'input_data/decays.csv       DECAYS   FILE\n' )
	file.write( 'input_data/gammas.csv       GAMMAS   FILE\n' )
	file.write( 'YIELDS:FILE\n' )
	file.write( 'yields/' + elements[Zp] + ',' + str(Ap) + ',' + energies[system[-1]] + '   YIELDS   FILE\n' )
	file.write( 'NONE  CHAINS OUTPUT\n' )
	file.write( 'output/decay_stems.csv   STEMS OUTPUT\n' )
	file.write( 'NONE  POPS OUTPUT\n' )
	file.write( 'NONE   GAMMAS OUTPUT\n' )
	file.write( 'NONE  ERROR LOG\n' )
	file.write( 'INITIALIZE\n' )
	file.write( 'IRRADIATION\n' )
	file.write( '1.0,1.0\n' )
	file.write( 'POPULATIONS\n' )
	file.write( 'COUNTS\n' )
	file.write( 'END\n' )
	file.close()
	os.system( 'cd FIER/FIER && make > /dev/null' )
	#-----------------------------------------------------------------------------------------


	#Read in the decay stems
	#-----------------------------------------------------------------------------------------
	file = open( 'FIER/FIER/output/decay_stems.csv', 'r' )
	lines = file.readlines()
	file.close()

	stems = []
	sublines = []
	for i in range(0,len(lines)):
	    line = lines[i]
	    if( not('------' in line) ):
	        sublines.append( line )
	    else:
	        sublines = sublines[1:]
	        subsublines = [ [] ]
	        cur = 0
	        for j in range(0,len(sublines)):
	            if( not('---' in sublines[j]) ):
	                subsublines[cur].append( sublines[j] )
	            else:
	                cur += 1
	                subsublines.append( [] )
	        
	        subsublines = subsublines[:-1]
	        for sub in subsublines:
	            stem = []
	            for subsub in sub:
	                Z = int( subsub.split(',')[0] )
	                A = int( subsub.split(',')[1] )
	                I = int( subsub.split(',')[2] )
	                stem.append( (Z,A,I) )
	            stems.append( stem )
	        
	        sublines = []
	#-----------------------------------------------------------------------------------------


	#Calculate the transmutation probability from one species to another
	#-----------------------------------------------------------------------------------------
	trans_probs = {}
	for stem in stems:
	    key = (stem[0][0],stem[0][1],stem[0][2],stem[len(stem)-1][0],stem[len(stem)-1][1],stem[len(stem)-1][2])
	    prob = 1.0
	    prob_unc = 0.0
	    for i in range(1,len(stem)):
	        decay = decays[stem[i-1][0],stem[i-1][1],stem[i-1][2]][stem[i][0],stem[i][1],stem[i][2]]
	        prob = prob * decay[0] 
	        prob_unc += ( decay[1] / decay[0] )**2.0
	    prob_unc = prob * sqrt( prob_unc )
	    trans_probs[key] = ( prob, prob_unc )
	#-----------------------------------------------------------------------------------------


	#Import the yields for the selected system
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
	    yields[Z,A,I] = Y
	    yields_unc[Z,A,I] = Y_unc
	#-----------------------------------------------------------------------------------------


	#Import the covariance matrix for the independent yields
	#-----------------------------------------------------------------------------------------
	file = open( 'matrices/' + system + '_corr.csv', 'r' )
	lines = file.readlines()
	file.close()

	inds = []
	Y = []
	parts = lines[0].split(',')[1:]
	for part in parts:
	    val = int(part)
	    Z = int(val/10000)
	    I = int( (val - Z*10000)/1000 )
	    A = val - Z*10000 - I*1000
	    inds.append( (Z,A,I) )
	    Y.append( yields[Z,A,I] )

	n = len( inds )
	Y_cov = numpy.zeros( (n,n) )

	for i in range(1,len(lines)):
	    parts = lines[i].split(',')[1:]
	    for j in range(0,len(parts)):
	        Y_cov[i-1,j] = float( parts[j] )
	#-----------------------------------------------------------------------------------------


	#Calculate the cumulative yields n = TRIALS times using resampled independent yields
	#-----------------------------------------------------------------------------------------
	numpy.random.seed(0)
	trials_res = numpy.zeros( (len(inds),TRIALS) )
	for n in range(0,TRIALS):
	    yields_vard = numpy.random.multivariate_normal( Y, Y_cov )
	    yields_vard_dict = {}
	    for i in range(0,len(inds)):
	        yields_vard_dict[inds[i]] = yields_vard[i]
	    
	    yields_cml = {}
	    for key in trans_probs:
	        prob = numpy.random.normal( trans_probs[key][0], trans_probs[key][1] )
	        try:
	            yields_cml[key[3],key[4],key[5]] += prob * yields_vard_dict[key[3],key[4],key[5]]
	        except KeyError as e:
	            try:
	                yields_cml[key[3],key[4],key[5]] = prob * yields_vard_dict[key[3],key[4],key[5]]
	            except KeyError as e:
	                pass
	    
	    for i in range(0,len(inds)):
	        trials_res[i,n] = yields_cml[inds[i]]
	#-----------------------------------------------------------------------------------------


	#Calculate the correlation and covariance matrix
	#-----------------------------------------------------------------------------------------
	yields_corr = numpy.corrcoef( trials_res )
	yields_cov = numpy.cov( trials_res )
	#-----------------------------------------------------------------------------------------


	#Import the England and Rider cumulative yields and generate a covariance matrix based 
	#on the evaluated cumulative yield uncertainties and the above correlation matrix
	#-----------------------------------------------------------------------------------------
	file = open( 'yields/cumulative/' + system + '_cml.csv', 'r' )
	lines = file.readlines()
	file.close()

	yields_cml = {}
	yields_cml_unc = {}
	for line in lines:
	    parts = line.split(',')
	    Z = int( parts[0] )
	    A = int( parts[1] )
	    I = int( parts[2] )
	    Y = float( parts[3] )
	    Y_unc = float( parts[4] )
	    yields_cml[Z,A,I] = Y
	    yields_cml_unc[Z,A,I] = Y_unc

	df_std = []
	df_stdT = []
	for ind in inds:
	    try:
	        df_std.append( yields_cml_unc[ind] )
	        df_stdT.append( [yields_cml_unc[ind]] )
	    except KeyError as e:
	        df_std.append( 0.0 )
	        df_stdT.append( [0.0] )
	df_std = numpy.array( df_std )
	df_stdT = numpy.array( df_stdT )
	yields_cov_normed = df_std * yields_corr * df_stdT
	#-----------------------------------------------------------------------------------------


	#Plot the correlation matrix
	#-----------------------------------------------------------------------------------------
	rcParams['figure.figsize'] = 14, 14
	rcParams['font.size'] = 22
	plt.matshow(yields_corr,cmap="RdBu",vmin=-1.0,vmax=1.0)
	cbar = plt.colorbar()
	cbar.set_label('Correlation Coeff.')
	plt.tight_layout()
	plt.xlabel('FY Index')
	plt.ylabel('FY Index')
	plt.savefig( 'figures/' + system + '_cml_corr.eps', format='eps', dpi=300 )
	plt.savefig( 'figures/' + system + '_cml_corr.png', dpi=500 )
	plt.clf()
	Image( 'figures/' + system + '_cml_corr.png' )
	#-----------------------------------------------------------------------------------------


	#Save the correlation and covariance matrices to file
	#-----------------------------------------------------------------------------------------
	#Covariance save
	file = open( 'matrices/' + system + '_cml_cov.csv', 'w' )
	for key in inds:
		file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
	file.write('\n')
	for i in range(0,len(inds)):
		key = inds[i]
		file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
		for j in range(0,len(inds)):
			key2 = inds[j]
			file.write( ', ' + str(yields_cov[i,j]) )
		file.write( '\n' )
	file.close()

	#Normed covariance save
	file = open( 'matrices/' + system + '_cml_normed_cov.csv', 'w' )
	for key in inds:
		file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
	file.write('\n')
	for i in range(0,len(inds)):
		key = inds[i]
		file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
		for j in range(0,len(inds)):
			key2 = inds[j]
			file.write( ', ' + str(yields_cov_normed[i,j]) )
		file.write( '\n' )
	file.close()

	#Correlation save
	file = open( 'matrices/' + system + '_cml_corr.csv', 'w' )
	for key in inds:
		file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
	file.write('\n')
	for i in range(0,len(inds)):
		key = inds[i]
		file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
		for j in range(0,len(inds)):
			key2 = inds[j]
			file.write( ', ' + str(yields_corr[i,j]) )
		file.write( '\n' )
	file.close()
	#-----------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------