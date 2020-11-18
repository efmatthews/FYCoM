#Eric Matthews
#September 23, 2019
#Fission Yield Covariance Matrix Generation (FYCoM)

#Constants
#---------------------------------------------------------------------------------------------
TRIALS = 10000
evaluations = [ 'ENDF.B.VII.1', 'JEFF.3.1' ]
NUM_PRECS = { 'ENDF.B.VII.1':0.0001, 'JEFF.3.1':0.0005 }
NUM_PRECS_EXCEPTIONS = { 'ENDF.B.VII.1:Am241H':0.0005, 'JEFF.3.1:U233H':0.0007 }
#---------------------------------------------------------------------------------------------

#Import Statements
#---------------------------------------------------------------------------------------------
import random
import numpy
import matplotlib.pyplot as plt
from pylab import rcParams
import warnings
warnings.filterwarnings("ignore")
import time
#---------------------------------------------------------------------------------------------



#Iterate over systems, generate a covariance matrix for each
#---------------------------------------------------------------------------------------------
for evaluation in evaluations:
	#Import the list of fissioning systems to evaluate covariance matrices for
	#-----------------------------------------------------------------------------------------
	file = open('yields/' + evaluation + '/systems.txt','r')
	lines = file.readlines()
	file.close()
	systems = []
	ZAs = []
	for line in lines:
		systems.append( line.split(',')[0].strip() )
		Z = int( line.split(',')[1] )
		A = int( line.split(',')[2] )
		ZAs.append( (Z,A) )
	#-----------------------------------------------------------------------------------------



	for p in range(0,len(systems)):
		time_start = time.time()
		system = systems[p]
		print( system )
		Zp = ZAs[p][0]
		Ap = ZAs[p][1]

		if( evaluation+':'+system in NUM_PRECS_EXCEPTIONS.keys() ):
			NUM_PREC = NUM_PRECS_EXCEPTIONS[evaluation+':'+system]
		else:
			NUM_PREC = NUM_PRECS[evaluation]

		#Set random seed for reproducability
		#-------------------------------------------------------------------------------------
		random.seed(0)
		numpy.random.seed(0)
		#-------------------------------------------------------------------------------------

		#Read in the yields for the system
		#-------------------------------------------------------------------------------------
		file = open( 'yields/' + evaluation + '/independent/' + system + '.csv', 'r' )
		lines = file.readlines()
		file.close()
		A_min = 300
		A_max = 0
		yields = {}
		yields_unc = {}
		Y_tot_orig = 0.0
		Z_tot_orig = 0.0
		A_tot_orig = 0.0
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
		    Y_tot_orig += Y
		    Z_tot_orig += Y * Z
		    A_tot_orig += Y * A
		    yields_unc[Z,A,I] = Y_unc
		key_list = list( yields.keys() )
		key_list_Z = numpy.asarray( [sublist for sublist in key_list] ).T[0,:]
		key_list_A = numpy.asarray( [sublist for sublist in key_list] ).T[1,:]
		key_list.sort()
		#-------------------------------------------------------------------------------------



		#Generate a table of isomer-to-ground state ratios as these should be maintained
		#-------------------------------------------------------------------------------------
		isomer_ratios = {}
		for key in yields.keys():
		    if( not( (key[0],key[1]) in isomer_ratios.keys() ) ):
		        ratios = [ 0.0, 0.0, 0.0 ]
		        for key2 in yields.keys():
		            if( (key2[0] == key[0]) and (key2[1] == key[1]) ):
		                ratios[ key2[2] ] = yields[key2]

		        norm = 1.0 / sum(ratios)
		        for i in range(0,2):
		            ratios[i] = norm * ratios[i]

		        isomer_ratios[key[0],key[1]] = ratios
		#-------------------------------------------------------------------------------------



		#Read in the P(nu,A) distribution for the system
		#-------------------------------------------------------------------------------------
		file = open( 'yields/' + evaluation + '/P_nu_A/' + system + '_nu_data.csv', 'r' )
		lines = file.readlines()
		file.close()
		P_nu_A = {}
		for line in lines:
		    parts = line.split(',')
		    A = int( parts[0] )
		    P_nu = []
		    for part in parts[1:]:
		        P_nu.append( float( part ) )
		    #Ensure normalization
		    #---------------------------------------------------------------------------------
		    norm = 1.0 / sum( P_nu )
		    for i in range(0,len(P_nu)):
		        P_nu[i] = norm * P_nu[i]
		    #---------------------------------------------------------------------------------
		    P_nu_A[A] = P_nu

		nu_bar = 0.0
		for key in key_list:
		    for i in range(0,10):
		        nu_bar += P_nu_A[key[1]][i] * i * yields[key]/200.0
		#-------------------------------------------------------------------------------------



		#Calculate middle A of system
		#-------------------------------------------------------------------------------------
		A_mid = None
		Y_tot = 0.0
		for A in range(A_min,A_max):
		    for Z in range(0,100):
		        for I in range(0,3):
		            try:
		                Y_tot += yields[Z,A,I]
		            except KeyError as e:
		                pass
		    if( Y_tot >= 100.0 ):
		        A_mid = A
		        break
		#-------------------------------------------------------------------------------------



		#Perform n = TRIALS of resamplings of fission yield library
		#-------------------------------------------------------------------------------------
		numpy.random.seed(0)
		random.seed(0)
		trials_res = numpy.zeros( (len(key_list),TRIALS) )
		extra_products = {}
		n = 0
		while(n < TRIALS):
		    yields_vard = {}

		    #Choose whether to resample from the heavy or light side
		    #---------------------------------------------------------------------------------
		    side = int( n % 2 )

		    if( side == 0 ): #light side to be resampled
		        As = list( range( A_min, A_mid ) )
		        As2 = list( range( A_mid, A_max+1 ) )
		    else: #heavy side to be resampled
		        As = list( range( A_mid, A_max+1 ) )
		        As2 = list( range( A_min, A_mid ) )
		    #---------------------------------------------------------------------------------



		    #Resample yields on side
		    #Go along each A in As, selected a Z in a weighted fashion
		    #Increase all yields along the A chain by the percent A,Z was increased
		    #---------------------------------------------------------------------------------
		    for A in As:
		        #Get keys with this A 
		        #-----------------------------------------------------------------------------
		        keys = []
		        for key in yields.keys():
		            if( key[1] == A ):
		                keys.append( key )
		        #-----------------------------------------------------------------------------



		        #Choose which yield will be resampled
		        #-----------------------------------------------------------------------------
		        weights = []
		        for key in keys:
		            weights.append( 1.0 / (yields_unc[key]/yields[key])**2.0 )
		        weights = numpy.array( weights ) / sum(weights)

		        selected = numpy.random.choice( list(range(0,len(keys))), p=weights )
		        selected = keys[ selected ]
		        #-----------------------------------------------------------------------------



		        #Sample selected yield and scale all yields in A chain up by the same percentage
		        #-----------------------------------------------------------------------------
		        Y_vard = numpy.random.normal( yields[selected], yields_unc[selected] )
		        while( Y_vard <= 0.0 ):
		            Y_vard = numpy.random.normal( yields[selected], yields_unc[selected] )

		        percent_shift = Y_vard / yields[selected]

		        for key in keys:
		            yields_vard[key] = percent_shift * yields[key]
		        #-----------------------------------------------------------------------------
		    #---------------------------------------------------------------------------------


		    #Normalize yields that were resampled from one side to sum to 1.0
		    #---------------------------------------------------------------------------------
		    norm = 0.0
		    for key in yields_vard.keys():
		        norm += yields_vard[key]
		    norm = 100.0 / norm

		    for key in yields_vard.keys():
		        yields_vard[key] = norm * yields_vard[key]
		    #---------------------------------------------------------------------------------


		    #Use the P(nu,A) distribution to calculate yields on the other side of the distribution
		    #---------------------------------------------------------------------------------
		    for key in list( yields_vard.keys() ):
		        Z_comp =  Zp - key[0]
		        P_nu = P_nu_A[key[1]][:]
		        nus = list( range(0,len(P_nu)) )
		        
		        #Check if each nu P_nu leads to an existing isotope in the evaluation
		        #If not remove that nu from P_nu and renormalize P_nu
		        #-----------------------------------------------------------------------------
		        for nu in range(len(P_nu)-1,-1,-1):
		            A_comp = Ap - nu - key[1]
		            if( not( (Z_comp,A_comp) in isomer_ratios.keys() ) ):
		                del nus[nu]
		                del P_nu[nu]

		        if( (len(nus) == 0) or (sum(P_nu) == 0) ):
		            continue

		        norm = 1.0 / sum(P_nu)
		        for i in range(0,len(P_nu)):
		            P_nu[i] = norm * P_nu[i]
		        #-----------------------------------------------------------------------------



		        #Calculate partial yields
		        #-----------------------------------------------------------------------------
		        for j in range(0,len(nus)):
		            nu = nus[j]
		            A_comp = Ap - nu - key[1]
		            for i in range(0,2):
		                Y = yields_vard[key] * P_nu[j] * isomer_ratios[Z_comp, A_comp][i]
		                if( Y > 0.0 ):
		                    try:
		                        yields_vard[Z_comp,A_comp,i] += Y
		                    except KeyError as e:
		                        yields_vard[Z_comp,A_comp,i] = Y
		        #-----------------------------------------------------------------------------
		    #---------------------------------------------------------------------------------
		    

		    #Renormalize so that all yields sum to 2.0
		    #---------------------------------------------------------------------------------
		    Y_tot = 0.0
		    for key in yields_vard.keys():
		        Y_tot += yields_vard[key]
		    
		    norm = 200.0 / Y_tot
		    for key in yields_vard.keys():
		        yields_vard[key] = norm * yields_vard[key]
		    #---------------------------------------------------------------------------------
		    
		    
		    #Some products are produced from the P(nu,A) distribution that aren't in the evaluation
		    #Move the yield given to these products to the product with the nearest A
		    #Note: these all are isomer yields with value 0.0 that appear to get assigned a value due to the 
		    #application of `isomer_ratios' 
		    #---------------------------------------------------------------------------------
		    toRemove = []
		    toAdd = {}
		    for key in yields_vard.keys():
		        if( not(key in key_list) and not(key in extra_products.keys()) ):
		            closest_key = None
		            diff = 100000
		            for key2 in key_list:
		                if( key[0] == key2[0] ):
		                    diff_cur = abs( key[1] - key2[1] )
		                    if( diff_cur < diff ):
		                        closest_key = key2
		                        diff = diff_cur
		            extra_products[key] = key2
		            
		        if( key in extra_products.keys() ):
		            try:
		                yields_vard[ extra_products[key] ] = yields_vard[ extra_products[key] ] + yields_vard[ key ]
		            except KeyError as e:
		                toAdd[ extra_products[key] ] = yields_vard[ key ]
		            toRemove.append( key )
		    
		    for key in toRemove:
		        del yields_vard[key]
		    for key in toAdd:
		        yields_vard[key] = toAdd[key]
		    #---------------------------------------------------------------------------------
		    
		    
		    #Check that this trial meets the numerical precision criteria, if not break and repeat
		    #---------------------------------------------------------------------------------
		    Y_tot = 0.0
		    for key in yields_vard.keys():
		        Y_tot += yields_vard[key]
		    Y_tot = Y_tot / 200.0
		    Z_tot = 0.0
		    A_tot = 0.0
		    for key in yields_vard.keys():
		        Z_tot += yields_vard[key] * key[0]
		        A_tot += yields_vard[key] * key[1]
		    Z_tot = abs( ((Z_tot/100.0) / (Zp*Y_tot)) - 1.0 )
		    A_tot = abs( ((A_tot/100.0) / ((Ap-nu_bar)*Y_tot)) - 1.0 )
		    if( ( Z_tot > NUM_PREC ) or ( A_tot > NUM_PREC ) ):
		        continue
		    #---------------------------------------------------------------------------------


		    #Place resampled yields into a numpy array for rapid covariance matrix calculation
		    #---------------------------------------------------------------------------------
		    for i in range(0,len(key_list)):
		        key = key_list[i]
		        try:
		            trials_res[i,n] = yields_vard[key]
		        except KeyError as e:
		            trials_res[i,n] = yields[key]
		    n += 1 #Increase index
		    print( '   ' + str(n) + '/' + str(TRIALS), end='\r')
		    #---------------------------------------------------------------------------------
		#-------------------------------------------------------------------------------------



		#Calculate the primary covariance and correlation matrices
		#-------------------------------------------------------------------------------------
		yields_corr = numpy.corrcoef( trials_res )
		yields_cov = numpy.cov( trials_res )
		#-------------------------------------------------------------------------------------



		#Calculate the normalized covariance and correlation matrices
		#-------------------------------------------------------------------------------------
		df_std = []
		df_stdT = []
		for key in key_list:
		    df_std.append( yields_unc[key] )
		    df_stdT.append( [yields_unc[key]] )
		df_std = numpy.array( df_std )
		df_stdT = numpy.array( df_stdT )
		yields_cov_prenorm = df_std * yields_corr * df_stdT

		total_sum = yields_cov_prenorm.sum()
		if( total_sum > 0.0 ):
		    non_diag = yields_cov_prenorm[:,:]
		    numpy.fill_diagonal( non_diag, 0.0 )
		    negative_sum = abs( non_diag[non_diag < 0.0].sum() )
		    norm = (total_sum + negative_sum) / negative_sum
		    yields_corr_norm = yields_corr[:,:]
		    yields_corr_norm[non_diag < 0.0] *= norm
		    yields_corr_norm[yields_corr_norm < -1.0] = -1.0
		    yields_cov_norm = df_std * yields_corr_norm * df_stdT
		elif( total_sum < 0.0 ):
		    non_diag = yields_cov_prenorm[:,:]
		    numpy.fill_diagonal( non_diag, 0.0 )
		    positive_sum = non_diag[non_diag > 0.0].sum()
		    norm = (total_sum + positive_sum) / positive_sum
		    yields_corr_norm = yields_corr[:,:]
		    yields_corr_norm[non_diag > 0.0] *= norm
		    yields_corr_norm[yields_corr_norm > 1.0] = 1.0
		    yields_cov_norm = df_std * yields_corr_norm * df_stdT
		#-------------------------------------------------------------------------------------



		#Verify conserved quantities 
		#-------------------------------------------------------------------------------------
		file = open( 'matrices/' + evaluation + '/verification_records/verification_' + system + '.csv', 'w' )
		file.write( 'Y_tot/Y_tot_orig, Z_tot/Z_CN, A_tot/(A_CN-nu_bar)\n' )
		for n in range(0,TRIALS):
		    Y_tot = sum( trials_res[:,n] ) / Y_tot_orig
		    Z_tot = sum( trials_res[:,n] * key_list_Z ) / (Zp*Y_tot*100.0)
		    A_tot = sum( trials_res[:,n] * key_list_A ) / ((Ap-nu_bar)*Y_tot*100.0)
		    file.write( str(Y_tot) + ', ' + str(Z_tot) + ', ' + str(A_tot) + '\n' )
		file.write( 'Matrix normalization:\n' )
		if( total_sum > 0.0 ):
			file.write( 'sign = -\n')
		else:
			file.write( 'sign = +\n' )
		file.write( 'norm = ' + str(norm) + '\n' )
		file.write( 'cov. sum = ' + str(yields_cov_norm.sum()) + '\n' )
		file.close()
		#-------------------------------------------------------------------------------------



		#Save the matrices to file
		#-------------------------------------------------------------------------------------
		#Primary covariance save
		file = open( 'matrices/' + evaluation + '/independent/' + system + '_cov.csv', 'w' )
		for key in key_list:
			file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
		file.write('\n')
		for i in range(0,len(key_list)):
			key = key_list[i]
			file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
			for j in range(0,len(key_list)):
				key2 = key_list[j]
				file.write( ', ' + str(yields_cov[i,j]) )
			file.write( '\n' )
		file.close()

		#Primary correlation save
		file = open( 'matrices/' + evaluation + '/independent/' + system + '_corr.csv', 'w' )
		for key in key_list:
			file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
		file.write('\n')
		for i in range(0,len(key_list)):
			key = key_list[i]
			file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
			for j in range(0,len(key_list)):
				key2 = key_list[j]
				file.write( ', ' + str(yields_corr[i,j]) )
			file.write( '\n' )
		file.close()

		#Normalized covariance save
		file = open( 'matrices/' + evaluation + '/independent/' + system + '_normed_cov.csv', 'w' )
		for key in key_list:
			file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
		file.write('\n')
		for i in range(0,len(key_list)):
			key = key_list[i]
			file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
			for j in range(0,len(key_list)):
				key2 = key_list[j]
				file.write( ', ' + str(yields_cov_norm[i,j]) )
			file.write( '\n' )
		file.close()

		#Normalized correlation save
		file = open( 'matrices/' + evaluation + '/independent/' + system + '_normed_corr.csv', 'w' )
		for key in key_list:
			file.write( ', ' + str(key[0]*10000 + key[2]*1000 + key[1]) )
		file.write('\n')
		for i in range(0,len(key_list)):
			key = key_list[i]
			file.write( str(key[0]*10000 + key[2]*1000 + key[1]) + ' ' )
			for j in range(0,len(key_list)):
				key2 = key_list[j]
				file.write( ', ' + str(yields_corr_norm[i,j]) )
			file.write( '\n' )
		file.close()
		#-------------------------------------------------------------------------------------



		#Plot the correlation matrix
		#-------------------------------------------------------------------------------------
		rcParams['figure.figsize'] = 14, 14
		rcParams['font.size'] = 22
		plt.matshow(yields_corr,cmap="RdBu",vmin=-1.0,vmax=1.0)
		cbar = plt.colorbar()
		cbar.set_label('Correlation Coeff.')
		plt.tight_layout()
		plt.xlabel('FY Index')
		plt.ylabel('FY Index')
		#plt.savefig( 'figures/' + system + '_corr.eps', format='eps', dpi=300 )
		plt.savefig( 'figures/' + evaluation + '/' + system + '_corr.png', dpi=500 )
		plt.clf()
		#-------------------------------------------------------------------------------------



		time_end = time.time()
		print( '   - ' + str(round(time_end-time_start)) + ' s.                            ' )
#---------------------------------------------------------------------------------------------