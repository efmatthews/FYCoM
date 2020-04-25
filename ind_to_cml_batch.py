#Eric Matthews
#September 23, 2019
#Fission Yield Covariance Matrix Generation (FYCoM)
#Independent covariances to cumulative covariances



#Constants
#---------------------------------------------------------------------------------------------
TRIALS = 10000
energies = { 'T':'thermal', 'F':'fission', 'H':'DT', 'D':'DD', 'SF':'SF' }
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
import time
from multiprocessing import Pool, cpu_count
from ind_to_cml import main
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



#Run FIER to generate decay stems for each system
#-----------------------------------------------------------------------------------------
for p in range(0,len(systems)):
	system = systems[p]
	Zp = ZAs[p][0]
	Ap = ZAs[p][1]
	file = open( 'FIER/FIER/deck.txt', 'w' )
	file.write( 'MODE:SINGLE\n' )
	file.write( 'ON DECAY PREDICTION\n' )
	file.write( 'input_data/isotopes.csv     ISOTOPES FILE\n' )
	file.write( 'input_data/decays.csv       DECAYS   FILE\n' )
	file.write( 'input_data/gammas.csv       GAMMAS   FILE\n' )
	file.write( 'YIELDS:ER\n' )
	if( system[-2:] == 'SF' ):
		energy = 'SF'
		A_in = Ap
	else:
		energy = energies[system[-1]]
		A_in = Ap-1
	file.write( elements[Zp] + ',' + str(A_in) + ',' + energy + '   YIELDS   FILE\n' )
	file.write( 'NONE  CHAINS OUTPUT\n' )
	file.write( 'output/decay_stems_' + str(system) + '.csv   STEMS OUTPUT\n' )
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


#Iterate over systems, covert covariance matrix for each
#---------------------------------------------------------------------------------------------
args = []
for p in range(0,len(systems)):
	args.append( ( systems[p], ZAs[p][0], ZAs[p][1], TRIALS ) )
p = Pool( cpu_count() )
p.map( main, args )
#---------------------------------------------------------------------------------------------