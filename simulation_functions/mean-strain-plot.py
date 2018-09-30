####################################################################################################################################
#               Plotting the Mean Strain Plot 
####################################################################################################################################
#importing all the libraries
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
#pd.set_option('precision',10)
#pd.set_option("display.max_columns",200)
#pd.set_option("display.max_rows",2000)
import sys 
sys.path.append("/home/jkhadka/plantdev")
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
import quadedge as qd
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
#from IPython.display import display, HTML
import nlopt #non linear optimizer
import argparse #argument parser, handles the arguments passed by command line
import os # to make directory
import time
#importing functions that the simulations rely on -- make sure to copy simulation_functions.py to the same directory as this file
from simulation_functions import *
import re
######setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#adding arguments
parser.add_argument("-l","--layer",help = "number of layers of cells to simulate",
                                   default = 5, type = int)
args = parser.parse_args()#arguments passed
numOfLayer = args.layer
#############################################################
#key to sort the file_names in order
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
########################################################################################

#listing all the coordinate files
coordinate_files = sorted((fn for fn in os.listdir('.') if fn.startswith('coordinates_')), key = numericalSort)
#################################################
#Final Step Number of Simulation
finalnumber = int(re.findall('\d+',coordinate_files[-1])[0])
#print finalnumber
#print coordinate_files
#####################################################################################################################
#####                 MAKING THE INITIAL QUADEDGE TISSUE
#####################################################################################################################
#    Function to make the tissue
#		a. makeCenteredHexagonLattice(numOfLayer, printcondition = False, projection = True) 
#	 Projection condition : it is to have the external vertex projected on to an external circle
##############################################################################################
cell, newedge = latdev.makeCenteredHexagonLattice(numOfLayer, True, False) 

######################################################################
## Now iterating with coordinates and plotting the mean strain
######################################################################
meanstrain = np.zeros(finalnumber+1)
#
for stepcount in range(finalnumber+1):
	setCartesianCoordinate(cell,stepcount)
	setTargetFormMatrix(cell,stepcount)
	####################################
	cell.calculateStressStrain()
	meanstrain[stepcount] = cell.getMeanStrainDeterminant()
	####################################
plt.figure()
plt.plot(range(finalnumber+1), meanstrain,'-x')
plt.savefig("mean_strain_determinant.png",transparent=True)
plt.clf()
plt.close()
