####################################################################################################################################
#                      To plot Strain & Stress vectors From given coordinates and target Form Matrix
####################################################################################################################################
import matplotlib
import numpy as np
matplotlib.use('agg')
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
##########################################################################################
# if simulation_functions.py is not present, copying the file
##########################################################################################
if not os.path.isfile('simulation_functions.py'):
	import shutil as st
	st.copy2("/home/jkhadka/transferdata/scripts/simulation_functions/simulation_functions.py",os.getcwd())
import simulation_functions as sf
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--steps", help="The number of files there are (number of steps in simulation)",type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 1.",
								   default = 1., type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0.",
								   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 1.0",
								   default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.0",
								   default = 0.0, type = float)
parser.add_argument("-k","--kappa",help = "value to set for kappa growth Rate of Mo, default = 1.0",
								   default = 0.01, type = float)
parser.add_argument("-n","--angle",help = "value to set for convex angle threshold, default = 360",
								   default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)
## Getting the arguments 
args = parser.parse_args()
#location = args.location
simulationStep = args.steps
cylinder = args.cylinder
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
kappa = args.kappa
anglethreshold = args.angle
#################################################################################
#location = "" #the location of files # EDIT THIS EVERYTIME (at least for now :( )
#os.chdir("location")
#################################################################################
#####################################################################################################################
#####                 MAKING THE INITIAL QUADEDGE TISSUE
#####################################################################################################################
cell, newedge = latdev.makeCenteredHexagonLattice(numOfLayer, True, False) 
####################################################################################
#put down the paramters to cell
####################################################################################
cell.setAlpha(alpha)
cell.setBeta(beta)
cell.setPressure(pressure)
cell.setZeta(zeta)
cell.setGamma(gamma)
cell.setKappa(kappa)
cell.setConvexAngleThreshold(anglethreshold)
## parameters printing
print "Alpha : ", alpha
print "beta : ", beta
print "pressure : ", pressure
print "zeta : ", zeta
print "gamma :", gamma
print "kappa :", kappa
print "angle threshold :" , anglethreshold
####################################################################################
# Making Initial Condition 
#			a. args.cylinder == False -> Make Dome 
#			b. args.cylinder == True -> Make Cylinder
####################################################################################
if args.cylinder: 
	cell = sf.makeCylinder(cell,numOfLayer)
else:
	cell = sf.makeDome(cell,numOfLayer)
####################################################################################
#cell.setCylindrical()#SETTING Cylindrical coordinates for all the vertices
#cell.setInitialParameters()
#numberofvertices = cell.countVertices()
#bendingThreshold = cell.getBendingThreshold()
#plotStressSurface(cell, numOfLayer, step=-1)
#plotStrainSurface(cell,numOfLayer, step =-1)
####################################################################################
# Now Plotting the Stress and Strain Plots
####################################################################################
for step in range(1,simulationStep+1):
		#print " Step : ", step
		cell = sf.setCartesianCoordinate(cell,step)## Loadding and setting coordinates
		cell = sf.setTargetFormMatrix(cell,step-1) ## Loadding and setting target From Matrix
		cell.setParameters()
		cell.calculateVertexForce()
		cell.calculateStressStrain()
		sf.plotPrimaryStrainSurface(cell,numOfLayer,step=step)
################################################################################
print "################################################################################"
print " DONE "
print "################################################################################"
