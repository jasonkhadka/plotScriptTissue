####################################################################################################################################
#                      Plotting Relation of Normal Force on Faces and their Area
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

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--steps", help="The number of files there are (number of steps in simulation)",type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-m', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
                                   default = 1., type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
                                   default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
                                   default = 0.0, type = float)
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
anglethreshold = args.angle
#################################################################################
#location = "" #the location of files # EDIT THIS EVERYTIME (at least for now :( )
#os.chdir("location")
#################################################################################
#######  Geometrical parameters of the tissue  #########
Length = 1.# length of the sides of the faces
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
numofface = 3*numOfLayer*(numOfLayer-1)+1
perimeterVertexNum = 6+12*(numOfLayer-1) #number of vertices in the perimeter of lattice
perimeterHexNum = (numOfLayer!=1)*(6*(numOfLayer-1)-1) + 1 # number of perimeter hexgons
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
#cell.setKappa(kappa)
cell.setConvexAngleThreshold(anglethreshold)
variancearray=[]
## parameters printing
print "Alpha : ", alpha
print "beta : ", beta
print "pressure : ", pressure
print "zeta : ", zeta
print "gamma :", gamma
print "angle threshold :" , anglethreshold
####################################################################################
# Making Initial Condition 
#			a. args.cylinder == False -> Make Dome 
#			b. args.cylinder == True -> Make Cylinder
####################################################################################
if args.cylinder: 
	cell = makeCylinder(cell,numOfLayer)
	cellradius = radius/1.5
	cell.setCylindrical()#SETTING Cylindrical coordinates for all the vertices
	for step in range(simulationStep):
		#print " Step : ", step
		cell = setCylindrical(cell,step)## Loadding and setting coordinates
		cell = setTargetFormMatrix(cell,step) ## Loadding and setting target From Matrix
		cell.setParameters()
else:
	cell = makeDome(cell,numOfLayer)
	cellradius = radius
	for step in range(simulationStep):
		print " Step : ", step
		cell = setCartesianCoordinate(cell,step)## Loadding and setting coordinates
		cell = setTargetFormMatrix(cell,step) ## Loadding and setting target From Matrix
		cell.setParameters()
		## PLOTTING THE SURFACE ##
		plotNormalForce(cell,numOfLayer,step = step)
		plotFaceArea(cell,numOfLayer,step = step)
		plotForceAreaRelation(cell,numOfLayer,step = step)
		variancearray = plotAreaDistribution(cell,numOfLayer,variancearray = variancearray,step = step)
################################################################################
print "################################################################################"
print " DONE "
print "################################################################################"
