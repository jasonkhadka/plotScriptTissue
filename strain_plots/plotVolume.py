########################################################################
# To be ran in the data folder
########################################################################

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jkhadka/plantdev")
import quadedge as qd
import gc
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
#sys.path.append("/home/jkhadka/transferfile/scripts/simulation_functions")
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
import simulation_functions as sf
import argparse #argument parser, handles the arguments passed by command line
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 18.
plt.rcParams['ytick.labelsize'] = 18.
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.titlesize'] = 18

####################################################################################################################
# Calculating the max time step for target surface area
####################################################################################################################
def plotVolume(endStep, startStep, targetarea,volplot1,volplot2,stepsize = 5):
	####################################################
	volumeArray = []
	surfaceAreaArray = []
	timeArray = []
	####################################################
	for step in range(startStep, endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return 
		################################################
		cell = sf.loadCellFromFile(step)
		################################################
		tissueSurfaceArea = sf.getSurfaceArea(cell)
		tissueVolume = cell.getCartesianVolume()
		################################################
		if (tissueSurfaceArea > targetArea):
			return
		################################################
		volumeArray.append(tissueVolume)
		surfaceAreaArray.append(tissueSurfaceArea)
		timeArray.append(step)
		################################################
	################################################
	gc.collect()
	################################################
	volplot1.plot(timeArray, volumeArray, linewidth = 3,c='salmon')
	volplot2.plot(surfaceAreaArray, volumeArray, linewidth = 3,c='c')
	################################################
	return 
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-d","--stepsize", help = "stepsize on plot", default = 1, type = int)
parser.add_argument('-c',"--targetsurfaceArea", help="The surface area for which the stress vs feedback plot would need to be plotStrainDifferenceSurface",
					default =None, type = float)

parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int,default = 8)
parser.add_argument("-a","--saveall",help = "to save all plot individually", action= "store_true")
parser.add_argument("-n","--eta", help = "if option is used, plot is varying eta, else plot is vary fastkappa", action= "store_false")
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
											  default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = None, type = int)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)

## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
stepsize = args.stepsize
targetarea = args.targetsurfaceArea
#alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
anglethreshold = args.angle
targetface = args.target
azim = args.azimuthal
elev = args.elevation
#norm = args.nonNormalize
ploteta = args.eta
targetid = args.target
saveall = args.saveall
# For surpressing err
class NullDevice():
	def write(self, s):
		pass


#print " start "
#original_stdout = sys.stderr # keep a reference to STDOUT

#sys.stderr = NullDevice()  # redirect the real STDOUT

################################################################################
import sys
import os
#################################
cwd = os.getcwd()#getting the current working directory
saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)
surfaceSaveDirectory = saveDirectory+r"/volumePlot"
if not os.path.exists(surfaceSaveDirectory):
		os.makedirs(surfaceSaveDirectory)
#################################################################################
#     Making the Plots
#################################################################################
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
########################################################


Length = 1.
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on


#######################################################
fig1 = plt.figure(1,figsize=(5,5))
ax1 = fig1.add_subplot(111)
ax1.set_ylabel('Volume of tissue, $V_T$')
ax1.set_xlabel('time step, $t$')
#######################################################
fig2 = plt.figure(2,figsize=(5,5))
ax2 = fig2.add_subplot(111)
ax2.set_ylabel('Volume of tissue, $V_T$')
ax2.set_xlabel('Surface area, $A_T$')
#######################################################
plotVolume(endStep, startStep, targetarea=targetarea,volplot1 = ax1,volplot2 = ax2)
################################################################################
fig1.savefig(surfaceSaveDirectory+'/'+'volumePlot_time=%03d.png'%endStep,
	bbox_inches='tight',dpi=100,transparent = True)
fig2.savefig(surfaceSaveDirectory+'/'+'volumePlot_area=%d.png'%targetarea,
	bbox_inches='tight',dpi=100,transparent = True)
################################################################################
print "\n"+ "#"*100
print " "*45+"DONE "
print "#"*100
################################################################################
plt.close('all')