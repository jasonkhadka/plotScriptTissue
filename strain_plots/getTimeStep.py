import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jkhadka/plantdev")
import quadedge as qd
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
#sys.path.append("/home/jkhadka/transferfile/scripts/simulation_functions")
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
sys.path.append('/home/jkhadka/transferdata/scripts/plotscript/')
import simulation_functions as sf
import pandas as pd
import argparse #argument parser, handles the arguments passed by command line
import gc
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 22

import numpy as n, pylab as p, time
import scipy.interpolate as interpolate
from scipy.spatial import ConvexHull

####################################################################################################################
# Calculating the max time step for target surface area
####################################################################################################################
def getTimeStep(targetArea, endStep, startStep=1, stepsize = 10):
	####################################################
	for step in range(startStep, endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return endStep,0.
		################################################
		cell = sf.loadCellFromFile(step,resetids = True)
		################################################
		tissueSurfaceArea = sf.getSurfaceArea(cell)
		if (tissueSurfaceArea > targetArea):
			gc.collect()
			for calstep in range(step-1,step-stepsize-1,-1):
					cell = sf.loadCellFromFile(calstep,resetids = True)
					tissueSurfaceArea = sf.getSurfaceArea(cell)
					if (tissueSurfaceArea <= targetArea):
						gc.collect()
						cell = sf.loadCellFromFile(calstep+1,resetids = True)
						tissueSurfaceArea = sf.getSurfaceArea(cell)
						return calstep+1,tissueSurfaceArea
		################################################
		gc.collect()
	return endStep,tissueSurfaceArea

####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--startarea", help="Start of simulation area",default =None, type = int)
parser.add_argument('-e',"--endarea", help="End of simulation area",default =850,  type = int)
parser.add_argument("-m","--maxeta", help = "if this is given, then eta is only cacluated till this value", type = float, default = 0.0)
parser.add_argument("-c","--surfacearea", help = "calculate time step for the surfacearea", type = float, default = 850.)
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int,default = 8)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
											  default = 1., type = float)
parser.add_argument("-n","--nonNormalize", help = "if option is used, the figures are not normalised", action= "store_false")
parser.add_argument("-L","--Large", help = "if option is used, calculation is done for larger primordia", action= "store_true")
parser.add_argument("-f","--fastkappa", help = "if option is used, the figures are made with respect to chaning fast kappa", action= "store_true")
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
parser.add_argument('-d',"--areastep", help="area step for calculating the growth in cell area", type = int,default = 20)
parser.add_argument('-j',"--jobid", help="jobid", type = int,default = None)
## Getting the arguments 
args = parser.parse_args()
#location = args.location
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
anglethreshold = args.angle
targetface = args.target
azim = args.azimuthal
elev = args.elevation
norm = args.nonNormalize
targetid = args.target
areastep = args.areastep
maxeta = args.maxeta
fastkappaOption = args.fastkappa
large  = args.Large
stepsize = 10
surfacearea = args.surfacearea
startarea = args.startarea
endarea = args.endarea
jobid = args.jobid


startStep = 1
endStep = 2000
# For surpressing err
class NullDevice():
	def write(self, s):
		pass
if targetface == None:
	if numOfLayer == 8:
		targetface = 135
	elif numOfLayer == 10:
		targetface = 214

#print " start "
#original_stdout = sys.stderr # keep a reference to STDOUT

#sys.stderr = NullDevice()  # redirect the real STDOUT

################################################################################
import sys
import os
#################################
cwd = os.getcwd()#getting the current working directory
DIRNAMES=1
listdir = sorted(os.walk('.').next()[DIRNAMES])#all the list of directory in cwd
directoryName = "plots"
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['fk']) for folder in listdir]
	column1 = 'fk'
else:
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
	column1 = 'eta'
#################################################################################
#     Making the table
########################################################
counter = 0
totalfolders = len(listdir)
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	timedata = pd.DataFrame(columns = ['fk','timestep', 'area'])
else:
	timedata = pd.DataFrame(columns = ['eta','timestep', 'area'])
#print listdir
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	else:
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	########################################################
	if (maxeta != 0) and (etacurrent > maxeta):
		continue
	########################################################
	percentStep = int((counter)/float(totalfolders)*100)
	sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
	sys.stdout.flush()
	###########################################################
	#plotting and Saving
	############################################################
	#print os.listdir('.')
	os.chdir(folder)
	############################################################
	step,tissueSurfaceArea = getTimeStep(surfacearea, endStep, startStep, stepsize = stepsize)
	timedata = timedata.append({column1:etacurrent,
					'timestep':step,
					'area':tissueSurfaceArea},ignore_index=True)
	############################################################
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#color bar fig
###############################################################################
timedata = timedata.sort_values(by=column1)
###############################################################################
print '\n','#'*100
print timedata