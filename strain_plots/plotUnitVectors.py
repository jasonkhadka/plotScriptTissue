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
import simulation_functions as sf
import argparse #argument parser, handles the arguments passed by command line
import gc

#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20.
plt.rcParams['ytick.labelsize'] = 20.
plt.rcParams['axes.labelsize'] =  20.
plt.rcParams['legend.fontsize'] = 20.
plt.rcParams['axes.titlesize'] =  20.

def getNeighbourFaces(cell,faceid):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	faceidarray = []
	while face != None:
		if face.getID() == faceid : 
			edges = qd.FaceEdgeIterator(face)
			edge = edges.next()
			while edge != None:
				faceidarray.append(edge.Right().getID())
				edge = edges.next()
			break
		face =faces.next()
	return faceidarray
########################################################################
def plotUnitVec(endstep,startstep=0,numOfLayer = 8,directory = '.'):

	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	# Getting Initial Area
	######################################################
	# Gathering face area
	######################################################
	totalMeanFaceArea = []
	m0determinantArray = []
	counter = startstep
	for i in range(startstep,endstep+1):
		if not os.path.isfile("qdObject_step=%03d.obj"%i):#check if file exists
			break
		percentStep = int((counter-startstep)/float(endstep-startstep)*100)
		sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
		sys.stdout.flush()
		cell = sf.loadCellFromFile(i)
		cell.setParameters()
		facearea = 0.
		tfmDet = 0.
		latdev.plotUnitVectors(cell,numOfLayer,  step = i, save = True, directory = directory)
		counter += 1
	return
	
#######################################################################
def plotEnergy(endstep,eta):
	energyarray = []
	for i in range(endstep+1):
		if not os.path.isfile("qdObject_step=%03d.obj"%i):#check if file exists
			break
		cell = sf.loadCellFromFile(i)
		cell.setOmega(.1)
		cell.setAlpha(1.)
		cell.setBeta(0.)
		cell.setGamma(0.0126)
		cell.setZeta(0.)
		cell.setSigma(0.)
		cell.setPressure(0.)
		cell.setConvexAngleThreshold(360)
		cell.setParameters()
		energyarray.append(cell.getEnergyCartesianVolume())
		#print eta,i, energyarray[-1], cell.getFirstTerm(), cell.getCartesianVolume(), cell.getBendingEnergy()
	tempfig = plt.figure(20)
	tax = tempfig.add_subplot(111)
	tax.plot(energyarray)
	#print energyarray
	os.chdir('..')
	tempfig.savefig('energy_eta=%.3f.png'%etacurrent,transparent='True')
	plt.close(20)
	return
	
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int,default = 8)
parser.add_argument("-a","--saveall",help = "to save all plot individually", action= "store_true")
parser.add_argument("-n","--eta", help = "target eta, default = 0",default = 0., type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
											  default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = 135, type = int)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)

## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
cylinder = args.cylinder
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
DIRNAMES=1
listdir = sorted(os.walk('.').next()[DIRNAMES])#all the list of directory in cwd
directoryName = r"unitvecplot_eta=%.2f"%ploteta
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]


totalfolders = len(listdir)
#print listdir
counter = 0
savedict = {}
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	if etacurrent != ploteta : continue
	#scaled color for plotting
	#print fkcurrent
	###########################################################
	#plotting and Saving
	############################################################
	#print os.listdir('.')
	os.chdir(folder)
	
	plotUnitVec(endStep,directory = saveDirectory)
	os.chdir("..")
	gc.collect()
	break

plt.close('all')

################################################################################
print "\n ################################################################################"
print " "*15l," DONE "
print "################################################################################"

