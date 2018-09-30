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
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20.
plt.rcParams['ytick.labelsize'] = 20.
plt.rcParams['axes.labelsize'] =  20.
plt.rcParams['legend.fontsize'] = 20.
plt.rcParams['axes.titlesize'] =  20.
def getNeighbourFaces(cell,faceid):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	faceidarray = [faceid]
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

def plotStrainDeterminant(cell,saveDirectory,step,targetID = 135,numOfLayer = 8):
	cell.calculateStrain()
	faceidarray = getNeighbourFaces(cell, targetID)
	#################################################
	faces =qd.CellFaceIterator(cell)
	face = faces.next()
	straindet = []
	faceids = []
	strainsum = 0.
	straincount = 0.
	while face != None:
	    if face.getID() ==1:
	        face = faces.next()
	        continue
	    straindet.append(face.getStrainDeterminant())
	    faceids.append(face.getID())
	    if not face.getExternalPosition():
	    	# excluding the boundary faces in mean
	    	strainsum += face.getStrainDeterminant()
	    	straincount += 1.
	    face = faces.next()
	strainMean = strainsum/straincount
	####################################
	#        Plotting Part
	#################################### 
	fig = plt.figure(1,figsize=(16,8))
	# Histogram
	ax1 = fig.add_subplot(121)
	ax1.hist(straindet,60,alpha =0.8)
	ax1.set_xlabel("strain determinant")
	ax1.set_ylabel("Frequency")
	ax1.set_title("Determinant of Strain, mean = %.5f"%strainMean)
	# face id vs straindeterminant
	ax2 = fig.add_subplot(122)
	ax2.set_ylabel("strain determinant")
	ax2.set_xlabel("faceid")
	ax2.set_title("Determinant of Strain vs Face ID")
	for f,s in zip(faceids,straindet):
	    if f in faceidarray:
	        s1 = ax2.scatter(f,s,c = 'g',marker = '*',s = 40, label = "fast growing")
	    else:
	        s2 = ax2.scatter(f,s,c = 'r',marker = '<',s = 40, label = "slow growing")
	fig.legend((s1,s2),("fast growing","slow growing"),bbox_to_anchor=(0.87,0.35))
	plt.tight_layout()
	plt.savefig(saveDirectory+r"/strainHistogram=%03d.png"%(step),transparent=True)
	plt.close('all')
	return


################################################################################################

################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int,default = 1)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
											  default = 0.8, type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-n","--eta", help = "target eta, default = 0",default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
											  default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = 10, type = float)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)

## Getting the arguments 
args = parser.parse_args()
#location = args.location
step = args.end
startStep = args.start
cylinder = args.cylinder
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
ploteta = args.eta

#################################################################################
import sys
import os
cwd = os.getcwd()#getting the current working directory
DIRNAMES=1
listdir = sorted(os.walk('.').next()[DIRNAMES])#all the list of directory in cwd
##################################################
directoryName = r"strainHistogram_eta=%.2f"%ploteta
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
######################################################################################
listdir = [d for d in listdir if d[0] == 'a']
folder = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
#######################################################
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	if etacurrent != ploteta : continue
	os.chdir(folder)
	if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
				print "File Not Found : qdObject_step=%03d.obj"%step
				quit()
	cell = sf.loadCellFromFile(step)
	cell.calculateStrain()
	plotStrainDeterminant(cell,saveDirectory,step,targetID = 135)
	break

################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

