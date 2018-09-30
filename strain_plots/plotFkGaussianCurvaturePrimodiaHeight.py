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
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 22.
plt.rcParams['ytick.labelsize'] = 22.
plt.rcParams['axes.labelsize'] = 30.
plt.rcParams['legend.fontsize'] = 25.
plt.rcParams['axes.titlesize'] = 35.


def getEdge(cell, vert1id, vert2id):
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	while vertex != None:
		if vertex.getID() == vert1id:
			break
		vertex= vertices.next()
	###################################
	edges = qd.VertexEdgeIterator(vertex)
	edge = edges.next()
	while edge != None:
		if edge.Dest().getID() == vert2id:
			break
		edge = edges.next()
	return edge
###################################################################
def plotMinGaussianCurvaturePrimodiaHeight(numOfLayer, targetid,endstep,fastkappa, mincurvatureplot, heightplot,color):
	import numpy as np
	import matplotlib.pyplot as plt
	########################################################################
	# Loading the Cell #
	########################################################################
	#######################################################################
	# Checking if the files exists if not going to step down
	#######################################################################
	meanGaussianCurvature = []
	primodialheight = []
	for step in range(1,endStep+1):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
			break
		cell = sf.loadCellFromFile(step)
		########################################################################
		#                 Starting the Calcuation of Curvature                 #
		########################################################################
		gaussianCurvature = []
		meanpointx = []
		meanpointy = []
		meanpointz = []
		##########
		edge = getEdge(cell,211,210)
		edgenext = edge
		####grabbing the origin of edge####
		vertexDest = edge.Dest()
		meanpointx.append(vertexDest.getXcoordinate())
		meanpointy.append(vertexDest.getYcoordinate())
		meanpointz.append(vertexDest.getZcoordinate())
		# Get the Gaussian Curvature
		gaussianCurvature.append(vertexDest.getGaussianCurvature())
		while True:
			edgenext = edgenext.Rprev()
			vertexDest = edgenext.Dest()
			meanpointx.append(vertexDest.getXcoordinate())
			meanpointy.append(vertexDest.getYcoordinate())
			meanpointz.append(vertexDest.getZcoordinate())
			gaussianCurvature.append(vertexDest.getGaussianCurvature())
			########################################################
			edgenext = edgenext.Lnext()
			vertexDest = edgenext.Dest()
			meanpointx.append(vertexDest.getXcoordinate())
			meanpointy.append(vertexDest.getYcoordinate())
			meanpointz.append(vertexDest.getZcoordinate())
			gaussianCurvature.append(vertexDest.getGaussianCurvature())
			########################################################
			edgenext = edgenext.Lnext()
			vertexDest = edgenext.Dest()
			meanpointx.append(vertexDest.getXcoordinate())
			meanpointy.append(vertexDest.getYcoordinate())
			meanpointz.append(vertexDest.getZcoordinate())
			gaussianCurvature.append(vertexDest.getGaussianCurvature())
			###############################################################
			edgenext = edgenext.Rprev()
			vertexDest = edgenext.Dest()
			meanpointx.append(vertexDest.getXcoordinate())
			meanpointy.append(vertexDest.getYcoordinate())
			meanpointz.append(vertexDest.getZcoordinate())
			gaussianCurvature.append(vertexDest.getGaussianCurvature())
			########################################################
			edgenext = edgenext.Lnext()
			vertexDest = edgenext.Dest()
			meanpointx.append(vertexDest.getXcoordinate())
			meanpointy.append(vertexDest.getYcoordinate())
			meanpointz.append(vertexDest.getZcoordinate())
			gaussianCurvature.append(vertexDest.getGaussianCurvature())
			#############################################################
			if edgenext.Org().getID() == edge.Org().getID():
				break            
		facetarget = sf.getFace(cell,targetid)
		targetx = facetarget.getXCentralised()
		targety = facetarget.getYCentralised()
		targetz = facetarget.getZCentralised()
		meanx = np.mean(meanpointx)
		meany = np.mean(meanpointy)
		meanz = np.mean(meanpointz)
		height = np.sqrt((meanx-targetx)**2+(meany-targety)**2+(meanz-targetz)**2)
		##########################################################################
		primodialheight.append(height)
		meanGaussianCurvature.append(np.mean(np.array(gaussianCurvature)))
	#################################################################################
	#                         Plotting 
	#################################################################################
	timestep = range(len(primodialheight))
	#print primodialheight, meanGaussianCurvature
	# Min Gaussian curvature
	mincurvatureplot.plot(timestep,meanGaussianCurvature,c=color,lw = 3)
	###################################
	# Height of Primodia
	heightplot.plot(timestep, primodialheight, c=color,lw = 3)
	########################################################
	return


####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int,default = 8)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
											  default = 0.8, type = float)
parser.add_argument("-n","--nonNormalize", help = "if option is used, the figures are not normalised", action= "store_false")
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
directoryName = "gaussianPlot"
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
fklist = [float(dict(item.split("=") for item in folder.split("_"))['fk']) for folder in listdir]
#################################################################################
#     Making the Plots
#################################################################################
import matplotlib.colors as colors
import matplotlib.cm as cmx
jet = cm = plt.get_cmap('viridis') 
maxvalue = max(fklist)
minvalue = min(fklist)
#print "Max value", maxvalue, " minvalue", minvalue
cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#fig = plt.figure(frameon=False,figsize=(20,16))
fig = plt.figure(figsize=(20,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
ax1.set_title("Mean Gaussian Curvature")
ax1.set_xlabel("time")
ax1.set_ylabel("Mean Gaussian Curvature")
###################################
# Height of Primodia
###################################
ax2.set_title("Height of Primodia")
ax2.set_xlabel("time")
ax2.set_ylabel("Height of Primodia")
########################################################
counter = 0
totalfolders = len(listdir)
#print listdir
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	#scaled color for plotting
	#print fkcurrent
	fkcolor = scalarMap.to_rgba(fkcurrent)
	percentStep = int((counter)/float(totalfolders)*100)
	sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
	sys.stdout.flush()
	###########################################################
	#plotting and Saving
	############################################################
	#print os.listdir('.')
	os.chdir(folder)
	#print float(folderdict['n'])
	#print "\n",os.getcwd()
	plotMinGaussianCurvaturePrimodiaHeight(numOfLayer = numOfLayer, targetid = targetid,endstep = endStep,fastkappa = fkcurrent,mincurvatureplot = ax1, heightplot=ax2,color = fkcolor)
	#plotMeanGrowthRate(step=endStep,azim = azim, elev = elev,eta = float(folderdict['n']),save=True)
	os.chdir("..")
	counter+= 1
###############################################################################
#color bar
scalarMap._A = []
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
clrbar.set_label("Fast Kappa")
#plt.tight_layout()
#plt.tight_layout( rect=[0, 0, 1, 1])
plt.savefig(saveDirectory+r"/plot_curvature_and_height.png",transparent = True)
plt.close('all')
################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

