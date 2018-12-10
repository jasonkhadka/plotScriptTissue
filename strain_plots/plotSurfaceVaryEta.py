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
import argparse #argument parser, handles the arguments passed by command line
import gc
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22

####################################################################################################################
# Calculating the max time step for target surface area
####################################################################################################################
def getTimeStep(targetArea, endStep, startStep=1, stepsize = 10):
	####################################################
	for step in range(startStep, endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return endStep,0.
		################################################
		cell = sf.loadCellFromFile(step)
		################################################
		tissueSurfaceArea = sf.getSurfaceArea(cell)
		if (tissueSurfaceArea > targetArea):
			gc.collect()
			for calstep in range(step-1,step-stepsize-1,-1):
					cell = sf.loadCellFromFile(calstep)
					tissueSurfaceArea = sf.getSurfaceArea(cell)
					if (tissueSurfaceArea <= targetArea):
						gc.collect()
						cell = sf.loadCellFromFile(calstep+1)
						tissueSurfaceArea = sf.getSurfaceArea(cell)
						return calstep+1,tissueSurfaceArea
		################################################
		gc.collect()
	return endStep,tissueSurfaceArea
##########################################################################################
#       Function to Plot Magnitude of Normal Forces on the Faces of the Cell
##########################################################################################
def plotSurface(cell,ax,color,surface=False,alpha = 0.6):
	 from mpl_toolkits.mplot3d import Axes3D
	 import matplotlib as mpl
	 from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	 import numpy as np
	 import matplotlib.pyplot as plt
	 ###############################################################
	 #                 Plotting the Cell                           #
	 ###############################################################
	 faces = qd.CellFaceIterator(cell)
	 ###################
	 face = faces.next()
	 faceCounter = 0
	 xcenarray = []
	 ycenarray = []
	 zcenarray = []
	 while (face != None):
		  if face.getID()==1:
				face  = faces.next()
				continue
		  faceid = face.getID()#grabbing face id
		  xlist = []
		  ylist = []
		  zlist = []
		  #print "== Face ID : ", faceid, "=="
		  edges = qd.FaceEdgeIterator(face)
		  edge = edges.next()
		  while edge != None:
				####grabbing the origin of edge####
				#centralised coordiante
				vertex = edge.Org()
				#print vertex.getID()
				xCoord1 = vertex.getXcoordinate()
				yCoord1 = vertex.getYcoordinate()
				zCoord1 = vertex.getZcoordinate()
				xlist.append(xCoord1)
				ylist.append(yCoord1)
				zlist.append(zCoord1)
				edge = edges.next()
		  xlist.append(xlist[0])
		  ylist.append(ylist[0])
		  zlist.append(zlist[0])
		  verts = [zip(xlist, ylist,zlist)]
		  if surface:
		  	pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		  	pc.set_edgecolor(color)
		  	ax.add_collection3d(pc)
		  else:
		  	ax.plot(xlist,ylist,zlist, c = color, lw = 3)
		  face = faces.next()
		  faceCounter+= 1
	 return
################################################################################################
#        Plotting Part 
################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step",default = 2000, type = int)
parser.add_argument('-d',"--stepsize", help="step size",default =10, type = int)
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-c","--targetarea",help = "targetarea for surface",
											  default = 800, type = float)
parser.add_argument("-n","--nonNormalize", help = "if option is used, the figures are not normalised", action= "store_false")
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 0., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
											  default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = 10, type = float)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = 70, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 0, type = float)
parser.add_argument("-f","--fastkappa", help = "if option is used, the figures are made with respect to chaning fast kappa", action= "store_true")
parser.add_argument("-x","--surface", help = "if option is used, surface will be plotted", action= "store_true")

## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
targetArea = args.targetarea
stepsize = args.stepsize
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
surface  = args.surface
fastkappaOption = args.fastkappa

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as ticker
#import the libraries
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
#limits of the plot
Length = 1.
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
#plotting part

#################################################################
# making the plot
#################################################################
fig = plt.figure(frameon=False,figsize=(10,10))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = Axes3D(fig)
xlim = 0.5
ax.set_xlim((-xlim*radius,xlim*radius))
ax.set_ylim((-xlim*radius,xlim*radius))
ax.set_zlim((-0.1*radius,0.7*radius))
ax.axis('off')
ax.xaxis.pane.set_edgecolor('black')
ax.yaxis.pane.set_edgecolor('black')
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False


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
else:
	etalist = sorted([float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir])

################################################################################
##################################################
counter = 0
totalfolders = len(listdir)
plotData = {}
#print listdir
colors = ['darkgreen','navy']
targeteta = [0.,8.]
from matplotlib.lines import Line2D
legenditems = []
for folder in listdir:
	if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	   etacurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	   etacurrent= growthRatio[etacurrent]
	else:
	   etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	########################################################
	if not etacurrent in targeteta:
		continue
	color = colors[counter]
	########################################################
	percentStep = int((counter)/float(totalfolders)*100)
	sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
	sys.stdout.flush()
	###########################################################
	#plotting and Saving
	###########################################################
	os.chdir(folder)
	step,tissueSurfaceArea = getTimeStep(targetArea, endStep, startStep, stepsize = stepsize)
	cell = sf.loadCellFromFile(step)
	plotSurface(cell,ax,color,surface,zorder = 10-counter)
	legenditems.append(Line2D([0], [0], color=color, label=r"$\eta=%d$"%etacurrent,lw='3'))
	os.chdir("..")
	gc.collect()
	counter+= 1

ax.legend(handles = legenditems)
ax.view_init(azim = azim, elev = elev)
if surface:
	fig.savefig(saveDirectory+r"/plot_surface_surf_eta=%d_%d.png"%(targeteta[0],targeteta[1]),transparent = True, bbox_inches="tight")
	fig.savefig(saveDirectory+r"/plot_surface_surf_eta=%d_%d.eps"%(targeteta[0],targeteta[1]),transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plot_surface_eta=%d_%d.png"%(targeteta[0],targeteta[1]),transparent = True, bbox_inches="tight")
	fig.savefig(saveDirectory+r"/plot_surface_eta=%d_%d.eps"%(targeteta[0],targeteta[1]),transparent = True, bbox_inches="tight")



plt.close('all')

################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"


