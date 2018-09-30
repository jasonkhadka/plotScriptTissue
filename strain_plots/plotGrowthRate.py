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
def getFaceAreaArray(totalstep):
	# Preparing the facearea dictionary
	step = 0
	cell = sf.loadCellFromFile(step)
	facearea = {}
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		"""if sf.checkExternalFace(face):
			face  = faces.next()
			continue
		"""
		facearea[face.getID()] = []
		face = faces.next()
	#############################################
	# Loading the cells for steps till totalstep
	#############################################
	for step in range(totalstep+1):
		cell = sf.loadCellFromFile(step)
		cell.setParameters()
		cell.calculateStrain()
		########################################################
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == 1:
				face  = faces.next()
				continue
			"""if sf.checkExternalFace(face):
				face  = faces.next()
				continue
			"""
			facearea[face.getID()].append(face.getAreaOfFace())
			face = faces.next()
	return facearea

def getRates(array):
	rates = []
	for i in range(len(array)-1):
		rates.append((array[i+1]-array[i])/array[i])
	return rates

def plotMeanGrowthRate(step, numOfLayer=8,eta = 0,targetface =10, alpha = 0.8, Length=1.0, save=False,azim = -70, elev=50):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(10,8))
	#fig = plt.figure(frameon=False)
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	#ax.xaxis.pane.fill = False
	#ax.yaxis.pane.fill = False
	#ax.zaxis.pane.fill = False
	#######################################################################
	# Checking if the files exists if not going to step down
	try:
		#cell = qd.objReadCell("qdObject_step=%03d.obj"%step)
		cell = sf.loadCellFromFile(step)
		#print "success"
	except:
		#print " exception reached"
		step -= 1
		#sys.stdout.write("step : ", step , os.getcwd())
		#sys.stdout.flush()
		#print "step : ", step , os.getcwd()
		plotMeanGrowthRate(step,azim = azim, elev = elev,eta = eta,save = save)
		return
	#print "Step : ", step
	########################################################################
	#                    Loading the step-1 step              #
	########################################################################
	completeFaceArray = getFaceAreaArray(step)
	rateDict = {}
	for key in completeFaceArray:
		value = np.array(completeFaceArray[key])
		rateDict[key] = np.mean(np.array(getRates(value)))
	####################################################################################
	# PLotting the faces now
	####################################################################################
	cell = qd.objReadCell("qdObject_step=%03d.obj"%step)
	# Flipping the Face ID after Loading to cell so that the face id before and after matches
	faces = qd.CellFaceIterator(cell)
	facecount = cell.countFaces()
	face= faces.next()
	while face != None:
		faceid = face.getID()
		face.setID(facecount-faceid+1)
		#print face.getID()
		face = faces.next()
	######################################################
	#settig target Form Matrix : 
	#TMF step corresponds to coordinate step
	####
	cell = sf.setTargetFormMatrix(cell,step)
	cell.setParameters()
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStrain()
	########    ########    ########    ########    ########
	#                 Plotting the Cell                    #
	########    ########    ########    ########    ########
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	maxvalue = 0.006
	minvalue = 0.
	#print "Max value", maxvalue, " minvalue", minvalue
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	counter = 0
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
		#print "== Face ID : ", faceid, "=="
		xmean = face.getXCentralised()
		ymean = face.getYCentralised()
		zmean = face.getZCentralised()
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
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print "face ID : ", face.getID(), " ratio : ", ratio
		#ratio = (face.getStrainTrace()+minTrace)/(minTrace+maxTrace)
		"""if sf.checkExternalFace(face):
			ratio = 0.
			color = scalarMap.to_rgba(ratio)
		else:
		"""
		ratio = rateDict[face.getID()]
		#print ratio
		color = scalarMap.to_rgba(ratio)
		#print ratio
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	ax.view_init(azim = azim, elev = elev)
	scalarMap._A = []
	ax.set_title("Mean Growth Rate Step %d"%step)
	cbar_ax = fig.add_axes([0.873, 0.2, 0.04, 0.55])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	clrbar.ax.tick_params(labelsize=20)
	clrbar.set_label("Growth Rate step", fontsize = 30)
	#ax.set_title("Time %d : Magnitude of Strain : Max %.4f   ; Min %.4f"%(step,maxTrace,minTrace), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	if save:
		saveDirectory = "../meanGrowthRate"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		plt.savefig(saveDirectory+r"/meanGrowthRate_eta=%.3f_time=%03d.png"%(eta,step),transparent=True)
		plt.close()
	return






####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
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

parser.add_argument("-t","--target", help = "Target face for faster growth", default = 10, type = float)
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


# For surpressing err
class NullDevice():
    def write(self, s):
        pass



original_stdout = sys.stderr # keep a reference to STDOUT

sys.stderr = NullDevice()  # redirect the real STDOUT

################################################################################
import sys
import os
#################################
cwd = os.getcwd()#getting the current working directory
listdir = sorted(os.listdir('.'))#all the list of directory in cwd
saveDirectory = cwd+"/meanGrowthRate"
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
else:
	listdir.remove("meanGrowthRate")# in case the folder was already there, need to drop it from iterating list
counter = 0
totalfolders = len(listdir)
for folder in listdir:
	# Converting folder name to dictionary
	folderdict = dict(item.split("=") for item in folder.split("_"))
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
	plotMeanGrowthRate(step=endStep,azim = azim, elev = elev,numOfLayer = numOfLayer, eta = float(folderdict['n']),save=True)
	os.chdir("..")
	counter+= 1
	plt.close('all')
################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

