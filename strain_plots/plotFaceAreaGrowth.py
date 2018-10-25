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
plt.rcParams['xtick.labelsize'] = 24.
plt.rcParams['ytick.labelsize'] = 24.
plt.rcParams['axes.labelsize'] = 24.
plt.rcParams['legend.fontsize'] = 24.
plt.rcParams['axes.titlesize'] = 30

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
########################################################################
def plotAverageFaceArea(endstep,areaplot,ax2 ,ax3,color,startstep=1,norm=True,fastid = 0):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	# Getting Initial Area
	if not os.path.isfile("qdObject_step=001.obj"):
		return
	cell = sf.loadCellFromFile(1)
	#neighbourhood array 
	faceidarray = getNeighbourFaces(cell,fastid)
	initialarea = {}
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1 : 
			face = faces.next()
			continue
		initialarea[face.getID()] = face.getAreaOfFace()
		face =faces.next()
	######################################################
	# Gathering face area
	######################################################
	fastfaceareamean = []
	fastfaceareastd = []
	slowfaceareamean = []
	slowfaceareastd = []
	totalMeanFaceArea = []
	for i in range(startstep,endstep+1):
		if not os.path.isfile("qdObject_step=%03d.obj"%i):#check if file exists
			break
		cell = sf.loadCellFromFile(i)
		fastfaceareaarray = []
		slowfaceareaarray = []
		######################################
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == 1 : 
				face = faces.next()
				continue
			area0 =initialarea[face.getID()]
			#print face.getID(), area0
			area = face.getAreaOfFace()
			if norm:
				facearea = area/area0
			else:
				facearea = area
			if face.getID() in faceidarray:
				fastfaceareaarray.append(facearea)
			else:
				slowfaceareaarray.append(facearea)
			face =faces.next()
		totalMeanFaceArea.append(np.mean(np.array(slowfaceareaarray + fastfaceareaarray)))
		fastfaceareamean.append(np.mean(fastfaceareaarray))
		fastfaceareastd.append(np.std(fastfaceareaarray))
		slowfaceareamean.append(np.mean(slowfaceareaarray))
		slowfaceareastd.append(np.std(slowfaceareaarray))
	#print facearea
	######################################################
	# making plot
	######################################################
	########################
	# plotting
	########################
	fastfaceareamean = np.array(fastfaceareamean)
	fastfaceareastd = np.array(fastfaceareastd)
	slowfaceareamean = np.array(slowfaceareamean)
	slowfaceareastd = np.array(slowfaceareastd)
	areaplot.fill_between(range(len(fastfaceareamean)),fastfaceareamean-fastfaceareastd,fastfaceareamean+fastfaceareastd,alpha = 0.3,color = color)
	areaplot.plot(range(len(fastfaceareamean)),fastfaceareamean,color = color)
	areaplot.fill_between(range(len(slowfaceareamean)),slowfaceareamean-slowfaceareastd,slowfaceareamean+slowfaceareastd,alpha = 0.5,color = color)
	areaplot.plot(range(len(slowfaceareamean)),slowfaceareamean,color = color)
	ax2.plot(range(len(totalMeanFaceArea)),totalMeanFaceArea,color = color)
	ax3.plot(np.log(range(1,len(totalMeanFaceArea)+1)),np.log(totalMeanFaceArea),color = color)
	#ax2.fill_between(range(len(fastfaceareamean)),fastfaceareamean-fastfaceareastd,fastfaceareamean+fastfaceareastd,alpha = 0.5)
	#ax2.plot(range(len(fastfaceareamean)),fastfaceareamean)
	#ax2.fill_between(range(len(slowfaceareamean)),slowfaceareamean-slowfaceareastd,slowfaceareamean+slowfaceareastd,alpha = 0.5)
	#ax2.plot(range(len(slowfaceareamean)),slowfaceareamean)
	#ax2.plot(range(len(value)),value,c = scalarMap.to_rgba(key))
	#ax1.set_xlim(0,100)
	#plt.show()
	return
########################################################################
def getAreaGrowthData(cell, areaCellDict, surfaceAreaArray,dAreaCellDict,counter):
	surfaceAreaArray[counter] = cell.getSurfaceArea()
	dareaTissue = surfaceAreaArray[counter]-surfaceAreaArray[counter-1]
	################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face 	!= None:
		faceid = face.getID()
		if faceid == 1:
			face = faces.next()
			continue
		########################################
		areaCellDict[faceid][counter] = face.getAreaOfFace() 
		dareaCell = areaCellDict[faceid][counter]-areaCellDict[faceid][counter-1]
		dAreaCellDict[faceid] = dareaCell/dareaTissue
		########################################
		face = faces.next()
	########################################
	return
########################################################################
def plotFaceAreaDerivative(faceAreaDerivativePlot,cell,dAreaCellDict,targetid,colormap = 'cool',alpha = 0.8,
	azim = -60, elev = 50):
	###############################################################
	# Average area growth rate
	###############################################################
	faceidarray = getNeighbourFaces(cell,targetid)
	averagedDArea = {}
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	minMagnitude = 0.
	maxMagnitude = 0.
	while (face != None):
		if ((face.getID()==1) or (face.getID() in faceidarray)):
			face  = faces.next()
			continue
		averagedDArea[face.getID()] = np.mean(dAreaCellDict[face.getID()])
		if averagedDArea[face.getID()] > maxMagnitude:
			maxMagnitude = averagedDArea[face.getID()]
		elif averagedDArea[face.getID()] < minMagnitude:
			minMagnitude = averagedDArea[face.getID()]
		###########################################################
		face = faces.next()
	###############################################################
	#                 Plotting the Cell                          #
	##############################################################
	jet = cm = plt.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	#########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID()==1:
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
		########################################################################################
		if face.getID() in faceidarray:
			color = 'w'
		else:
			color = scalarMap.to_rgba(averagedDArea[face.getID()])
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		faceAreaDerivativePlot.add_collection3d(pc)
		#faceAreaDerivativePlot.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
	########################################################################################
	scalarMap._A = []
	clrbar1 = plt.colorbar(scalarMap, ax=faceAreaDerivativePlot,shrink = 0.5,aspect = 10,ticks=np.linspace(minMagnitude,maxMagnitude,3))
	clrbar1.set_label(r"Area Growth Rate")
	faceAreaDerivativePlot.view_init(azim = azim, elev = elev)
	#######################################################
	return 
########################################################################
def plotAverageGrowthRate(endStep,areaDerivativePlot, faceAreaDerivativePlot,targetid, startStep=1,norm=True,
	fastid = 0,azim = -60, 
	elev = 50,stepsize = 1):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	######################################################
	# Getting the first step 
	######################################################
	if not os.path.isfile("qdObject_step=%03d.obj"%startStep):
		return
	cell = sf.loadCellFromFile(startStep)
	######################################################
	# dict of area
	######################################################
	dAreaCellDict = {}
	areaCellDict= {}
	surfaceAreaArray = np.zeros(endStep-startStep)
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	surfaceAreaArray[0] = cell.getSurfaceArea()
	######################################################
	while face != None:
		if face.getID() == 1 : 
			face = faces.next()
			continue
		areaCellDict[face.getID()] = np.zeros(int((endStep-startStep)/stepsize)+1)
		dAreaCellDict[face.getID()] = np.zeros(int((endStep-startStep)/stepsize))
		areaCellDict[face.getID()][0] = face.getAreaOfFace() 
		face =faces.next()
	######################################################
	# Gathering face area
	######################################################
	stepcounter = 0
	########################################
	for i in range(startStep+1,endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%i):#check if file exists
			break
		cell = sf.loadCellFromFile(i)
		######################################
		getAreaGrowthData(cell, areaCellDict, surfaceAreaArray,dAreaCellDict,stepcounter)
		######################################################
		stepcounter += 1
		######################################################
	########################
	# plotting
	########################
	plotFaceAreaDerivative(faceAreaDerivativePlot,cell,dAreaCellDict,targetid,azim = azim, 
		elev = elev)
	########################
	return
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
						return calstep,tissueSurfaceArea
		################################################
		gc.collect()
	return endStep,tissueSurfaceArea
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-d","--stepsize", help = "stepsize on plot", default = 1, type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
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
saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)
surfaceSaveDirectory = saveDirectory+r"/areaGrowthRate"
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
if targetface == None:
    if numOfLayer == 8:
        targetface = 135
    elif numOfLayer == 10:
        targetface = 214


Length = 1.
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(111,projection='3d')
ax2.set_xlim((-0.6*radius,0.6*radius))
ax2.set_ylim((-0.6*radius,0.6*radius))
ax2.set_zlim((-0.4*radius,0.8*radius))
ax2.axis('off')
#################################################################################
if targetarea:
	endStep,surfacearea = getTimeStep(targetArea, endStep, startStep=startStep, stepsize = stepsize)
	ax2.set_title("Surface Area = %d"%surfacearea)

plotAverageGrowthRate(endStep,areaDerivativePlot=None, 
	faceAreaDerivativePlot=ax2,targetid = targetface,
	startStep=startStep,
	norm=True,fastid = 0,azim = azim, elev = elev,stepsize = stepsize)
################################################################################
fig.savefig(surfaceSaveDirectory+'/'+'faceAreaGrowthRate_time=%03d.png'%endStep,
	bbox_inches='tight',dpi=100,transparent = True)
################################################################################
print "\n"+ "#"*100
print " "*45+"DONE "
print "#"*100
################################################################################
plt.close('all')