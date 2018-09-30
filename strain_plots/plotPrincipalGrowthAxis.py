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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
sys.path.append('/home/jkhadka/transferdata/scripts/strain_plots/')
import simulation_functions as sf
import argparse #argument parser, handles the arguments passed by command line
import gc
import ellipse as ep
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

def getCurrentFormMatrix(cell):
	#getting all the target form matrix
	formmatrixDictionary = {}#dictionary to save all the targetformmatrix
	## Getting and saving matrix in array
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	matrixArray = np.zeros((2,2))
	#starting Iteration
	while face != None:
		faceid = face.getID()
		if faceid == 1:
			face = faces.next()
			continue
		matrixArray[0,0] = qd.getCurrentFormMatrix(face, 0,0)
		matrixArray[0,1] = qd.getCurrentFormMatrix(face, 0,1)
		matrixArray[1,0] = qd.getCurrentFormMatrix(face, 1,0)
		matrixArray[1,1] = qd.getCurrentFormMatrix(face, 1,1)
		#saving this in dictionary
		#print matrixArray
		formmatrixDictionary[faceid] = np.copy(matrixArray)
		face = faces.next()
	return formmatrixDictionary
def matrixDifference(formMatrix0,formMatrix):
	diffmatrix = {}
	matrixArray = np.zeros((2,2))
	for key in formMatrix0.keys():
		matrixArray[0,0] = formMatrix[key][0,0]-formMatrix0[key][0,0]
		matrixArray[0,1] = formMatrix[key][0,1]-formMatrix0[key][0,1]
		matrixArray[1,0] = formMatrix[key][1,0]-formMatrix0[key][1,0]
		matrixArray[1,1] = formMatrix[key][1,1]-formMatrix0[key][1,1]
		diffmatrix[key] = np.copy(matrixArray)
	return diffmatrix
################################################################################################
def getMatrixDifferenceEllipsePoints(cell,targetface, diffmatrix,primordialfaceid = 135,zoomfactor=1.):
	#getting the Target Form Matrix
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while True:
		if face.getID() == targetface:
			break
		face = faces.next()
	#print "Face Id : ", face.getID()
	targetformmatrix = diffmatrix
	unitx = face.getUnitx()
	unity = face.getUnity()
	unitz = face.getUnitz()
	unit_mat = np.matrix([[qd.doublearray_getitem(unitx,0),qd.doublearray_getitem(unitx,1),qd.doublearray_getitem(unitx,2)],
						 [qd.doublearray_getitem(unity,0),qd.doublearray_getitem(unity,1),qd.doublearray_getitem(unity,2)],
						 [qd.doublearray_getitem(unitz,0),qd.doublearray_getitem(unitz,1),qd.doublearray_getitem(unitz,2)]])
	#transposing unitmatrix
	transpose_unitmat = np.matrix(np.transpose(unit_mat))
	#Getting Centroid of face
	xcent = face.getXCentralised()
	ycent = face.getYCentralised()
	zcent = face.getZCentralised()
	##### getting data from ellipse & getting transformed coordinate to 3d Cartesian
	data = ep.plot_ellipse(cov=targetformmatrix, data_out=True,norm = True)
	points = zoomfactor*np.matrix(np.vstack((data,np.zeros(len(data[0])))))
	transformedpoints = transpose_unitmat*points
	transformedpoints[0]+= xcent
	transformedpoints[1]+= ycent
	transformedpoints[2]+= zcent
	return transformedpoints


def plotGrowthDirection(cell, diffmatrix, plotdiff,step,targetface, azim = -60, elev = 30, saveDirectory = None):
	# Iterating over Face
	# getting neighbourhood of fast groing cell
	######### Color Map ####################################
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	jet = cm = plt.get_cmap('viridis') 
	maxvalue = 1.#np.pi/2
	minvalue = 0.#-1.*np.pi
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	########################################################
	faceidarray = getNeighbourFaces(cell,targetface)
	primordialface = sf.getFace(cell,targetface)
	px = primordialface.getXCentralised()
	py = primordialface.getYCentralised()
	pz = primordialface.getZCentralised()
	plotdiff.scatter(px,py,pz,marker='*', c= 'm',s = 30)
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1: 
			face = faces.next()
			continue
		faceid = face.getID()
		if faceid in faceidarray:
			c = 'g'
			lw = 3
			zoomfactor = 1
		else:
			c = 'k'
			lw = 1.5
			zoomfactor = 1
		#######################################################
		# Calculation of eigen vecs 
		#######################################################
		facediffmatrix = diffmatrix[face.getID()]
		#########################################
		#Getting points to plot ellipse
		ellipsepoints = getMatrixDifferenceEllipsePoints(cell,face.getID(),facediffmatrix,primordialfaceid = targetface,zoomfactor = zoomfactor)
		###############################
		eigvec, eigval, u = np.linalg.svd(facediffmatrix)
		eigval = np.divide(eigval,np.max(eigval))
		eigval = 10*eigval
		#######################################################
		# Getting the unit matrix for transformation
		#######################################################
		unitx = face.getUnitx()
		unity = face.getUnity()
		unitz = face.getUnitz()
		transformMatrix = np.transpose(np.matrix(
							[[qd.doublearray_getitem(unitx,0),qd.doublearray_getitem(unitx,1),qd.doublearray_getitem(unitx,2)],
							 [qd.doublearray_getitem(unity,0),qd.doublearray_getitem(unity,1),qd.doublearray_getitem(unity,2)],
							 [qd.doublearray_getitem(unitz,0),qd.doublearray_getitem(unitz,1),qd.doublearray_getitem(unitz,2)]]))
		#################################
		#print "#################################"
		#print eigvec
		eigvec = np.matrix(np.vstack((eigvec.T, [0,0])))
		#print eigvec
		eigvec = np.matmul(transformMatrix,eigvec).T
		#print eigvec
		#print transformMatrix
		#################################
		# Now plotting the eigenvecs
		#################################
		xmean = face.getXCentralised()
		ymean = face.getYCentralised()
		zmean = face.getZCentralised()
		plotdiff.scatter(xmean,ymean,zmean,marker='*', c= 'b',s = 30)
		################################################################################################
		# Plotting the angle between the radial and primary direction
		################################################################################################
		################################################################################################
		# Getting the angle for between radial and primary growth direction
		################################################################################################
		#print px, py, pz
		############################################################
		radialvec = np.array([px-xmean,py-ymean,pz-zmean])
		try:
			radialvec /= np.linalg.norm(radialvec)
		except RuntimeWarning:
			print px, py, pz
			print xmean, ymean, zmean
			radialvec = np.array([0.,0.,0.])
		primaryvec = np.array([eigvec[0,0],eigvec[0,1],eigvec[0,2]])
		primaryvec /= np.linalg.norm(primaryvec)
		############################################################
		# angle : here it means dot product
		############################################################
		angle= np.abs(np.dot(radialvec,primaryvec))
		#print px, py, pz, xmean, ymean, zmean, "radialvec ", radialvec
		#print "="*20
		#print faceid, angle
		color = scalarMap.to_rgba(angle)
		#print eigvec, eigvec.shape, eigvec[0],eigval
		xlist = []
		ylist = []
		zlist = []
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
		plotdiff.text(xmean,ymean,zmean, "%.2f"%angle)
		verts = [zip(xlist, ylist,zlist)]
		#plotdiff.text(xmean, ymean, zmean, "%.2f"%(angle))
		#verts = [zip(np.array(ellipsepoints[0])[0],np.array(ellipsepoints[1])[0], np.array(ellipsepoints[2])[0])]
		pc = Poly3DCollection(verts,alpha = .7,linewidths=1, facecolor = color,zorder = 10)
		pc.set_edgecolor(color)
		plotdiff.add_collection3d(pc)
		############################################################
		# For now just plotting Principal direction of Growth
		############################################################
		"""plotdiff.quiver(xmean,ymean,zmean,
													radialvec[0],radialvec[1],radialvec[2],
							color='b',length = 1,pivot='tail',zorder=-1)"""
		if eigval[0]> eigval[1]:
									colorvec1 = 'r'
									vec1 = plotdiff.quiver(xmean,ymean,zmean,
												eigvec[0,0],eigvec[0,1],eigvec[0,2],
												color=colorvec1,length = 1,pivot='tail',zorder=-1)
									#colorvec2 = 'b'
		else:
									colorvec2 = 'r'	
									vec2 =plotdiff.quiver(xmean,ymean,zmean,
												eigvec[1,0],eigvec[1,1],eigvec[1,2],
												color=colorvec2,length = 1,pivot='tail',zorder=-1)
		############################################################
		# Plotting face
		############################################################
		xlist = []
		ylist = []
		zlist = []
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
		if faceid in faceidarray:
			line, = plotdiff.plot(xlist,ylist,zlist,c=c,lw = lw)
		else:
			plotdiff.plot(xlist,ylist,zlist,c=c,lw = lw)
		face = faces.next()
	plotdiff.view_init(azim = azim, elev= elev)
	
	#plotdiff.legend((vec1,vec2,line),("major direction","minor direction","Fast Growing"),bbox_to_anchor=(0.87, 0.95))
	########################################################################
	plt.savefig(saveDirectory+r"/principalGrowth_step=%03d.png"%(step),transparent=True)
	del plotdiff.collections[:]
	plotdiff.lines = []
	return plotdiff
##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotPrincipalGrowthAxis(endStep, targetface, startStep = 0, numOfLayer=8,step = None, alpha = 0.8, Length=1.0,azim = -70, elev=50,saveDirectory = '.'):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	jet = cm = plt.get_cmap('viridis') 
	maxvalue = 1.#np.pi/2
	minvalue = 0.#-1.*np.pi
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	#########################################################
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(10,8))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	scalarMap._A = []
	clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar.set_label("Angle", fontsize = 20)
	########################################################################
	#               Getting the Form Matrix Difference                     #
	########################################################################
	cell0 = sf.loadCellFromFile(startStep) 
	formMatrix0 = getCurrentFormMatrix(cell0)
	for step in range(startStep+1,endStep+1):
		percentStep = int((step-startStep)/float(endStep - startStep)*100)
		sys.stdout.write('\r'+"step : "+ str(step) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
		sys.stdout.flush()
		##################################################################
		cell = sf.loadCellFromFile(step)
		formMatrix = getCurrentFormMatrix(cell)
		## Calculating the difference ##
		diffmatrix = matrixDifference(formMatrix0,formMatrix)
		## Plot the growth direction 
		ax = plotGrowthDirection(cell, diffmatrix, ax, step,targetface=targetface,azim = azim, elev = elev,saveDirectory = saveDirectory)
	plt.close()
	return



	
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =0, type = int)
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

parser.add_argument("-t","--target", help = "Target face for faster growth", default = None, type = int)
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
############################
if targetface == None:
	# get the central face (Top face at the dome)
	targetface = 3*numOfLayer*(numOfLayer-1)+1

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
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]

########################################################
directoryName = r"principalGrowthAxis_eta=%.2f"%ploteta
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
#######################################################
counter = 0
savedict = {}
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	if etacurrent != ploteta : continue
	#print os.listdir('.')
	os.chdir(folder)
	plotPrincipalGrowthAxis(endStep=endStep,startStep = startStep, targetface = targetface, numOfLayer=numOfLayer,step = None, alpha = 0.8, Length=1.0,azim = -70, elev=50,saveDirectory = saveDirectory)
	os.chdir("..")
	gc.collect()
	break

plt.close('all')

################################################################################
print "\n ################################################################################"
print " "*15," DONE "
print "################################################################################"

