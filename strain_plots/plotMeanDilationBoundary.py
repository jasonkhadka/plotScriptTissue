import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/jkhadka/plantdev")
import quadedge as qd
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
import gc

#sys.path.append("/home/jkhadka/transferfile/scripts/simulation_functions")
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
import simulation_functions as sf
plt.rcParams['xtick.labelsize'] = 20.
plt.rcParams['ytick.labelsize'] = 20.
plt.rcParams['axes.labelsize'] = 25.
plt.rcParams['legend.fontsize'] = 25.
plt.rcParams['axes.titlesize'] = 30.

import re
#from PIL import Image
#############################################################
#key to sort the file_names in order
numbers = re.compile(r'(\d+)')
def numericalSort(value):
	parts = numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts
############################################################
###############################################################################################################
# Function to calculate the height of primordia
###############################################################################################################
def calculatePrimiordiaHeight(cell,targetid,large = False):
	# getting the vertexlist of primiordia boundary
	vertexlist = sf.getPrimordiaBoundaryVertexList(cell,targetid, large = large)
	vertexNum = len(vertexList)
	meanx, meany, meanz = 0.,0.,0.
	for vertex in vertexlist:
		meanx += vertex.getXcoordinate()
		meany += vertex.getYcoordinate()
		meanz += vertex.getZcoordinate()
	meanx /= vertexNum
	meany /= vertexNum
	meanz /= vertexNum
	##########################################################################
	targetface = sf.getFace(cell,targetid)
	#####################################
	targetx = facetarget.getXCentralised()
	targety = facetarget.getYCentralised()
	targetz = facetarget.getZCentralised()
	#####################################
	height = np.sqrt((meanx-targetx)**2+(meany-targety)**2+(meanz-targetz)**2)
	#####################################
	return height
###############################################################################################################
# Function to calculate the height of primordia
###############################################################################################################
def calculatePrimiordiaHeight(cell,targetid,large = False):
	# getting the vertexlist of primiordia boundary
	vertexlist = sf.getPrimordiaBoundaryVertexList(cell,targetid)#, large = large)
	vertexNum = len(vertexlist)
	meanx, meany, meanz = 0.,0.,0.
	for vertex in vertexlist:
		meanx += vertex.getXcoordinate()
		meany += vertex.getYcoordinate()
		meanz += vertex.getZcoordinate()
	meanx /= vertexNum
	meany /= vertexNum
	meanz /= vertexNum
	##########################################################################
	facetarget = sf.getFace(cell,targetid)
	#####################################
	targetx = facetarget.getXCentralised()
	targety = facetarget.getYCentralised()
	targetz = facetarget.getZCentralised()
	#####################################
	height = np.sqrt((meanx-targetx)**2+(meany-targety)**2+(meanz-targetz)**2)
	#####################################
	return height
###############################################################################################################
# Get 2 dictionaries one with radialDirection vectors for each face and 
# another with orthoradialDirection vectors for each face
###############################################################################################################
def getRadialOrthoradialDict(cell,targetid, large = False):
	primordiafacelist= sf.getSeparatePrimordiaBoundaryFaceList(cell, targetid, large=large)
	orthoradialDict = {}
	radialDict = {}
	targetface = sf.getFace(cell,targetid)
	targetcentroid = np.array([targetface.getXCentralised(), targetface.getYCentralised(),targetface.getZCentralised()])
	############################################################
	for listiter in range(len(primordiafacelist)):
		facelist = primordiafacelist[listiter]
		########################################################
		for count in range(len(facelist)):
			face = facelist[count]
			#nextface = facelist[(count+1)%len(facelist)]
			normal = face.getNormal()
			normalvec = np.array([qd.doublearray_getitem(normal,0),
							qd.doublearray_getitem(normal,1),
							qd.doublearray_getitem(normal,2)])
			########################################################
			facecentroid = np.array([face.getXCentralised(), face.getYCentralised(),face.getZCentralised()])
			#Calculating direction Vector to primordia centroid
			direcvec = np.subtract(targetcentroid,facecentroid)
			direcvec = direcvec/np.linalg.norm(direcvec)
			#Radial Vec to on-plane unitvector
			radialvec = direcvec - normalvec*(np.dot(direcvec,normalvec))
			radialvec = radialvec/np.linalg.norm(radialvec)
			########################################################
			crossvec = np.cross(radialvec,normalvec)
			crossvec = crossvec/np.linalg.norm(crossvec)
			orthoradialDict[face.getID()] = crossvec
			radialDict[face.getID()]= radialvec
	return radialDict,orthoradialDict
###############################################################################################################
# for a face, projecting its Growth eigendecomposed vectors onto radial-orthoradial direction
###############################################################################################################
def getRadialOrthoradialGrowth(face, radialDict,orthoradialDict, vectors = False):
	eigenvec1 = face.getRotGrowthEigenVector1()
	eigenvec2 = face.getRotGrowthEigenVector2()
	eigenvalue1 = face.getRotGrowthEigenValue1()
	eigenvalue2 = face.getRotGrowthEigenValue2()
	vec1 =eigenvalue1*np.array([qd.doublearray_getitem(eigenvec1,0),
					qd.doublearray_getitem(eigenvec1,1),
					qd.doublearray_getitem(eigenvec1,2)])
	vec2 = eigenvalue2*np.array([qd.doublearray_getitem(eigenvec2,0),
					qd.doublearray_getitem(eigenvec2,1),
					qd.doublearray_getitem(eigenvec2,2)])
	radialvec = np.copy(radialDict[face.getID()])
	orthoradialvec = np.copy(orthoradialDict[face.getID()])
	#print radialvec,orthoradialvec
	############################################
	radialComp = np.dot(radialvec,vec1)+np.dot(radialvec,vec2)
	orthoradialComp = np.dot(orthoradialvec,vec1)+np.dot(orthoradialvec,vec2)
	############################################
	if vectors:
		radialvec = radialComp*radialvec
		orthoradialvec = orthoradialComp*orthoradialvec
		return radialvec, orthoradialvec
	############################################
	return radialComp, orthoradialComp#radialvec,orthoradialvec
###############################################################################################################
# Function to calculate the height of primordia
###############################################################################################################
def plotGrowthAgainstFeedbackPoint(cell1,cell2, targetid,eta,plot,color='r',large = False,otherplot=None):
	################################################
	sf.calculateDilation(cell1,cell2)
	################################################
	faceList = sf.getPrimordiaBoundaryFaceList(cell1,targetid,large= large)
	radialDict, orthoradialDict = getRadialOrthoradialDict(cell1,targetid,large = large)
	maximalGrowth = []
	minimalGrowth = []
	traceGrowth = []
	absSumGrowth = []
	detGrowth = []
	radialGrowth = []
	orthoradialGrowth = []
	sumRadialOrthoradial = []
	sumAbsRadialOrthoradial = []
	######################################################
	for face in faceList:
		maximalGrowth.append(max(face.getRotGrowthEigenValue1(),
								face.getRotGrowthEigenValue2()
								)
							)
		minimalGrowth.append(min(face.getRotGrowthEigenValue1(),
								face.getRotGrowthEigenValue2()
								)
			)
		absSumGrowth.append((abs(face.getRotGrowthEigenValue1())+
									abs(face.getRotGrowthEigenValue2())
									)
							)
		traceGrowth.append((face.getRotGrowthEigenValue1()+face.getRotGrowthEigenValue2()))
		detGrowth.append((qd.getEigenMatrix(face.rotGrowth,0,0)*qd.getEigenMatrix(face.rotGrowth,1,1)-
						(qd.getEigenMatrix(face.rotGrowth,1,0)*qd.getEigenMatrix(face.rotGrowth,0,1))
						))
		radGrowth, orthGrowth = getRadialOrthoradialGrowth(face,radialDict,orthoradialDict)
		radialGrowth.append(radGrowth)
		orthoradialGrowth.append(orthGrowth)
		sumRadialOrthoradial.append(radGrowth+orthGrowth)
		sumAbsRadialOrthoradial.append(abs(radGrowth)+abs(orthGrowth))
		######################################################
	maximalGrowth = np.array(maximalGrowth)
	traceGrowth = np.array(traceGrowth)
	absSumGrowth = np.array(absSumGrowth)
	detGrowth = np.array(detGrowth)
	minimalGrowth = np.array(minimalGrowth)
	radialGrowth = np.array(radialGrowth)
	orthoradialGrowth = np.array(orthoradialGrowth)
	sumRadialOrthoradial = np.array(sumRadialOrthoradial)
	sumAbsRadialOrthoradial = np.array(sumAbsRadialOrthoradial)
	plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
	############################################################
	N = np.sqrt(len(maximalGrowth))
	if otherplot:
		#print maximalGrowth,np.mean(maximalGrowth), np.std(maximalGrowth)
		plot.errorbar(eta,np.mean(maximalGrowth), yerr = np.std(maximalGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[0].errorbar(eta,np.mean(traceGrowth), yerr = np.std(traceGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[1].errorbar(eta,np.mean(absSumGrowth), yerr = np.std(absSumGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[2].errorbar(eta,np.mean(minimalGrowth), yerr = np.std(minimalGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[3].errorbar(eta,np.mean(sumAbsRadialOrthoradial), yerr = np.std(sumAbsRadialOrthoradial)/N,fmt='o', color = color,**plotargs)
		otherplot[4].errorbar(eta,np.mean(radialGrowth), yerr = np.std(radialGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[5].errorbar(eta,np.mean(orthoradialGrowth), yerr = np.std(orthoradialGrowth)/N,fmt='o', color = color,**plotargs)
		otherplot[6].errorbar(eta,np.mean(sumRadialOrthoradial), yerr = np.std(sumRadialOrthoradial)/N,fmt='o', color = color,**plotargs)
		#otherplot[3].errorbar(eta,np.mean(absSumGrowth), yerr = np.std(absSumGrowth)/np.sqrt(len(maximalGrowth)),fmt='o', color = color)
	elif areaplot:
		plot.errorbar(eta,np.mean(maximalGrowth), yerr = np.std(maximalGrowth)/np.sqrt(len(maximalGrowth)),fmt='o', color = color,**plotargs)
	
	################################################
	return zip([np.mean(radialGrowth), np.mean(orthoradialGrowth), np.mean(sumAbsRadialOrthoradial), np.mean(absSumGrowth)],
			[np.std(radialGrowth)/N, np.std(orthoradialGrowth)/N, np.std(sumAbsRadialOrthoradial)/N, np.std(absSumGrowth)/N])

###############################################################################################################
###	Plotting the Growth magnitude vs feedback
###############################################################################################################
def plotGrowthAgainstFeedback(targetid, targetHeight, targetArea, eta,endStep, 
	areaplot, heightplot,large = False,otherplot = None,stepsize = 5):
	####################################################
	heightPlotStatus = False
	areaPlotStatus = True
	####################################################
	radialGrowthData = []
	orthoradialGrowthData = []
	sumAbsRadialOrthoradialData = []
	absSumGrowthData = []
	####################################################
	for step in range(1, endStep,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return
		################################################
		percentStep = int((step)/endStep*100)
		sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
		sys.stdout.flush()
		################################################
		cell1 = sf.loadCellFromFile(step)
		cell2 = sf.loadCellFromFile(step+stepsize)
		areaGrowthPoints = plotGrowthAgainstFeedbackPoint(cell1,cell2,targetid,eta,None,color='rebeccapurple' ,large = large,otherplot = None)
		radialGrowthData.append(areaGrowthPoints[0])
		orthoradialGrowthData.append(areaGrowthPoints[1])
		sumAbsRadialOrthoradialData.append(areaGrowthPoints[2])
		absSumGrowthData.append(areaGrowthPoints[3])
		################################################
		gc.collect()
	return [radialGrowthData, orthoradialGrowthData,sumAbsRadialOrthoradialData,absSumGrowthData]
################################################################################################
#        Plotting Part 
################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-c',"--surfaceArea", help="The surface area for which the Growth vs feedback plot would need to be plotStrainDifferenceSurface",
					default =0., type = float)
parser.add_argument('-i',"--height", help="The height of primiordia for which the Growth vs feedback plot would need to be plotStrainDifferenceSurface",
					default =0., type = float)
parser.add_argument("-m","--maxeta", help = "if this is given, then eta is only cacluated till this value", type = float, default = 0.0)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument("-f","--fastkappa", help = "if option is used, the figures are made with respect to chaning fast kappa", action= "store_true")
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument('-d',"--stepsize", help="stepsize for progressing on looking for match of facearea and height", type = int,default = 20)
#parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",default = 10, type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
											  default = 1., type = float)
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

parser.add_argument("-t","--target", help = "Target face for faster growth", default = None, type = float)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)
parser.add_argument("-L","--Large", help = "if option is used, calculation is done for larger primordia", action= "store_true")
## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
#cylinder = args.cylinder
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
anglethreshold = args.angle
fastkappaOption = args.fastkappa
maxeta = args.maxeta
stepsize = stepsize
#############################################
azim = args.azimuthal
elev = args.elevation
norm = args.nonNormalize

targetArea = args.surfaceArea
targetHeight = args.height
targetface = args.target
large = args.Large


if targetface == None:
	if numOfLayer == 8:
		targetface = 135
	elif numOfLayer == 10:
		targetface = 214

#################################################################################
fig = plt.figure(figsize=(12,6))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
plot1 = fig.add_subplot(121,,projection='3d')
plot2 = fig.add_subplot(122,projection='3d')

#################################################################################
# Area plot
##################################

import sys
import os
cwd = os.getcwd()#getting the current working directory
DIRNAMES=1
listdir = sorted(os.walk('.').next()[DIRNAMES])#all the list of directory in cwd
directoryName = "plots"
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
totalfolders = len(listdir)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
# if true calculate with respect to changing fastkappa, else Eta
if fastkappaOption:
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['fk']) for folder in listdir]
else:
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
#################################################################################
counter = 0
plotData = {}
########################################################
plotMeanDilation(endStep,plot1 = plot1, plot2 = plot2,large = large,stepsize = stepsize)
########################################################
### Saving figure and Data
########################################################
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	np.save('meanDilationBoundary_fk_a=%d_h=%.1f.npy'%(targetArea,targetHeight),plotData)
else:
	np.save('meanDilationBoundary_eta_a=%d_h=%.1f.npy'%(targetArea,targetHeight),plotData)
##############################################################################
plt.tight_layout()
if large:# larger primiordia
	fig.savefig(saveDirectory+r"/plotlarge_dilationmagnitude_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plotsmall_dilationmagnitude_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),transparent = True, bbox_inches="tight")

plt.close('all')


################################################################################
print "\n"+"#"*55
print " "*25+"DONE "
print "#"*55

