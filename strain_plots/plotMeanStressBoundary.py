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
# for a face, projecting its stress eigendecomposed vectors onto radial-orthoradial direction
###############################################################################################################
def getRadialOrthoradialStress(face, radialDict,orthoradialDict, vectors = False):
	eigenvec1 = face.getStressEigenVector1()
	eigenvec2 = face.getStressEigenVector2()
	eigenvalue1 = face.getStressEigenValue1()
	eigenvalue2 = face.getStressEigenValue2()
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
def plotStressAgainstFeedbackPoint(cell,targetid,eta,plot,color='r',large = False,otherplot=None,savefig = None):
	################################################
	cell.calculateStressStrain()
	################################################
	faceList = sf.getPrimordiaBoundaryFaceList(cell,targetid,large= large)
	radialDict, orthoradialDict = getRadialOrthoradialDict(cell,targetid,large = large)
	maximalStress = []
	minimalStress = []
	traceStress = []
	absSumStress = []
	detStress = []
	radialStress = []
	orthoradialStress = []
	sumRadialOrthoradial = []
	sumAbsRadialOrthoradial = []
	######################################################
	for face in faceList:
		maximalStress.append(max(face.getStressEigenValue1(),
								face.getStressEigenValue2()
								)
							)
		minimalStress.append(min(face.getStressEigenValue1(),
								face.getStressEigenValue2()
								)
			)
		absSumStress.append((abs(face.getStressEigenValue1())+
									abs(face.getStressEigenValue2())
									)
							)
		traceStress.append((face.getStressEigenValue1()+face.getStressEigenValue2()))
		detStress.append((qd.getEigenMatrix(face.stress,0,0)*qd.getEigenMatrix(face.stress,1,1)-
						(qd.getEigenMatrix(face.stress,1,0)*qd.getEigenMatrix(face.stress,0,1))
						))
		radstress, orthstress = getRadialOrthoradialStress(face,radialDict,orthoradialDict)
		radialStress.append(radstress)
		orthoradialStress.append(orthstress)
		sumRadialOrthoradial.append(radstress+orthstress)
		sumAbsRadialOrthoradial.append(abs(radstress)+abs(orthstress))
		######################################################
	maximalStress = np.array(maximalStress)
	traceStress = np.array(traceStress)
	absSumStress = np.array(absSumStress)
	detStress = np.array(detStress)
	minimalStress = np.array(minimalStress)
	radialStress = np.array(radialStress)
	orthoradialStress = np.array(orthoradialStress)
	sumRadialOrthoradial = np.array(sumRadialOrthoradial)
	sumAbsRadialOrthoradial = np.array(sumAbsRadialOrthoradial)
	plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
	############################################################
	N = np.sqrt(len(maximalStress))
	if otherplot:
		#print maximalStress,np.mean(maximalStress), np.std(maximalStress)
		plot.errorbar(eta,np.mean(maximalStress), yerr = np.std(maximalStress)/N,fmt='o', color = color,**plotargs)
		otherplot[0].errorbar(eta,np.mean(traceStress), yerr = np.std(traceStress)/N,fmt='o', color = color,**plotargs)
		otherplot[1].errorbar(eta,np.mean(absSumStress), yerr = np.std(absSumStress)/N,fmt='o', color = color,**plotargs)
		otherplot[2].errorbar(eta,np.mean(minimalStress), yerr = np.std(minimalStress)/N,fmt='o', color = color,**plotargs)
		otherplot[3].errorbar(eta,np.mean(sumAbsRadialOrthoradial), yerr = np.std(sumAbsRadialOrthoradial)/N,fmt='o', color = color,**plotargs)
		otherplot[4].errorbar(eta,np.mean(radialStress), yerr = np.std(radialStress)/N,fmt='o', color = color,**plotargs)
		otherplot[5].errorbar(eta,np.mean(orthoradialStress), yerr = np.std(orthoradialStress)/N,fmt='o', color = color,**plotargs)
		otherplot[6].errorbar(eta,np.mean(sumRadialOrthoradial), yerr = np.std(sumRadialOrthoradial)/N,fmt='o', color = color,**plotargs)
		#otherplot[3].errorbar(eta,np.mean(absSumStress), yerr = np.std(absSumStress)/np.sqrt(len(maximalStress)),fmt='o', color = color)
		################################
		# saveplot
		################################
		savefig.errorbar(eta,np.mean(radialStress), yerr = np.std(radialStress)/N,fmt='o', color = 'm',label=r'\sigma_r',**plotargs)
		savefig.errorbar(eta,np.mean(orthoradialStress), yerr = np.std(orthoradialStress)/N,fmt='<', color = 'g',label=r'\sigma_o',**plotargs)
	else:
		plot.errorbar(eta,np.mean(maximalStress), yerr = np.std(maximalStress)/np.sqrt(len(maximalStress)),fmt='o', color = color,**plotargs)
	
	################################################
	return zip([np.mean(radialStress), np.mean(orthoradialStress), np.mean(sumAbsRadialOrthoradial), np.mean(absSumStress)],
			[np.std(radialStress)/N, np.std(orthoradialStress)/N, np.std(sumAbsRadialOrthoradial)/N, np.std(absSumStress)/N])

###############################################################################################################
###	Plotting the Stress magnitude vs feedback
###############################################################################################################
def plotStressAgainstFeedback(targetid, targetHeight, targetArea, eta,endStep, 
	areaplot, heightplot,large = False,otherplot = None,savefig = None,stepsize = 20):
	####################################################
	heightPlotStatus = False
	areaPlotStatus = True
	####################################################
	radialStressData = []
	orthoradialStressData = []
	sumAbsRadialOrthoradialData = []
	absSumStressData = []
	####################################################
	for step in range(1, endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return
		################################################
		cell = sf.loadCellFromFile(step)
		################################################
		if heightPlotStatus:
			primordialHeight = calculatePrimiordiaHeight(cell, targetid, large = large)
			if (primordialHeight > targetHeight):
				for calstep in range(step-1,step-11,-1):
					cell = sf.loadCellFromFile(calstep)
					primordialHeight = calculatePrimiordiaHeight(cell, targetid, large = large)
					if (primordialHeight<=targetHeight):
						heightStressPoints = plotStressAgainstFeedbackPoint(cell,targetid,eta,heightplot,color='salmon' ,large = large,savefig = savefig)
						heightPlotStatus = False
						break
		################################################
		if (areaPlotStatus):
			tissueSurfaceArea = sf.getSurfaceArea(cell)
			if (tissueSurfaceArea > targetArea):
				for calstep in range(step-1,step-stepsize-1,-1):
					cell = sf.loadCellFromFile(calstep)
					tissueSurfaceArea = sf.getSurfaceArea(cell)
					if (tissueSurfaceArea <= targetArea):
						areaStressPoints = plotStressAgainstFeedbackPoint(cell,targetid,eta,areaplot,color='rebeccapurple' ,large = large,otherplot = otherplot,savefig = savefig)
						radialStressData.append(areaStressPoints[0])
						orthoradialStressData.append(areaStressPoints[1])
						sumAbsRadialOrthoradialData.append(areaStressPoints[2])
						absSumStressData.append(areaStressPoints[3])
						areaPlotStatus = False
						break
		################################################
		if not (heightPlotStatus or areaPlotStatus):
			return
		################################################
		gc.collect()

	return [radialStressData, orthoradialStressData,sumAbsRadialOrthoradialData,absSumStressData]
################################################################################################
#        Plotting Part 
################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-c',"--surfaceArea", help="The surface area for which the stress vs feedback plot would need to be plotStrainDifferenceSurface",
					default =0., type = float)
parser.add_argument('-i',"--height", help="The height of primiordia for which the stress vs feedback plot would need to be plotStrainDifferenceSurface",
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
stepsize = args.stepsize

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
fig = plt.figure(figsize=(18,18))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
areaplot = fig.add_subplot(331)
areaplot1 = fig.add_subplot(332)
areaplot2 = fig.add_subplot(333)
areaplot3 = fig.add_subplot(334)

heightplot = fig.add_subplot(335)

areaplot4 = fig.add_subplot(336)
areaplot5 = fig.add_subplot(337)
areaplot6 = fig.add_subplot(338)
areaplot7 = fig.add_subplot(339)




#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Area plot
##################################
areaplot.set_title(r"$A_t =  %d$"%targetArea)
areaplot.set_ylabel(r"$\langle \sigma_{max} \rangle$")
areaplot.set_xlabel(r"$\eta$")

areaplot1.set_title(r"$A_t =  %d$"%targetArea)
areaplot1.set_ylabel(r"$\langle tr (\sigma) \rangle$")
areaplot1.set_xlabel(r"$\eta$")

areaplot2.set_title(r"$A_t =  %d$"%targetArea)
areaplot2.set_ylabel(r"$\langle \sum_i |\sigma_{ii}| \rangle$")
areaplot2.set_xlabel(r"$\eta$")

areaplot3.set_title(r"$A_t =  %d$"%targetArea)
areaplot3.set_ylabel(r"$\langle \sigma_{min} \rangle$")
areaplot3.set_xlabel(r"$\eta$")

areaplot4.set_title(r"$A_t =  %d$"%targetArea)
areaplot4.set_ylabel(r"$\langle |\sigma_o|+|\sigma_r| \rangle$")
areaplot4.set_xlabel(r"$\eta$")

# Radial stress
areaplot5.set_title(r"$A_t =  %d$"%targetArea)
areaplot5.set_ylabel(r"$\langle \sigma_{r} \rangle$")
areaplot5.set_xlabel(r"$\eta$")

# Orthoradial stress
areaplot6.set_title(r"$A_t =  %d$"%targetArea)
areaplot6.set_ylabel(r"$\langle \sigma_{o} \rangle$")
areaplot6.set_xlabel(r"$\eta$")

# sum of the radial and orthoradial stress (should be equal to tr(\sigma))
areaplot7.set_title(r"$A_t =  %d$"%targetArea)
areaplot7.set_ylabel(r"$\langle \sigma_{r}+\sigma_{o} \rangle$")
areaplot7.set_xlabel(r"$\eta$")




#ax1.set_ylim(-0.2,0.)
###################################
# Height of Primodia
###################################
heightplot.set_title(r"$h =  %.2f $"%targetHeight)
heightplot.set_ylabel(r"$\langle \sigma_{max} \rangle$")
heightplot.set_xlabel(r"$\eta$")
###################################

#################################################################################
fig1 = plt.figure(2, figsize=(6,6))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
saveplot = fig1.add_subplot(111)
#########################################
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
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	else:
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	#################################################################################
	if (maxeta != 0) and (etacurrent > maxeta):
		continue
	########################################################
	os.chdir(folder)
	########################################################
	percentStep = int((counter)/float(totalfolders)*100)
	sys.stdout.write('\r'+"step : "+ str(counter) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
	sys.stdout.flush()
	########################################################
	file_name = sorted((fn for fn in os.listdir('.') if fn.startswith('surface')), key = numericalSort)[-1]
	endStep = int(numbers.split(file_name)[1])
	########################################################
	plotData[etacurrent] = plotStressAgainstFeedback(targetface, targetHeight, targetArea,etacurrent, endStep,areaplot=areaplot,
							 heightplot=heightplot,large = large, otherplot = [areaplot1,areaplot2,areaplot3,areaplot4,areaplot5,areaplot6,areaplot7],savefig= saveplot,stepsize = stepsize)
	########################################################
	os.chdir("..")
	counter += 1
########################################################
### Saving figure and Data
########################################################
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	np.save('meanStressBoundary_fk_a=%d_h=%.1f.npy'%(targetArea,targetHeight),plotData)
else:
	np.save('meanStressBoundary_eta_a=%d_h=%.1f.npy'%(targetArea,targetHeight),plotData)
##############################################################################
fig.tight_layout()
if large:# larger primiordia
	fig.savefig(saveDirectory+r"/plotlarge_stressmagnitude_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plotsmall_stressmagnitude_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),transparent = True, bbox_inches="tight")




##########################################################################################
# saveplot saving the orthoradial and radial plot together
##########################################################################################
from collections import OrderedDict


fig1.tight_layout()
handles, labels = saveplot.get_legend_handles_labels()
by_label = OrderedDict(zip(labels,handles))
#print by_label
saveplot.legend(by_label.values(),by_label.keys())


if large:# larger primiordia
	fig1.savefig(saveDirectory+r"/plotlarge_radialOrthoradialStress_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),dpi = 300,transparent = True, bbox_inches="tight")
	fig1.savefig(saveDirectory+r"/plotlarge_radialOrthoradialStress_faceArea=%d_height=%.2f.eps"%(targetArea,targetHeight),format = eps, dpi = 300,transparent = True, bbox_inches="tight")
else:
	fig1.savefig(saveDirectory+r"/plotsmall_radialOrthoradialStress_faceArea=%d_height=%.2f.png"%(targetArea,targetHeight),dpi = 300, transparent = True, bbox_inches="tight")
	fig1.savefig(saveDirectory+r"/plotsmall_radialOrthoradialStress_faceArea=%d_height=%.2f.eps"%(targetArea,targetHeight),format = eps, dpi = 300,transparent = True, bbox_inches="tight")




#fig1.close('all')
#fig.close('all')

################################################################################
print "\n"+"#"*55
print " "*25+"DONE "
print "#"*55

