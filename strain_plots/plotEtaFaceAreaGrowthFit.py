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
def plotAverageFaceArea(endstep,areaplot,ax2 ,ax3,m0plot,m0log,color,startstep=1,norm=False,fastid = 0):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	# Getting Initial Area
	######################################################
	# Gathering face area
	######################################################
	totalMeanFaceArea = []
	m0determinantArray = []
	for i in range(startstep,endstep+1):
		if not os.path.isfile("qdObject_step=%03d.obj"%i):#check if file exists
			break
		cell = sf.loadCellFromFile(i)
		cell.setParameters()
		facearea = 0.
		tfmDet = 0.
		######################################
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		numfaces = 0.
		while face != None:
			if face.getID() == 1 : 
				face = faces.next()
				continue
			tfmDet += face.getTargetFormMatrixDeterminant()
			#print face.getID(), area0
			facearea+= face.getAreaOfFace()
			if np.isnan(face.getAreaOfFace()):
				print "Step :", i, "face ID : ", face.getID()
			numfaces += 1.
			face  = faces.next()
		totalMeanFaceArea.append(facearea/numfaces)
		m0determinantArray.append(tfmDet/numfaces)
	#print facearea
	######################################################
	# making plot
	######################################################
	########################
	# plotting
	########################
	#areaplot.fill_between(range(len(fastfaceareamean)),fastfaceareamean-fastfaceareastd,fastfaceareamean+fastfaceareastd,alpha = 0.3,color = color)
	#areaplot.plot(range(len(fastfaceareamean)),fastfaceareamean,color = color)
	#areaplot.fill_between(range(len(slowfaceareamean)),slowfaceareamean-slowfaceareastd,slowfaceareamean+slowfaceareastd,alpha = 0.5,color = color)
	#areaplot.plot(range(len(slowfaceareamean)),slowfaceareamean,color = color)
	time = np.arange(startstep, len(totalMeanFaceArea)+startstep)
	ax2.plot(time,totalMeanFaceArea,color = color)
	import math
	logMeanArea = [math.log(i) for i in totalMeanFaceArea]
	ax3.plot(time,logMeanArea,color = color)
	m0plot.plot(time,m0determinantArray,color = color)
	logm0 = []
	for i in m0determinantArray:
		try:
			math.log(i)
		except ValueError:
			print os.getcwd(), 7*" ",i
			logm0.append(0)
		else:
			logm0.append(math.log(i))
	m0log.plot(time,logm0,color = color)
	#ax2.fill_between(range(len(fastfaceareamean)),fastfaceareamean-fastfaceareastd,fastfaceareamean+fastfaceareastd,alpha = 0.5)
	#ax2.plot(range(len(fastfaceareamean)),fastfaceareamean)
	#ax2.fill_between(range(len(slowfaceareamean)),slowfaceareamean-slowfaceareastd,slowfaceareamean+slowfaceareastd,alpha = 0.5)
	#ax2.plot(range(len(slowfaceareamean)),slowfaceareamean)
	#ax2.plot(range(len(value)),value,c = scalarMap.to_rgba(key))
	#ax1.set_xlim(0,100)
	#plt.show()
	return [totalMeanFaceArea, m0determinantArray]
	
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
directoryName = "gaussianPlot"
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
#################################################################################
#     Making the Plots
#################################################################################
import matplotlib.colors as colors
import matplotlib.cm as cmx
jet = cm = plt.get_cmap('viridis') 
maxvalue = max(etalist)
minvalue = min(etalist)
#print "Max value", maxvalue, " minvalue", minvalue
cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#fig = plt.figure(frameon=False,figsize=(20,16))
fig = plt.figure(1,figsize=(10,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
#ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(211)
ax3 = fig.add_subplot(212)
#ax2 = fig.add_subplot(122)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
#ax1.set_title("Face Area")
#ax1.set_xlabel("Time")
#ax1.set_ylabel(r"$<A>$")
##################################
ax2.set_title("Mean Face Area")
ax2.set_xlabel("Time")
ax2.set_ylabel(r"$<A>$")
##################################
ax3.set_title("$log(<A>)$ vs time")
ax3.set_xlabel(r"$t$")
ax3.set_ylabel(r"$log(<A>)$")
#ax1.set_ylim(0.,80.)
fig1 = plt.figure(figsize=(10,10))
m0plot = fig1.add_subplot(211)
m0plot.set_title("Determinant of TargetFormmatrix")
m0plot.set_xlabel("Determinant")
m0plot.set_ylabel("time")
m0log = fig1.add_subplot(212)
m0log.set_title("$log(determinant)$ of TargetFormmatrix")
m0log.set_ylabel("$log(determinant)$")
m0log.set_xlabel("time")

#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
###################################
# Height of Primodia
###################################
#ax2.set_title("Height of Primodia")
#ax2.set_xlabel("time")
#ax2.set_ylabel("Height of Primodia")
########################################################
counter = 0
totalfolders = len(listdir)
#print listdir
savedict = {}
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	if etacurrent != 0:
		continue
	#scaled color for plotting
	#print fkcurrent
	etacolor = scalarMap.to_rgba(etacurrent)
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
	
	savedict[etacurrent] = plotAverageFaceArea(endStep,0,ax2,ax3,m0plot, m0log,color = etacolor,startstep=1,fastid = targetface)
	#plotMeanGrowthRate(step=endStep,azim = azim, elev = elev,eta = float(folderdict['n']),save=True)
	plotEnergy(endStep,etacurrent)
	#os.chdir("..")
	gc.collect()
	counter+= 1
#################################################################################
# fitting exponential function to the growth
#################################################################################
def fitExpFunc(t,A0,l):
    return A0*np.exp(l*t)
def fitLinFunc(t,m,c):
    return m*t + c
import scipy.optimize as sop


#######################################################
# Array shape : 
# [totalMeanFaceArea, m0determinantArray]
#######################################################
noFeedbackAreaGrowth = savedict[0][0]
time = np.arange(len(noFeedbackAreaGrowth))
logNoFeedback = np.log(noFeedbackAreaGrowth)
ppot, pcov = sop.curve_fit(fitLinFunc,time, logNoFeedback)
print noFeedbackAreaGrowth, np.array(noFeedbackAreaGrowth).shape
#Plotting 
fitfig = plt.figure(3,figsize=(10,20))
fitfig.suptitle("No Feedback Exponential Plot",fontsize = 40)
fitplot1 = fitfig.add_subplot(211)
fitplot1.set_title("Area vs Time")
fitplot1.set_xlabel("Time")
fitplot1.set_ylabel("Area")
fitplot2 = fitfig.add_subplot(212)
fitplot2.set_title("Area vs Time")
fitplot2.set_xlabel("Time")
fitplot2.set_ylabel("log(Area)")

from matplotlib.offsetbox import AnchoredText
at = AnchoredText(r"$A = %.5fe^{%.7ft}$"%(np.exp(ppot[1]),ppot[0]),loc = 2,prop=dict(size=30))
fitplot1.add_artist(at)
fitplot1.plot(time, noFeedbackAreaGrowth,'g',lw= 3)
fitplot1.plot(time, np.exp(fitLinFunc(time,ppot[0],ppot[1])),'r:',lw= 3)

#ax2.plot(time, fitfunc(time,ppot[0],ppot[1]),'k:')
fitplot2.plot(time, logNoFeedback,lw= 3)
fitplot2.plot(time, fitLinFunc(time,ppot[0],ppot[1]),'k:',lw = 3)


fitfig.tight_layout(rect=[0, 0.03, 1, 0.95])
fitfig.savefig(saveDirectory+r"/plot_eta_face_area_fit_exp_function.png",transparent = True)
plt.close(3)


###############################################################################
### Saving Data Dictionary ###
np.save('facegrowth_m0growth.npy',savedict)

#color bar
scalarMap._A = []
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.87, 0.1, 0.02, 0.7])
clrbar = plt.colorbar(scalarMap,cax = cbar_ax)
clrbar.set_label("$\eta$")
#plt.tight_layout()
#plt.tight_layout( rect=[0, 0, 1, 1])

fig.tight_layout()
fig.savefig(saveDirectory+r"/plot_eta_face_area.png",transparent = True)

#color bar
scalarMap._A = []
fig1.subplots_adjust(bottom=0.2)
cbar_ax = fig1.add_axes([0.87, 0.1, 0.02, 0.7])
clrbar = plt.colorbar(scalarMap,cax = cbar_ax)
clrbar.set_label("$\eta$")
#plt.tight_layout()
#plt.tight_layout( rect=[0, 0, 1, 1])
fig1.tight_layout()
fig1.savefig(saveDirectory+r"/plot_eta_targetform_determinant.png",transparent = True)
plt.close('all')

################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

