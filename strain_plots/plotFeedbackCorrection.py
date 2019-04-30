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
plt.rcParams['xtick.labelsize'] = 18.
plt.rcParams['ytick.labelsize'] = 18.
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.titlesize'] = 18

####################################################################################################################
# calculating the roundness of the faces
####################################################################################################################
def calculateMeanRoundness(facelist):
	totalRoundness = 0.
	for face in facelist:
		#########################
		#getting perimeter
		#########################
		face.setSumEdgeLength()
		perimeter = face.getSumEdgeLength()
		areaofface = face.getAreaOfFace()
		#########################
		totalRoundness += 2.*np.sqrt(np.pi*areaofface)/perimeter
		#########################
	return totalRoundness/len(facelist)

####################################################################################################################
# Calculating and plotting mean stress and growth
####################################################################################################################
def plotFeedbackCorrection(targetid, targetsurfacearea,endStep = 2000, 
	startStep=1,stepsize= 20,maxarea = None, areastep = 10,resetids = False,
	saveName = None):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	########################################################################
	# faceidarray for Primordia
	if not os.path.isfile("qdObject_step=001.obj"):
		return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
	########################################################################
	# getting the time step for computation
	########################################################################
	step,tissueSurfaceArea = sf.getTimeStep(targetsurfacearea, endStep, startStep, stepsize = stepsize,resetids = resetids)
	#######################################################################
	# Starting the Calculation
	#######################################################################
	cell = sf.loadCellFromFile(step,resetids = resetids)#loading cell
	#######################################################################
	targetface = sf.getFace(cell, targetid)#getting top priomrdial face
	cell.calculateStressStrain()#calculating stress and strain
	cell.setRadialOrthoradialVector(targetface)#setting the rad/orthorad vectors
	#######################################################################
	# calculating the rad/orthorad Component of Correction terms
	#######################################################################
	cell.setRadialOrthoradialFeedbackCorrection()
	cell.setRadialOrthoradialStress()
	cell.setPrincipalDeformationVector()
	cell.setPrincipalDeformationFeedbackCorrection()
	#######################################################################
	# Plotting
	#######################################################################
	sf.plotRadialOrthoradialFeedbackCorrectionSurface(cell, numOfLayer,
		azim = -60, elev = 60,
		name =saveName+"_feedbackcorrection_radialOrthoradial")
	sf.plotPrincipalDeformationFeedbackCorrectionSurface(cell, numOfLayer,
		azim = -60, elev = 60,
		name =saveName+"_feedbackcorrection_principalDeformation")
	sf.plotStressSurface(cell, numOfLayer,
		azim = -60, elev = 60,
		name =saveName+"_stressSurface_principal")
	return
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int,default = 2000)
parser.add_argument("-m","--maxeta", help = "if this is given, then eta is only cacluated till this value", type = float, default = 0.0)
parser.add_argument("-x","--maxarea", help = "if this is given, then plot is only made till this area value value", type = float, default = None)
parser.add_argument('-c',"--surfaceArea", help="The surface area for which the stress vs feedback plot would need to be plotStrainDifferenceSurface",
					default =850, type = int)
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int,default = 8)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
											  default = 1., type = float)
parser.add_argument("-n","--nonNormalize", help = "if option is used, the figures are not normalised", action= "store_false")
parser.add_argument("-L","--Large", help = "if option is used, calculation is done for larger primordia", action= "store_true")
parser.add_argument("-f","--fastkappa", help = "if option is used, the figures are made with respect to chaning fast kappa", action= "store_true")
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 1., type = float)

parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = None, type = int)
parser.add_argument("-o","--othersidetissue",help = "target cell for other side of the tissue (non-primordia)", type = int,default = 226)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)
parser.add_argument('-d',"--areastep", help="area step for calculating the growth in cell area", type = int,
						default = 20)
parser.add_argument('-j',"--jobid", help="jobid", type = int,default = None)
parser.add_argument("-r","--resetids", help = "if option is used, the figures are not normalised", action= "store_false")
## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
targetface = args.target
azim = args.azimuthal
elev = args.elevation
norm = args.nonNormalize
################################
targetsurfacearea= args.surfaceArea
targetid = args.target
othertargetid = args.othersidetissue
################################
areastep = args.areastep
maxeta = args.maxeta
fastkappaOption = args.fastkappa
large  = args.Large
stepsize = 20
maxarea = args.maxarea

jobid = args.jobid
resetids = args.resetids

# For surpressing err
class NullDevice():
	def write(self, s):
		pass
if targetface == None:
	if numOfLayer == 8:
		targetface = 135
	elif numOfLayer == 10:
		targetface = 214

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
directoryName = "plots"
saveDirectory = cwd+"/"+directoryName
if not os.path.exists(saveDirectory):
	os.makedirs(saveDirectory)
# only taking directories that contain data
listdir = [d for d in listdir if d[0] == 'a']
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['fk']) for folder in listdir]
else:
	etalist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
#################################################################################
#     Making the Plots
#################################################################################
import matplotlib.colors as colors
import matplotlib.cm as cmx
etathreshold = 0
jet = cm = plt.get_cmap('viridis') 
##################################################
if maxeta == 0.:
	maxvalue = max(etalist)
else:
	maxvalue = maxeta
minvalue = min(etalist)
cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#fig = plt.figure(frameon=False,figsize=(20,16))
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(111)
##########################
ax1.set_ylabel(r"Mean Roundness, $\langle R \rangle_c$")
ax1.set_xlabel(r"Surface Area, $A_T$")
#################################################################################
counter = 0
totalfolders = len(listdir)
#print listdir
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	else:
		etacurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
	etacolor = scalarMap.to_rgba(etacurrent)
	########################################################
	if (maxeta != 0) and (etacurrent > maxeta):
		continue
	if  (etacurrent %2 != 0):#only even eta analysed for now
		continue
	########################################################
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
	saveName = saveDirectory+'/plot_eta%d_area%d'%(etacurrent,targetsurfacearea)
	plotFeedbackCorrection(targetid, targetsurfacearea,endStep = endStep, 
	startStep=startStep,stepsize= stepsize,maxarea = None, areastep = 10,resetids = resetids,
	saveName = saveName)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	plt.close('all')
	counter+= 1

################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"

