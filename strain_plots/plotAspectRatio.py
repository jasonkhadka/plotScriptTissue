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
def calculateMeanAspectRatio(facelist):
	totalRoundness = 0.
	for face in facelist:
		#########################
		totalAspectRatio += sf.getAspectRatio(face)
		#########################
	return totalAspectRatio/len(facelist)
####################################################################################################################
# Calculating and plotting mean stress and growth
####################################################################################################################
def plotAspectRatio(targetid,othertargetid, targetsurfacearea,
	startStep=1,stepsize= 1,largerCondition =True ,maxarea = None, areastep = 10,
	startarea = None,
	endarea = 850):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	########################################################################
	# faceidarray for Primordia
	if not os.path.isfile("qdObject_step=001.obj"):
		return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
	cell = sf.loadCellFromFile(1,resetids=True)
	initialTissueSurfaceArea = sf.getSurfaceArea(cell)
	#######################################################################
	# Starting the Calculation
	#######################################################################
	#######################################################################
	laststep = 1
	plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
	#######################################################################
	primordialroundnessArray = []
	primordialBoundaryroundnessArray = []
	othertissueroundnessArray = []
	tissueSurfaceAreaArray = []
	if not startarea:#no startarea given
		startarea = int(initialTissueSurfaceArea)
	#######################################################################
	listsurfacearea = np.linspace(startarea,endarea,10)
	###################################################
	for steparea in listsurfacearea:
		step,tissueSurfaceArea = getTimeStep(steparea, endStep, laststep, stepsize = 10)
		########################################################################
		if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
			break
		cell = sf.loadCellFromFile(step,resetids=True)
		################################################
		primordialfaces = sf.getPrimordiaFaces(cell,targetid, large = False)
		primordialBoundaryfaces = sf.getPrimordiaBoundaryFaceList(cell,targetid,large= True)
		othertissuefacelist1 = sf.getPrimordiaFaces(cell,othertargetid, large = False)
		othertissuefacelist2 = sf.getPrimordiaBoundaryFaceList(cell,othertargetid, large = True)
		othertissuefacelist = othertissuefacelist1 + othertissuefacelist2
		################################################
		########################################################################
		# calculate the roundness
		########################################################################
		primordiaRoundness = calculateMeanAspectRatio(primordialfaces)
		primordialBoundaryRoundness = calculateMeanAspectRatio(primordialBoundaryfaces)
		othertissueRoundness = calculateMeanAspectRatio(othertissuefacelist)
		######################################################
		tissueSurfaceAreaArray.append(tissueSurfaceArea)
		primordialroundnessArray.append(primordiaRoundness)
		primordialBoundaryroundnessArray.append(primordialBoundaryRoundness)
		othertissueroundnessArray.append(othertissueRoundness)
		########################################################################
		laststep = step
		########################################################################
		print tissueSurfaceArea, step 
	return [tissueSurfaceAreaArray,primordialroundnessArray,
			primordialBoundaryroundnessArray,othertissueroundnessArray]

####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--startarea", help="Start of area",default =None, type = int)
parser.add_argument('-e',"--endarea", help="area end",default = 850, type = int)
parser.add_argument("-m","--maxeta", help = "if this is given, then eta is only cacluated till this value", type = float, default = 0.0)
parser.add_argument("-x","--maxarea", help = "if this is given, then plot is only made till this area value value", type = float, default = None)
parser.add_argument('-c',"--surfaceArea", help="The surface area for which the stress vs feedback plot would need to be plotStrainDifferenceSurface",
					default =800, type = int)
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
parser.add_argument("-o","--othersidetissue",help = "target cell for other side of the tissue (non-primordia)", type = int)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)
parser.add_argument('-d',"--areastep", help="area step for calculating the growth in cell area", type = int,
						default = 10)
parser.add_argument('-j',"--jobid", help="jobid", type = int,default = None)
## Getting the arguments 
args = parser.parse_args()
#location = args.location
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
targetarea = args.surfaceArea
targetid = args.target
othertargetid = args.othersidetissue
################################
areastep = args.areastep
maxeta = args.maxeta
fastkappaOption = args.fastkappa
large  = args.Large
stepsize = 10
maxarea = args.maxarea
startarea = args.startarea
endarea =args.endarea
jobid = args.jobid

endStep = 2000
startStep = 1


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
jet = cm = plt.get_cmap('plasma') 
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
ax1.set_ylabel(r"Mean Aspect Ratio, $\langle R \rangle_c$")
ax1.set_xlabel(r"Surface Area, $A_T$")
#################################################################################
counter = 0
totalfolders = len(listdir)
plotData = {}
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
	plotData[etacurrent] = plotAspectRatio(targetid=targetid,othertargetid=othertargetid, targetsurfacearea=targetarea,
	startStep=startStep,areastep = areastep,startarea = startarea,
				endarea = endarea)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#plotting
###############################################################################
plotargs = {"linewidth":1}
for key,data in plotData.iteritems():
	color = scalarMap.to_rgba(key)
	##################################
	#mean stress
	##################################
	#data = [tissueSurfaceAreaArray,primordialroundnessArray,
	#		primordialBoundaryroundnessArray,othertissueroundnessArray]
	ax1.plot(data[0], data[1],"*-" ,markersize = 10,label=r"primordium",c=color,**plotargs)
	#ortho Stress
	ax1.plot(data[0], data[2],'<--',markersize = 10, label=r"boundary",c=color,**plotargs)
	ax1.plot(data[0], data[3],'s:', markersize = 10,label=r"meristem",c=color,**plotargs)
	##################################
############################################################
# Legend of the plot
############################################################
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker = '*', color='k', label=r"primordium",**plotargs),
				   Line2D([0], [0],  marker = '<',linestyle = "--", color='k', label=r"boundary",**plotargs),
				   Line2D([0], [0], marker = 's',linestyle = ":", color='k', label=r"meristem",**plotargs),
				   ]
ax1.legend(handles = legend_elements)
ax1.set_ylim(0.8,1.0)
###############################################################################
#color bar fig
###############################################################################
scalarMap._A = []
fig.subplots_adjust(bottom=0.235)
cbar_ax = fig.add_axes([0.15, 0.07, 0.7, 0.03])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
plt.tight_layout()
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	clrbar.set_label(r"Fast Growth Rate")
else:
	clrbar.set_label(r"Mechanical Feedback, $\eta$")

"""##########################################################################################
# Making the legend
##########################################################################################
from collections import OrderedDict
handles, labels = ax1.get_legend_handles_labels()
by_label = OrderedDict(zip(labels,handles))
ax1.legend(by_label.values(),by_label.keys(),prop={'size': 14})


handles, labels = ax2.get_legend_handles_labels()
by_label = OrderedDict(zip(labels,handles))
ax2.legend(by_label.values(),by_label.keys(),prop={'size': 14})
"""
"""###############################################################################
#color bar fig1
scalarMap._A = []
fig1.subplots_adjust(bottom=0.2)
cbar_ax = fig1.add_axes([0.15, 0.05, 0.7, 0.02])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
clrbar.set_label("$\eta$")
###############################################################################
#color bar fig2
scalarMap._A = []
fig2.subplots_adjust(bottom=0.2)
cbar_ax = fig2.add_axes([0.15, 0.05, 0.7, 0.02])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
clrbar.set_label("$\eta$")"""
#plt.tight_layout()
#plt.tight_layout( rect=[0, 0, 1, 1])
#fig.tight_layout(rect=[0.1,0.1,1.,0.9])
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	fig.savefig(saveDirectory+r"/plot_cell_aspectratio.png",transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plot_cell_aspectratio.png",transparent = True, bbox_inches="tight")



#fig1.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_areaPrimodia_%d.png"%endStep,transparent = True)
#fig2.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_surfaceratio_%d.png"%endStep,transparent = True)
#fig3.savefig(saveDirectory+r"/plot_eta_vs_sphericity_%d.png"%endStep,transparent = True)
plt.close('all')
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	if jobid:
		np.save('job=%d_aspectratio_fk_time=%d_targetface=%d.npy'%(jobid,endStep,targetid),plotData)
	else:
		np.save('aspectratio_fk_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)
else:
	if jobid:
		np.save('job=%d_aspectratio_eta_time=%d_targetface=%d.npy'%(jobid,endStep,targetid),plotData)
	else:
		np.save('aspectratio_eta_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)



################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"

