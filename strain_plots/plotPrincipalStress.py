import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sop
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
import string
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['legend.fontsize'] = 22
plt.rcParams['axes.titlesize'] = 22
####################################################################################################################
# Add subplot annotation
####################################################################################################################
def addAnnotation(subplot,n = 1):
	subplot.text(-0.1,1.1,string.ascii_lowercase[n],transform = ax.transAxes,size = 20, weight = 'bold')
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
                        cell = sf.loadCellFromFile(calstep+1)
                        tissueSurfaceArea = sf.getSurfaceArea(cell)
                        return calstep+1,tissueSurfaceArea
        ################################################
        gc.collect()
    return endStep,tissueSurfaceArea
###############################################################################################################
# Calculate primordia height
###############################################################################################################
def getPrimordiaHeight(cell, targetid):
    ###################################################################
    def addMeanVertex(vertex,meanx,meany,meanz):
        meanx += vertex.getXcoordinate()
        meany += vertex.getYcoordinate()
        meanz += vertex.getZcoordinate()
        return meanx,meany,meanz
    ########################################################################
    # Getting the primordial boundary
    ########################################################################
    facetarget = sf.getFace(cell, targetid)
    ##########################################
    # Vertex on primordial boundary
    ##########################################
    vertexList = sf.getPrimordiaBoundaryVertexList(cell, targetid)
    vertexNum = len(vertexList)
    ####################################################
    # Calculation of primordial height starts here
    # This is for smaller primordia
    ####################################################
    meanx = 0.
    meany = 0.
    meanz = 0.
    for vert in vertexList:#while edge.Dest().getID() != targetedge.Dest().getID():
        meanx,meany,meanz = addMeanVertex(vert,meanx,meany,meanz)
    ######################################
    targetx = facetarget.getXCentralised()
    targety = facetarget.getYCentralised()
    targetz = facetarget.getZCentralised()
    meanx /= vertexNum
    meany /= vertexNum
    meanz /= vertexNum
    height = np.sqrt((meanx-targetx)**2+(meany-targety)**2+(meanz-targetz)**2)
    ######################################
    return height
####################################################################################################################
# fit function
####################################################################################################################
def fitLinFunc(t,m,c):
	return m*t + c
####################################################################################################################
# calculating growth ratio
####################################################################################################################
def getGrowthRatio(numOfLayer, targetid ,endStep ,startStep = 1,stepsize = 5):
	if not os.path.isfile("qdObject_step=001.obj"):
			return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
	cell = sf.loadCellFromFile(1)
	########################################################################
	meanprimordiaArray = []
	meanrestArray = []
	timeArray = []
	########################################################################
	fitlen = 50
	finalstep = startStep + stepsize*fitlen
	for step in range(startStep, finalstep, stepsize):
		########################################################################
		if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
			break
		cell = sf.loadCellFromFile(step)
		################################################
		primordiafacelist  = sf.getPrimordiaFaces(cell, targetid, large = False)
		primordiaarea = 0.
		for face in primordiafacelist:
			primordiaarea += face.getAreaOfFace()
		################################################
		tissueSurfaceArea = sf.getSurfaceArea(cell)
		################################################
		primordialface = sf.getFace(cell, targetid)
		restoftissuearea =  tissueSurfaceArea - primordiaarea
		################################################
		numOfPrimordialcell = len(primordiafacelist)
		numOfrestofcell = cell.countFaces() -1 -numOfPrimordialcell
		################################################
		meanprimordiafacearea = primordiaarea/numOfPrimordialcell
		meanrestoftissuefacearea = restoftissuearea/(numOfrestofcell)
		################################################
		meanprimordiaArray.append(meanprimordiafacearea)
		meanrestArray.append(meanrestoftissuefacearea)
		timeArray.append(step-1)
		################################################
	logfastarea = np.log(meanprimordiaArray)
	logslowarea = np.log(meanrestArray)
	################################################
	fastareafit, m = sop.curve_fit(fitLinFunc,timeArray[:fitlen],logfastarea[:fitlen],bounds=([-np.inf,logfastarea[0]-0.000001],[+np.inf,logfastarea[0]]))
	slowareafit, m = sop.curve_fit(fitLinFunc,timeArray[:fitlen],logslowarea[:fitlen],bounds=([-np.inf,logslowarea[0]-0.000001],[+np.inf,logslowarea[0]]))
	################################################
	return fastareafit[0]/slowareafit[0]
####################################################################################################################
# Calculating and plotting mean stress and growth
####################################################################################################################
def plotPrincipalStress(numOfLayer, targetid,endStep,eta, 
	meanstressplot, 
	color,startStep=0,stepsize= 1,largerCondition =True ,maxarea = None, areastep = 20,
	startarea = 680,
	endarea = 850):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	########################################################################
	# faceidarray for Primordia
	if not os.path.isfile("qdObject_step=001.obj"):
		return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
	cell = sf.loadCellFromFile(1)
	initialTissueSurfaceArea = sf.getSurfaceArea(cell)
	#######################################################################
	# Starting the Calculation
	#######################################################################
	laststep = 1
	plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
	#######################################################################
	tissueSurfaceAreaArray = []
	meanstressEigenvalue1Array = []
	meanstressEigenvalue2Array = []
	heightArray = []
	primordialAreaArray = []
	for steparea in range(startarea, endarea, int(areastep)):
		step,tissueSurfaceArea = getTimeStep(steparea, endStep, laststep, stepsize = 10)
		########################################################################
		if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
			break
		cell = sf.loadCellFromFile(step)
		################################################
		cell.calculateStressStrain()
		################################################
		primordialface = sf.getFace(cell, targetid)
		################################################
		cell.setRadialOrthoradialVector(primordialface)
		cell.setRadialOrthoradialStress()
		################################################
		########################################################################
		#  Starting the Calculation of mean growth and mean stress on boundary
		########################################################################
		# mean stress
		########################################################################
		faceList = sf.getPrimordiaBoundaryFaceList(cell,targetid,large= large)
		height = getPrimordiaHeight(cell,targetid)
		###############################################
		primordiafacelist  = sf.getPrimordiaFaces(cell, targetid, large = False)
		primordiaarea = 0.
		for face in primordiafacelist:
			primordiaarea += face.getAreaOfFace()
		######################################################
		#radialDict, orthoradialDict = getRadialOrthoradialDict(cell,targetid,large = large)
		######################################################
		stressEigenvalue1Array = []
		stressEigenvalue2Array = []
		for face in faceList:
			stresseigenvalue1 = face.getStressEigenValue1()
			stresseigenvalue2 = face.getStressEigenValue2()
			#######################################################
			stressEigenvalue1Array.append(stresseigenvalue1)
			stressEigenvalue2Array.append(stresseigenvalue2)
		######################################################
		tissueSurfaceAreaArray.append(tissueSurfaceArea)
		heightArray.append(height)
		######################################################
		meanstressEigenvalue1Array.append(np.mean(stressEigenvalue1Array))
		meanstressEigenvalue2Array.append(np.mean(stressEigenvalue2Array))
		primordialAreaArray.append(primordiaarea)
		########################################################################
		laststep = step
		########################################################################
	return [tissueSurfaceAreaArray,meanstressEigenvalue1Array, meanstressEigenvalue2Array, 
			np.add(meanstressEigenvalue1Array, meanstressEigenvalue2Array), 
			np.subtract(meanstressEigenvalue2Array,meanstressEigenvalue1Array),
			heightArray,primordialAreaArray]
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-m","--maxeta", help = "if this is given, then eta is only cacluated till this value", type = float, default = 0.0)
parser.add_argument("-x","--maxarea", help = "if this is given, then plot is only made till this area value value", type = float, default = None)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
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
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
											  default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = None, type = int)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)
parser.add_argument('-d',"--areastep", help="area step for calculating the growth in cell area", type = int,default = 20)

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
areastep = args.areastep
maxeta = args.maxeta
fastkappaOption = args.fastkappa
large  = args.Large
stepsize = 10
maxarea = args.maxarea
startarea = 680
endarea = 850
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
##################################################
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
fig = plt.figure(figsize=(10,15))
fig2 = plt.figure(figsize=(11,5))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)

ax8 = fig2.add_subplot(121)
ax7 = fig2.add_subplot(122)

##########################
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Stress
##################################
ax1.set_xlabel(r"Surface Area, $A_T$")
ax1.set_ylabel(r"Principal Stresses, $\sigma_1+\sigma_2$")


ax2.set_xlabel(r"Surface Area, $A_T$")
ax2.set_ylabel(r"Principal Stresses, $\sigma_1$")

ax3.set_xlabel(r"Surface Area, $A_T$")
ax3.set_ylabel(r"Principal Stresses, $\sigma_2$")

ax4.set_xlabel(r"Surface Area, $A_T$")
ax4.set_ylabel(r"Principal Stresses, $\sigma_2-\sigma_1$")

ax5.set_xlabel(r"Primordial height, $h$")
ax5.set_ylabel(r"Principal Stresses, $\sigma_1+\sigma_2$")

ax6.set_xlabel(r"Surface Area, $A_T$")
ax6.set_ylabel(r"Primordial height, $h$")


ax8.set_xlabel(r"Surface Area, $A_T$")
ax8.set_ylabel(r"Primordial Area, $A_P$")

ax7.set_xlabel(r"Surface Area, $A_T$")
ax7.set_ylabel(r"Primordial height, $h$")

########################################################
growthRatio = {}
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	for folder in listdir:
		# Converting folder name to dictionary
		fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
		########################################################
		if (maxeta != 0) and (fkcurrent > maxeta):
			continue
		########################################################
		os.chdir(folder)
		########################################################
		growthRatio[fkcurrent] = getGrowthRatio(numOfLayer = numOfLayer, targetid = targetid,
			endStep = endStep,startStep = startStep)
		#print sys.getsizeof(plotData)
		os.chdir("..")
		gc.collect()
	########################################################
	maxvalue= max(growthRatio.values())
	minvalue =int(min(growthRatio.values()))
	########################################################
	##################################################
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	##################################################
##################################################
counter = 0
totalfolders = len(listdir)
plotData = {}
#print listdir
for folder in listdir:
	# Converting folder name to dictionary
	#print folder
	if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	   etacurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
	   etacurrent= growthRatio[etacurrent]
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
	plotData[etacurrent] = plotPrincipalStress(numOfLayer, targetid,endStep,etacurrent, 
	meanstressplot=ax1, 
	color = etacolor,startStep=0,stepsize= 1,largerCondition =True ,maxarea = None, areastep = 20,
	startarea = startarea,
	endarea = endarea)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#plotting
###############################################################################
plotargs = {"linewidth":3}
for key,data in plotData.iteritems():
	color = scalarMap.to_rgba(key)
	##################################
	#mean stress
	##################################
	#rad Stress
	ax1.plot(data[0], data[3],"-." ,label=r"$\sigma$",c=color,**plotargs)
	ax2.plot(data[0], data[1],"-." ,label=r"$\sigma$",c=color,**plotargs)
	ax3.plot(data[0], data[2],"-." ,label=r"$\sigma$",c=color,**plotargs)
	ax4.plot(data[0], data[4],"-." ,label=r"$\sigma$",c=color,**plotargs)
	ax5.plot(data[5], data[3],"-." ,label=r"$\sigma$",c=color,**plotargs)
	ax6.plot(data[0], data[5],"-." ,label=r"$\sigma$",c=color,**plotargs)
	##################################
	ax7.plot(data[0], data[5] ,label=r"$\sigma$",c=color,**plotargs)
	ax8.plot(data[0], data[6],"--" ,label=r"$\sigma$",c=color,**plotargs)
	#ortho Stress
	#ax1.plot(data[0], data[2], label=r"$\sigma_{2}$",c=color,**plotargs)
############################################################
# Legend of the plot
############################################################
#from matplotlib.lines import Line2D
#legend_elements = [Line2D([0], [0], linestyle = "-.", color='k', label=r"$\sigma_{1}$",**plotargs),
#				   Line2D([0], [0],  color='k', label=r"$\sigma_{2}$",**plotargs)]
#ax1.legend(handles = legend_elements)

###############################################################################
#color bar fig
###############################################################################
#plt.tight_layout()
scalarMap._A = []
#clrbar = plt.colorbar(scalarMap, ax=ax1,shrink = 0.9,aspect = 8,ticks=np.linspace(minvalue, maxvalue, 4).astype('int'))#,orientation='horizontal',cax = cbar_ax)
#clrbar1 = plt.colorbar(scalarMap, ax=ax7,shrink = 0.9,aspect = 8,ticks=np.linspace(minvalue, maxvalue, 4).astype('int'))#,orientation='horizontal',cax = cbar_ax)

################################################################################
axpos1 = ax1.get_position()
axpos7 = ax7.get_position()
################################################################################
#print axpos
fig2.subplots_adjust(right=0.9)
fig2.tight_layout(rect=[0.,0.,.9,.9])
#fig2.tight_layout(rect=[0.,0.,.9,.9])
cbar_ax7 = fig2.add_axes([0.91, 0.15, 0.025, 0.7])
clrbar1 = plt.colorbar(scalarMap,cax = cbar_ax7,ticks=np.linspace(minvalue, maxvalue, 3).astype('int'))
clrbarpos1 = [axpos1.x0+axpos1.width,axpos1.y0,0.04,axpos1.height]
#clrbarpos7 = [axpos7.x0+axpos7.width,axpos7.y0,0.04,axpos7.height]
################################################################################
cbar_ax1 = fig.add_axes(clrbarpos1)
#cbar_ax7 = fig2.add_axes(clrbarpos7)
################################################################################
clrbar = plt.colorbar(scalarMap,cax = cbar_ax1,ticks=np.linspace(minvalue, maxvalue, 3).astype('int'))
clrbar1 = plt.colorbar(scalarMap,cax = cbar_ax7,ticks=np.linspace(minvalue, maxvalue, 3).astype('int'))
################################################################################
ax7.set_xticks(np.linspace(startarea,endarea,3))
ax8.set_xticks(np.linspace(startarea,endarea,3))
################################################################################

if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	clrbar.set_label(r"Growth ratio, $r_g$")
	clrbar1.set_label(r"Growth ratio, $r_g$")

else:
	clrbar.set_label(r"Mechanical Feedback, $\eta$")


fig.savefig(saveDirectory+r"/plot_principalStress_targetface=%d.png"%(targetid),transparent = True, bbox_inches="tight")
fig2.savefig(saveDirectory+r"/plot_growthRatio_primordialGrowth_targetid=%d.png"%(targetid),transparent = True, bbox_inches="tight")
fig2.savefig(saveDirectory+r"/plot_growthRatio_primordialGrowth_targetid=%d.eps"%(targetid),transparent = True, bbox_inches="tight")



#fig1.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_areaPrimodia_%d.png"%endStep,transparent = True)
#fig2.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_surfaceratio_%d.png"%endStep,transparent = True)
#fig3.savefig(saveDirectory+r"/plot_eta_vs_sphericity_%d.png"%endStep,transparent = True)
plt.close('all')
### Saving Data Dictionary ###
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	np.save('principalstress_fk_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)
else:
	np.save('principalstress_eta_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)


################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"

