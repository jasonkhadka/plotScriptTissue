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
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['legend.fontsize'] =14
plt.rcParams['axes.titlesize'] = 14.


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
def getRadialOrthoradialStress(face, radialDict=None,orthoradialDict=None, vectors = False):
    #eigenvec1 = face.getStressEigenVector1()
    #eigenvec2 = face.getStressEigenVector2()
    eigenvalue1 = face.getStressEigenValue1()
    eigenvalue2 = face.getStressEigenValue2()
    """
    vec1unit = np.array([qd.doublearray_getitem(eigenvec1,0),
                    qd.doublearray_getitem(eigenvec1,1),
                    qd.doublearray_getitem(eigenvec1,2)])
    vec2unit = np.array([qd.doublearray_getitem(eigenvec2,0),
                    qd.doublearray_getitem(eigenvec2,1),
                    qd.doublearray_getitem(eigenvec2,2)])
    vec1 =eigenvalue1*vec1unit
    vec2 = eigenvalue2*vec2unit
    radialvec = np.copy(radialDict[face.getID()])
    orthoradialvec = np.copy(orthoradialDict[face.getID()])
    #print radialvec,orthoradialvec
    ############################################
    ########################################################################################
    radsign1 = 1.#(np.dot(radialvec, vec1unit)<0.)*(-1)+(np.dot(radialvec, vec1unit)>0.)*(1)
    radsign2 = 1.#(np.dot(radialvec, vec2unit)<0.)*(-1)+(np.dot(radialvec, vec2unit)>0.)*(1)
    orthosign1 = 1.#(np.dot(orthoradialvec, vec1unit)<0.)*(-1)+(np.dot(orthoradialvec, vec1unit)>0.)*(1)
    orthosign2 = 1.#(np.dot(orthoradialvec, vec2unit)<0.)*(-1)+(np.dot(orthoradialvec, vec2unit)>0.)*(1)
    #print radsign1, radsign2, orthosign1, orthosign2
    ########################################################################################
    radialComp = eigenvalue1*np.dot(radialvec,vec1unit)*radsign1+eigenvalue2*np.dot(radialvec,vec2unit)*radsign2
    orthoradialComp = #eigenvalue1*np.dot(orthoradialvec,vec1unit)*orthosign1+eigenvalue2*np.dot(orthoradialvec,vec2unit)*orthosign2
    #radialComp = np.dot(radialvec,vec1)+np.dot(radialvec,vec2)
    #orthoradialComp = np.dot(orthoradialvec,vec1)+np.dot(orthoradialvec,vec2)
    ############################################
    if vectors:
        radialvec = radialComp*radialvec
        orthoradialvec = orthoradialComp*orthoradialvec
        return radialvec, orthoradialvec
    """
    ############################################
    radialComp = face.getRadialStress()
    orthoradialComp = face.getOrthoradialStress()
    ############################################
    return radialComp, orthoradialComp, eigenvalue1, eigenvalue2#radialvec,orthoradialvec
###################################################################
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
# for a face, projecting its Growth eigendecomposed vectors onto radial-orthoradial direction
###############################################################################################################
def getRadialOrthoradialGrowth(face, radialDict=None,orthoradialDict=None, vectors = False):
    #eigenvec1 = face.getRotGrowthEigenVector1()
    #eigenvec2 = face.getRotGrowthEigenVector2()
    eigenvalue1 = face.getRotGrowthEigenValue1()
    eigenvalue2 = face.getRotGrowthEigenValue2()
    """
    vec1unit = np.array([qd.doublearray_getitem(eigenvec1,0),
                    qd.doublearray_getitem(eigenvec1,1),
                    qd.doublearray_getitem(eigenvec1,2)])
    vec2unit = np.array([qd.doublearray_getitem(eigenvec2,0),
                    qd.doublearray_getitem(eigenvec2,1),
                    qd.doublearray_getitem(eigenvec2,2)])
    vec1 =eigenvalue1*vec1unit
    vec2 = eigenvalue2*vec2unit
    #print vec1unit
    radialvec = np.copy(radialDict[face.getID()])
    orthoradialvec = np.copy(orthoradialDict[face.getID()])
    #print radialvec,orthoradialvec
    ########################################################################################
    # sign change if needed : 
    #       if EigenVector and RadialVec are opposite facing
    ########################################################################################
    #print "rad.vec1 :",np.dot(radialvec, vec1unit), 'rad.vec2',np.dot(radialvec, vec2unit), np.dot(vec1unit,vec2unit)
    ########################################################################################
    radsign1 = 1.#(np.dot(radialvec, vec1unit)<0.)*(-1)+(np.dot(radialvec, vec1unit)>0.)*(1)
    radsign2 = 1.#(np.dot(radialvec, vec2unit)<0.)*(-1)+(np.dot(radialvec, vec2unit)>0.)*(1)
    orthosign1 = 1.#(np.dot(orthoradialvec, vec1unit)<0.)*(-1)+(np.dot(orthoradialvec, vec1unit)>0.)*(1)
    orthosign2 = 1.#(np.dot(orthoradialvec, vec2unit)<0.)*(-1)+(np.dot(orthoradialvec, vec2unit)>0.)*(1)
    #print radsign1, radsign2, orthosign1, orthosign2
    ########################################################################################
    radialComp = eigenvalue1*np.dot(radialvec,vec1unit)*radsign1+eigenvalue2*np.dot(radialvec,vec2unit)*radsign2
    orthoradialComp = eigenvalue1*np.dot(orthoradialvec,vec1unit)*orthosign1+eigenvalue2*np.dot(orthoradialvec,vec2unit)*orthosign2
    #radialComp = np.dot(radialvec,vec1)*radsign1+np.dot(radialvec,vec2)*radsign2
    #orthoradialComp = np.dot(orthoradialvec,vec1)*orthosign1+np.dot(orthoradialvec,vec2)*orthosign2
    ############################################
    
    if vectors:
        radialvec = radialComp*radialvec
        orthoradialvec = orthoradialComp*orthoradialvec
        return radialvec, orthoradialvec
    """
    ############################################
    radialComp = face.getRadialGrowth()
    orthoradialComp = face.getOrthoradialGrowth()
    ############################################
    ############################################
    return radialComp, orthoradialComp, eigenvalue1, eigenvalue2#radialvec,orthoradialvec
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
# Calculating and plotting mean stress and growth
####################################################################################################################
def plotHeightGrowthScatter(numOfLayer, targetid,endStep,eta, 
    stressscatter, growthscatter,stressscatter1, growthscatter1, 
    color,maxeta = 20,startStep=0,stepsize= 1,largerCondition =True ,maxarea = None, areastep = 20):
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
    #####################################
    #Getting initial area of primodial
    #####################################
    if not os.path.isfile("qdObject_step=001.obj"):
        return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
    cell = sf.loadCellFromFile(1)
    #######################################################################
    laststep = 1
    plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
    #######################################################################
    orthoradialStressArray = []
    radialStressArray = []
    radialGrowthArray = []
    orthoradialGrowthArray = []
    tissueSurfaceAreaArray = []
    primordiaAreaArray = []
    heightArray = []
    meanstressEigenvalue1Array = []
    meanstressEigenvalue2Array = []

    meangrowthEigenvalue1Array = []
    meangrowthEigenvalue2Array = []
    for steparea in range(680, 800, int(areastep)):
        step,tissueSurfaceArea = getTimeStep(steparea, endStep, laststep, stepsize = 5)
        ########################################################################
        step2,tissueSurfaceArea2 = getTimeStep(steparea+stepsize, endStep, step, stepsize = 5)
        ########################################################################
        if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
            break
        cell = sf.loadCellFromFile(step)
        cell2 = sf.loadCellFromFile(step2)
        ################################################
        cell.calculateStressStrain()
        ################################################
        primordialface = sf.getFace(cell, targetid)
        ################################################
        cell.setRadialOrthoradialVector(primordialface)
        cell.setRadialOrthoradialStress()
        ################################################
        sf.calculateDilation(cell,cell2)
        ########################################################################
        #  Starting the Calculation of mean growth and mean stress on boundary
        ########################################################################
        # mean stress
        ########################################################################
        faceList = sf.getPrimordiaBoundaryFaceList(cell,targetid,large= large)
        primordiafacelist  = sf.getPrimordiaFaces(cell, targetid, large = False)
        primordiaarea = 0.
        for face in primordiafacelist:
            primordiaarea += face.getAreaOfFace()
        ######################################################
        #radialDict, orthoradialDict = getRadialOrthoradialDict(cell,targetid,large = large)
        ######################################################
        radialStress = []
        orthoradialStress = []
        radialGrowth = []
        orthoradialGrowth = []
        stressEigenvalue1Array = []
        stressEigenvalue2Array = []
        growthEigenvalue1Array = []
        growthEigenvalue2Array = []
        dTissueSurfaceArea = tissueSurfaceArea2-tissueSurfaceArea
        height = getPrimordiaHeight(cell,targetid)
        height2 = getPrimordiaHeight(cell2,targetid)
        dhdA = (height2-height)/dTissueSurfaceArea
        for face in faceList:
            radstress, orthstress,stresseigenvalue1, stresseigenvalue2  = getRadialOrthoradialStress(face)
            radGrowth, orthGrowth, growtheigenvalue1, growtheigenvalue2 = getRadialOrthoradialGrowth(face)
            #######################################################
            # plotting
            #######################################################
            sigma2 = max(stresseigenvalue1,stresseigenvalue2)
            sigma1  = min(stresseigenvalue1,stresseigenvalue2)
            g2 = max(growtheigenvalue1,growtheigenvalue2)
            g1 = min(growtheigenvalue1,growtheigenvalue2)
            #######################################################
            stressscatter.scatter(sigma2-sigma1, dhdA, c = color,marker = 'o',alpha = 0.7,zorder= maxeta-eta)
            growthscatter.scatter(g2-g1, dhdA, c = color,marker = 'o',alpha = 0.7,zorder= maxeta-eta)
            #######################################################
            stressscatter1.scatter(stresseigenvalue2-stresseigenvalue1, dhdA, c = color,marker = 'o',alpha = 0.7,zorder= maxeta-eta)
            growthscatter1.scatter(growtheigenvalue2-growtheigenvalue1, dhdA, c = color,marker = 'o',alpha = 0.7,zorder= maxeta-eta)
            #######################################################
        ########################################################################
        laststep = step
        ########################################################################
        print tissueSurfaceArea, tissueSurfaceArea2,dTissueSurfaceArea, step , step2
    return 
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
fig = plt.figure(figsize=(10,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
##########################
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Stress Scatter plot
##################################
ax1.set_xlabel(r"$\sigma_2-\sigma_1$")
ax1.set_ylabel(r"$\frac{\Delta h}{\Delta A_T}$")

ax2.set_xlabel(r"$g_2-g_1$")
ax2.set_ylabel(r"$\frac{\Delta h}{\Delta A_T}$")

ax3.set_xlabel(r"$\sigma_2-\sigma_1$")
ax3.set_ylabel(r"$\frac{\Delta h}{\Delta A_T}$")

ax4.set_xlabel(r"$g_2-g_1$")
ax4.set_ylabel(r"$\frac{\Delta h}{\Delta A_T}$")


########################################################
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
	plotHeightGrowthScatter(numOfLayer = numOfLayer, targetid = targetid,endStep = endStep,eta = etacurrent,
				stressscatter = ax1,growthscatter = ax2,stressscatter1 =ax3,growthscatter1 = ax4,
                startStep = startStep,  maxeta = maxeta,
				color = etacolor,stepsize = stepsize,
                largerCondition = large,maxarea = maxarea, areastep = areastep)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#color bar fig
###############################################################################
plt.tight_layout()
scalarMap._A = []
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.15, 0.07, 0.7, 0.03])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	clrbar.set_label(r"Fast Growth Rate")
else:
	clrbar.set_label(r"$\eta$")

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
	fig.savefig(saveDirectory+r"/plot_anistropy_heightgrowth_scatterplot_time=%d_targetface=%d.png"%(endStep,targetid),transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plot_anistropy_heightgrowth_scatterplot_time=%d_targetface=%d.png"%(endStep,targetid),transparent = True, bbox_inches="tight")



#fig1.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_areaPrimodia_%d.png"%endStep,transparent = True)
#fig2.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_surfaceratio_%d.png"%endStep,transparent = True)
#fig3.savefig(saveDirectory+r"/plot_eta_vs_sphericity_%d.png"%endStep,transparent = True)
plt.close('all')
### Saving Data Dictionary ###

################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"
