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
import primordiaplotfunctions as ppf
import simulation_functions as sf
import argparse #argument parser, handles the arguments passed by command line
import gc
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 18.
plt.rcParams['ytick.labelsize'] = 18.
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.titlesize'] = 18

def getEdge(cell, vert1id, vert2id):
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	while vertex != None:
		if vertex.getID() == vert1id:
			break
		vertex= vertices.next()
	###################################
	edges = qd.VertexEdgeIterator(vertex)
	edge = edges.next()
	while edge != None:
		if edge.Dest().getID() == vert2id:
			break
		edge = edges.next()
	return edge
###################################################################
def addMean(edge,meanx,meany,meanz):
    meanx += edge.Dest().getXcoordinate()
    meany += edge.Dest().getYcoordinate()
    meanz += edge.Dest().getZcoordinate()
    return meanx,meany,meanz
###################################################################
def addMeanVertex(vertex,meanx,meany,meanz):
    meanx += vertex.getXcoordinate()
    meany += vertex.getYcoordinate()
    meanz += vertex.getZcoordinate()
    return meanx,meany,meanz

####################################################################################################
# Get array of all vertices on boundary of primordia
####################################################################################################
def getPrimordiaBoundaryVertexList(cell, targetface, large=False):
    ########################################################
    """
    Get a list of all vertices arround a primordia
    targetid : center of primordia
    large : if true, large primordia is calculated for (2 layers)
            if false, small primordia (1 layer)
    """
    face = sf.getFace(cell, targetface)
    edge = face.getEdge()
    vertexList = []
    ########################################################
    # tracing the primordial boundary
    ########################################################
    if large:
        #########################################################
        # Larger Proimordia
        #########################################################
        for _ in range(2):
            edge = edge.Rprev()
        #####################
        for _ in range(3):
            edge = edge.Lnext()
        ####################################################
        # listing Primordia vertex starts here
        ####################################################
        for _ in range(6):
            edge = edge.Rprev()
            vertexList.append(edge.Dest())
            #####################
            for _ in range(2):
                edge = edge.Lnext()
                vertexList.append(edge.Dest())
            #####################
            edge = edge.Rprev()
            vertexList.append(edge.Dest())
            #####################
            edge = edge.Lnext()
            vertexList.append(edge.Dest())
    else:
        #########################################################
        # Smaller Proimordia
        #########################################################
        for _ in range(1):
            edge = edge.Rprev()
        ###############################################################
        # listing Primordia vertex starts here
        ###############################################################
        for _ in range(3):
            edge = edge.Lnext()
            vertexList.append(edge.Dest())
        #####################
        for _ in range(5):
            edge = edge.Rprev()
            vertexList.append(edge.Dest())
            ####################
            edge = edge.Lnext()
            vertexList.append(edge.Dest())
            ####################
            edge = edge.Lnext()
            vertexList.append(edge.Dest())
            ####################
    return vertexList
###################################################################
# Calculating the Average face area
###################################################################
def getFaceAreaData(cell,faceidarray):
    cell.setParameters()
    facearea = 0.
    tfmDet = 0.
    ######################################
    slowfacearea = 0.
    fastfacearea = 0.
    ######################################
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    numfaces = 0.
    fastnum = 0.
    while face != None:
        if face.getID() == 1 : 
            face = faces.next()
            continue
        tfmDet += face.getTargetFormMatrixDeterminant()
        #print face.getID(), area0
        facearea+= face.getAreaOfFace()
        if np.isnan(face.getAreaOfFace()):
            print "Step :", i, "face ID : ", face.getID()
        ###################################################
        # Saving in slow and fast array
        ###################################################
        if face.getID() in faceidarray:
            fastfacearea += face.getAreaOfFace()
            fastnum += 1
        else:
            slowfacearea += face.getAreaOfFace()
        numfaces += 1.
        face  = faces.next()
    return tfmDet,slowfacearea,fastfacearea,fastfacearea

###################################################################
def plotMinGaussianCurvaturePrimodiaHeight(numOfLayer, targetid,endStep,eta, 
    mincurvatureplot, heightplot,ax3,ax4,ax5,ax6,ax7,
    color,startStep=0,stepsize= 1,largerCondition = False,maxarea = None,resetids = True):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    ########################################################################
    # faceidarray for Primordia
    if not os.path.isfile("qdObject_step=000.obj"):
        return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
    cell = sf.loadCellFromFile(0,resetids = resetids)
    faceList = sf.getPrimordiaFaces(cell,targetid, large = largerCondition)
    faceidarray = [xface.getID() for xface in faceList]
    primordialFaceNum  = len(faceList)
    slowFaceNum = cell.countFaces()-primordialFaceNum- 1.
    #print largerCondition, faceidarray
    #######################################################################
    # Checking if the files exists if not going to step down
    #######################################################################
    meanGaussianCurvature = []
    primodialheight = []
    primodialAreaArray = []
    tissueSurfaceAreaArray = []
    surfaceAreaRatio = []
    sphericityArray = []
    tissueVolumeArray = []
    timestep = []
    ######################################################
    # Gathering face area
    ######################################################
    m0determinantArray = []
    slowarray = []
    fastarray = []
    #######################################################################
    # Starting the Calculation
    #######################################################################
    #####################################
    #Getting initial area of primodial
    #####################################
    if not os.path.isfile("qdObject_step=001.obj"):
        return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
    cell = sf.loadCellFromFile(1,resetids = resetids)
    #################################
    # face area data
    #################################
    tfmdet, slowfacearea, fastfacearea,areaPrimodia = getFaceAreaData(cell,faceidarray)
    ##################################################################
    areaInitialPrimodia = areaPrimodia
    #######################################################################
    for step in range(startStep,endStep+1,stepsize):
        if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
            break
        cell = sf.loadCellFromFile(step,resetids = resetids)
        #################################
        # face area data
        #################################
        tfmdet, slowfacearea, fastfacearea,areaPrimodia = getFaceAreaData(cell,faceidarray)
        m0determinantArray.append(tfmdet)
        ################################################################################
        # saving the mean area 
        ################################################################################
        slowarray.append(slowfacearea/slowFaceNum)
        fastarray.append(fastfacearea/primordialFaceNum)
        ########################################################################
        #  Starting the Calcuation of Primordial Height & Curvature            #
        ########################################################################
        gaussianCurvature = []
        meanpointx = []
        meanpointy = []
        meanpointz = []
        tissueSurfaceArea = sf.getSurfaceArea(cell)
        tissueVolume = cell.getCartesianVolume()
        sphericity = (np.pi**(1./3.)*(2.**(1./2.)*3*tissueVolume)**(2./3.))/(tissueSurfaceArea)
        ########################################################################
        # Getting the primordial boundary
        ########################################################################
        facetarget = sf.getFace(cell, targetid)
        ##########################################
        # Vertex on primordial boundary
        ##########################################
        vertexList = getPrimordiaBoundaryVertexList(cell, targetface=targetid,large = largerCondition)
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
            gaussianCurvature.append(vert.getGaussianCurvature())
        ######################################
        targetx = facetarget.getXCentralised()
        targety = facetarget.getYCentralised()
        targetz = facetarget.getZCentralised()
        meanx /= vertexNum
        meany /= vertexNum
        meanz /= vertexNum
        height = np.sqrt((meanx-targetx)**2+(meany-targety)**2+(meanz-targetz)**2)
        ##########################################################################
        primodialheight.append(height)
        tissueSurfaceAreaArray.append(tissueSurfaceArea)
        primodialAreaArray.append(areaPrimodia)
        surfaceAreaRatio.append((areaPrimodia/(tissueSurfaceArea)))
        sphericityArray.append(sphericity)
        meanGaussianCurvature.append(np.mean(np.array(gaussianCurvature)))
        tissueVolumeArray.append(tissueVolume)
        timestep.append(step-1.)
    #################################################################################
    #                         Plotting 
    #################################################################################
    # calculating the plotlen
    if maxarea:
        plotlen = ppf.getPlotlenMaxArea(tissueSurfaceAreaArray,maxarea)
        #print timestep, primodialheight, meanGaussianCurvature
        # Min Gaussian curvature
        #print timestep
        mincurvatureplot.plot(tissueSurfaceAreaArray[:plotlen],meanGaussianCurvature[:plotlen],c=color,lw = 1.5)
        ###################################
        # Height of Primodia
        heightplot.plot(tissueSurfaceAreaArray[:plotlen], primodialheight[:plotlen], c=color,lw = 1.5)
        ###################################
        # primordial area vs surface area
        ax3.plot(tissueSurfaceAreaArray[:plotlen], primodialAreaArray[:plotlen], c = color, lw = 1.5)
        ###################################
        # surface area vs time
        ax4.plot(timestep[:plotlen],tissueSurfaceAreaArray[:plotlen], c = color, lw = 1.5)
    #print timestep, primodialheight, meanGaussianCurvature
    # Min Gaussian curvature
    #print timestep
    mincurvatureplot.plot(tissueSurfaceAreaArray,meanGaussianCurvature,c=color,lw = 1.5)
    ###################################
    # Height of Primodia
    heightplot.plot(tissueSurfaceAreaArray, primodialheight, c=color,lw = 1.5)
    ###################################
    # primordial area vs surface area
    ax3.plot(tissueSurfaceAreaArray, primodialAreaArray, c = color, lw = 1.5)
    ###################################
    # surface area vs time
    ax4.plot(timestep,tissueSurfaceAreaArray, c = color, lw = 1.5)
    ########################################################
    #ax3.plot(primodialAreaArray,meanGaussianCurvature,c=color,lw =1.5)
    #ax4.plot(primodialAreaArray,primodialheight,c=color,lw = 1.5)
    ########################################################
    #ax5.plot(surfaceAreaRatio,meanGaussianCurvature,c=color,lw =1.5)
    #ax6.plot(surfaceAreaRatio,primodialheight,c=color,lw = 1.5)
    ########################################################
    #ax7.plot(timestep, sphericityArray,c=color,lw = 1.5)
    return [primodialheight, meanGaussianCurvature,primodialAreaArray,
    tissueSurfaceAreaArray,tissueVolumeArray, 
    sphericityArray,m0determinantArray,
     fastarray, slowarray,timestep]
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
parser.add_argument("-d","--stepsize", help = "stepsize on plot", default = 1, type = int)
parser.add_argument('-j',"--jobid", help="jobid", type = int,default = None)
parser.add_argument("-r","--resetids", help = "if option is used, the figures are not normalised", action= "store_false")

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
stepsize = args.stepsize
maxeta = args.maxeta
fastkappaOption = args.fastkappa
large  = args.Large
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
fig = plt.figure(figsize=(10,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
ax1.set_title("Mean Gaussian Curvature")
ax1.set_xlabel("Tissue Surface Area")
ax1.set_ylabel(r"Mean $\kappa$ Curvature")
#ax1.set_ylim(-0.2,0.)
###################################
# Height of Primodia
###################################
ax2.set_title("Height of Primordia")
ax2.set_xlabel("Tissue Surface Area")
ax2.set_ylabel("Primordial Height")
#########################################
# Primordia area vs Total surface area
#########################################
ax3.set_title("Primordia area vs Surface Area")
ax3.set_xlabel("Tissue Surface Area")
ax3.set_ylabel("Primodia area")
###################################
# Surface area vs time
###################################
ax4.set_title("Tissue Surface area growth")
ax4.set_xlabel("time")
ax4.set_ylabel("Tissue Surface Area")
#################################################################################
fig1 = plt.figure(figsize=(20,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax31 = fig1.add_subplot(121)
ax41 = fig1.add_subplot(122)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
ax31.set_title("Mean Gaussian Curvature")
ax31.set_xlabel(r"$\frac{AreaOfPrimodia}{InitialAreaOfPrimodia}$")
ax31.set_ylabel("Mean Gaussian Curvature")
#ax1.set_ylim(-0.2,0.)
###################################
# Height of Primodia
###################################
ax41.set_title("Height of Primodia")
ax41.set_xlabel(r"$\frac{AreaOfPrimodia}{InitialAreaOfPrimodia}$")
ax41.set_ylabel("Height of Primodia")
#################################################################################
fig2 = plt.figure(figsize=(20,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax5 = fig2.add_subplot(121)
ax6 = fig2.add_subplot(122)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
ax5.set_title("Mean Gaussian Curvature")
ax5.set_xlabel(r"$\frac{AreaOfPrimodia}{SurfaceAreaOfTissue}$")
ax5.set_ylabel("Mean Gaussian Curvature")
#ax1.set_ylim(-0.2,0.)
###################################
# Height of Primodia
###################################
ax6.set_title("Height of Primodia")
ax6.set_xlabel(r"$\frac{AreaOfPrimodia}{SurfaceAreaOfTissue}$")
ax6.set_ylabel("Height of Primodia")
#ax2.set_ylim(0.1,1.)
########################################################
fig3 = plt.figure(figsize=(10,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax7 = fig3.add_subplot(111)
ax7.set_xlabel("time")
ax7.set_ylabel("sphericity")
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
	plotData[etacurrent] = plotMinGaussianCurvaturePrimodiaHeight(numOfLayer = numOfLayer, targetid = targetid,endStep = endStep,eta = etacurrent,
				mincurvatureplot = ax1,startStep = startStep,  
                heightplot=ax2,ax3= ax3, ax4 = ax4, ax5 = ax5, ax6 = ax6,ax7 = ax7,
				color = etacolor,stepsize = stepsize,
                largerCondition = large,maxarea = maxarea)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#color bar fig
scalarMap._A = []
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	clrbar.set_label(r"Fast Growth Rate")
else:
	clrbar.set_label(r"$\eta$")


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
	fig.savefig(saveDirectory+r"/plot_fk_vs_curvature_height_%d_targetface=%d.png"%(endStep,targetid),transparent = True)
else:
	fig.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_%d_targetface=%d.png"%(endStep,targetid),transparent = True,)



#fig1.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_areaPrimodia_%d.png"%endStep,transparent = True)
#fig2.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_surfaceratio_%d.png"%endStep,transparent = True)
#fig3.savefig(saveDirectory+r"/plot_eta_vs_sphericity_%d.png"%endStep,transparent = True)
plt.close('all')
### Saving Data Dictionary ###
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	np.save('meancurvature_height_fk_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)
else:
	np.save('meancurvature_height_eta_time=%d_targetface=%d.npy'%(endStep,targetid),plotData)


################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"

