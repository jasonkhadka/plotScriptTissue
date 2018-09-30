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
plt.rcParams['xtick.labelsize'] = 10.
plt.rcParams['ytick.labelsize'] = 10.
plt.rcParams['axes.labelsize'] = 10.
plt.rcParams['legend.fontsize'] = 10.
plt.rcParams['axes.titlesize'] = 10.

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
if ploteta:
		fklist = [float(dict(item.split("=") for item in folder.split("_"))['fk']) for folder in listdir]
else:
		fklist = [float(dict(item.split("=") for item in folder.split("_"))['n']) for folder in listdir]
#################################################################################
#     Making the Plots
#################################################################################
import matplotlib.colors as colors
import matplotlib.cm as cmx
jet = cm = plt.get_cmap('viridis') 
maxvalue = max(fklist)
minvalue = min(fklist)
#print "Max value", maxvalue, " minvalue", minvalue
cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#fig = plt.figure(frameon=False,figsize=(20,16))
fig = plt.figure(figsize=(10,10))
#fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
#ax2 = fig.add_subplot(122)
#fig.set_aspect(aspect='equal', adjustable='box')
#ax.axis('equal')
#################################################################################
# Min Gaussian curvature
##################################
ax1.set_title("Face Area")
ax1.set_xlabel("Time")
ax1.set_ylabel("Face Area")
##################################
ax2.set_title("Mean Face Area")
ax2.set_xlabel("Time")
ax2.set_ylabel("Mean Face Area")
##################################
ax3.set_title("Log-Log plot Area vs time")
ax3.set_xlabel(r"$log(t)$")
ax3.set_ylabel(r"$log(<A>)$")
#ax1.set_ylim(0.,80.)
###################################
# Height of Primodia
###################################
#ax2.set_title("Height of Primodia")
#ax2.set_xlabel("time")
#ax2.set_ylabel("Height of Primodia")
########################################################
if not saveall:
	counter = 0
	totalfolders = len(listdir)
	#print listdir
	for folder in listdir:
		# Converting folder name to dictionary
		#print folder
		if ploteta:
			fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
		else:
			fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
		#scaled color for plotting
		#print fkcurrent
		fkcolor = scalarMap.to_rgba(fkcurrent)
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
		plotAverageFaceArea(endStep,ax1,ax2,ax3,color = fkcolor,startstep=1,norm=True,fastid = targetface)
		#plotMeanGrowthRate(step=endStep,azim = azim, elev = elev,eta = float(folderdict['n']),save=True)
		os.chdir("..")
		counter+= 1
	###############################################################################
	#color bar
	scalarMap._A = []
	fig.subplots_adjust(bottom=0.2)
	cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.02])
	clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)
	if ploteta:
		clrbar.set_label("Fast Kappa")
	else:
		clrbar.set_label("$\eta$")
	#plt.tight_layout()
	#plt.tight_layout( rect=[0, 0, 1, 1])
	plt.savefig(saveDirectory+r"/plot_face_area.png",transparent = True)
	plt.close('all')
else:
	counter = 0
	ax1.set_xlabel("Time")
	ax1.set_ylabel("Face Area")
	#ax1.set_ylim(0.,80.)
	totalfolders = len(listdir)
	#print listdir
	for folder in listdir:
		# Converting folder name to dictionary
		#print folder
		if ploteta:
			fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['fk'])
		else:
			fkcurrent = float(dict(item.split("=") for item in folder.split("_"))['n'])
		#scaled color for plotting
		#print fkcurrent
		fkcolor = scalarMap.to_rgba(fkcurrent)
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
		plotAverageFaceArea(endstep = endStep,areaplot = ax1,ax2 = ax2, color = fkcolor,startstep=1,norm=True,fastid = targetface)
		if ploteta:
			plt.savefig(saveDirectory+r"/plot_face_area_fk=%.2f.png"%fkcurrent,transparent = True)
		else:
			ax1.set_title("Face Area $\eta$ = %.2f"%fkcurrent)
			plt.savefig(saveDirectory+r"/plot_face_area_n=%.2f.png"%fkcurrent,transparent = True)
		plt.cla()
		#plotMeanGrowthRate(step=endStep,azim = azim, elev = elev,eta = float(folderdict['n']),save=True)
		os.chdir("..")
		counter+= 1
plt.close('all')

################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

