import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append("/Users/jasonkhadka/Documents/git/plantdev")
sys.path.append("/home/jkhadka/plantdev/")
import quadedge as qd
sys.path.append('/Users/jasonkhadka/Documents/git/plantdev/python_quadedge')
sys.path.append('/home/jkhadka/plantdev/python_quadedge')
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
sys.path.append('/home/jkhadka/transferdata/scripts/strain_plots/')
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
from mpl_toolkits.mplot3d import Axes3D
#from simulation_functions import *
import simulation_functions as sf
import ellipse as ep
import os

#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20.
plt.rcParams['ytick.labelsize'] = 20.
plt.rcParams['axes.labelsize'] = 25.
plt.rcParams['legend.fontsize'] = 25.
plt.rcParams['axes.titlesize'] = 30.

##########################################################################################
#  Getting points of target form matrix
##########################################################################################
def getTargetFormEllipsePoints(cell,targetface = 10):
    #getting the Target Form Matrix
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while True:
        if face.getID() == targetface:
            break
        face = faces.next()
    #print "Face Id : ", face.getID()
    targetformmatrix = np.array([[qd.getTargetFormMatrix(face, 0,0),qd.getTargetFormMatrix(face, 0,1)],
                                 [qd.getTargetFormMatrix(face, 1,0),qd.getTargetFormMatrix(face, 1,1)]
                                 ])
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
    data = ep.plot_ellipse(cov=targetformmatrix, data_out=True)
    points = np.matrix(np.vstack((data,np.zeros(len(data[0])))))
    transformedpoints = transpose_unitmat*points
    transformedpoints[0]+= xcent
    transformedpoints[1]+= ycent
    transformedpoints[2]+= zcent
    ################################
    return transformedpoints
################################################################################################
def getCurrentFormEllipsePoints(cell,targetface = 10):
    #getting the Target Form Matrix
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while True:
        if face.getID() == targetface:
            break
        face = faces.next()
    #print "Face Id : ", face.getID()
    targetformmatrix = np.array([[qd.getCurrentFormMatrix(face, 0,0),qd.getCurrentFormMatrix(face, 0,1)],
                                 [qd.getCurrentFormMatrix(face, 1,0),qd.getCurrentFormMatrix(face, 1,1)]
                                 ])
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
    data = ep.plot_ellipse(cov=targetformmatrix, data_out=True)
    points = np.matrix(np.vstack((data,np.zeros(len(data[0])))))
    transformedpoints = transpose_unitmat*points
    transformedpoints[0]+= xcent
    transformedpoints[1]+= ycent
    transformedpoints[2]+= zcent
    ################################
    return transformedpoints
################################################################################################
def getStrainEllipsePoints(cell,targetface = 10):
    #getting the Target Form Matrix
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while True:
        if face.getID() == targetface:
            break
        face = faces.next()
    #print "Face Id : ", face.getID()
    targetformmatrix = ep.getStrainMatrix(cell,targetface)
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
    data = ep.plot_ellipse(cov=targetformmatrix, data_out=True)
    points = 2.*np.matrix(np.vstack((data,np.zeros(len(data[0])))))
    transformedpoints = transpose_unitmat*points
    transformedpoints[0]+= xcent
    transformedpoints[1]+= ycent
    transformedpoints[2]+= zcent
    ################################
    return transformedpoints
###################################################################
def plotStrainEllipseSurface(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0,save = False, ids = True, azim = -60, elev = 60, targetface = None):
    #calculating forces, stress-matrix and strain-matrix
    #cell.calculateVertexForce()
    cell.calculateStressStrain()
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    #import the libraries
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib as mpl
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np
    import matplotlib.pyplot as plt
    #limits of the plot
    radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
    #plotting part
    fig = plt.figure(frameon=False,figsize=(20,16))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    #fig.set_aspect(aspect='equal', adjustable='box')
    ax = Axes3D(fig)
    ax.set_xlim((-.5*radius,.5*radius))
    ax.set_ylim((-.5*radius,.5*radius))
    ax.set_zlim((-0.,1.*radius))
    ax.axis('off')
    #ax.axis('equal')
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ########################################################################
    #                 Plotting the Cell                    #
    ########################################################################
    faces = qd.CellFaceIterator(cell)
    ###################
    face = faces.next()
    while (face != None):
        if face.getID() == 1:
            face  = faces.next()
            continue
        faceid = face.getID()#grabbing face id
        xlist = []
        ylist = []
        zlist = []
        #print "== Face ID : ", faceid, "=="
        xmean = face.getXCentralised()
        ymean = face.getYCentralised()
        zmean = face.getZCentralised()
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
        ax.plot(xlist,ylist,zlist,'k')
        if ids:
                color = 'k'
                ax.text(xmean,ymean,zmean, face.getID(),fontsize = 20, color = color)
        if faceid == targetface:## if this face is the target face
                color = 'r'
                ax.text(xmean,ymean,zmean, face.getID(),fontsize = 20, color = color)
        face = faces.next()
        #if face.getID() == 1: break
    #plt.clf()
    ##PLOTTING ELLIPSE
    ##PLOTTING ELLIPSE
    numofface = cell.countFaces()
    for i in range(2,numofface+1):
        # Target Form
        TFellipsepoints = getTargetFormEllipsePoints(cell,i)
        verts = [zip(np.array(TFellipsepoints[0])[0],np.array(TFellipsepoints[1])[0], np.array(TFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = 0.5,linewidths=1, facecolor = 'r',zorder = 10,label="TFM" if i == 0 else "")
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # Current Form
        CFellipsepoints = getCurrentFormEllipsePoints(cell,i)
        verts = [zip(np.array(CFellipsepoints[0])[0],np.array(CFellipsepoints[1])[0], np.array(CFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = 0.5,linewidths=1,zorder = 9, facecolor = 'b')
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # Strain 
        strainellipsepoints = getStrainEllipsePoints(cell,i)
        verts = [zip(np.array(strainellipsepoints[0])[0],np.array(strainellipsepoints[1])[0], np.array(strainellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = 0.5,linewidths=1,zorder = 11, facecolor = 'k')
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
    ax.view_init(azim = azim,elev=elev)
    scatter1_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='r', marker = 's',ms=20)
    scatter2_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='b', marker = 's',ms=20)
    scatter3_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='k', marker = 's',ms=20)
    ax.legend([scatter1_proxy, scatter2_proxy,scatter3_proxy], ['Target Form', 'Current Form',"Strain"], numpoints = 1,loc = 0)
    plt.suptitle("Step =%03d"%step,fontsize = 30)
    ### Saving
    if save:
        saveDirectory = "strainFigures"
        import os
        if not os.path.exists(saveDirectory):
            os.makedirs(saveDirectory)
        plt.savefig(saveDirectory+r"/strainEllipse_Time=%03d.png"%(step))
        plt.close()
    return

############################################################################################################
################################################################################################
#			Plotting Part 
################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =0, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-i","--ids", help = "if option is used, ids on the faces will be, else by default it will not", action= "store_true")

parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0.8",
                                   default = .8, type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
                                   default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
                                   default = 0.0, type = float)
parser.add_argument("-n","--angle",help = "value to set for convex angle threshold, default = 360",
                                   default = 360., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)
parser.add_argument("-t","--target", help = "Target face for faster growth", default = 0, type = float)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)

## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
ids = args.ids
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
anglethreshold = args.angle
targetface = args.target
#################################################################################
import sys
for step in range(startStep,endStep+1):
	percentStep = int((step-startStep)/float(endStep-startStep)*100)
	sys.stdout.write('\r'+"step : "+ str(step) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
	sys.stdout.flush()
	cell = qd.objReadCell("qdObject_step=%03d.obj"%step)
	# Flipping the Face ID after Loading to cell so that the face id before and after matches
	faces = qd.CellFaceIterator(cell)
	facecount = cell.countFaces()
	face= faces.next()
	while face != None:
	    faceid = face.getID()
	    face.setID(facecount-faceid+1)
	    #print face.getID()
	    face = faces.next()
	######################################################
	#settig target Form Matrix : 
	#TMF step corresponds to coordinate step
	####
	cell = sf.setTargetFormMatrix(cell,step)
	cell.setParameters()
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStrain()
	#print " Energy: ", cell.getEnergyCartesianVolume()
	###########################################################
	#plotting and Savnig
	############################################################
	plotStrainEllipseSurface(cell,numOfLayer,step = step, alpha = alpha, save = True,ids = ids, 
		targetface = targetface,azim = args.azimuthal, elev = args.elevation)
	latdev.plotSurface(cell,numOfLayer,name = "strainFigures/step_%03d.png"%step, alpha = alpha, azim = args.azimuthal, elev = args.elevation)

################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"
