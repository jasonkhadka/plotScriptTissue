import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/jasonkhadka/Documents/git/plantdev")
sys.path.append('/home/jkhadka/plantdev')
sys.path.append('/home/jkhadka/plantdev/python_quadedge')
sys.path.append('/home/jkhadka/transferdata/scripts/strain_plots/')
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
import quadedge as qd
#sys.path.append('/Users/jasonkhadka/Documents/git/plantdev/python_quadedge')
#sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scripts/strain_plots/')
#sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scripts/simulation_functions/')
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import simulation_functions as sf
import ellipse as ep

###############################################################################################################
###         Loading the TargetFormMatrix from file and setting it to the faces for given growthstep
###############################################################################################################
def setTargetFormMatrix(cell, growthcounter):
    ### Loading the TargetFormMatrix
    loadedDictionary = np.load("TargetFormMatrix_step=%d.npy"%growthcounter).item()
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face != None:
        faceid = face.getID()
        matrixArray = loadedDictionary[faceid]
        # Setting Target Form Matrix
        qd.setTargetFormMatrix(face, 0,0, matrixArray[1][1])
        qd.setTargetFormMatrix(face, 1,1, matrixArray[0][0])
        qd.setTargetFormMatrix(face, 0,1, matrixArray[0][1])
        qd.setTargetFormMatrix(face, 1,0, matrixArray[1][0])
        face = faces.next()
    cell.setParameters()
    ######
    return cell
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
def getStrainEllipsePoints(cell,targetface = 10,factor = 1):
    #getting the Target Form Matrix
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while True:
        if face.getID() == targetface:
            break
        face = faces.next()
    #print "Face Id : ", face.getID()
    targetformmatrix = factor*ep.getStrainMatrix(cell,targetface)
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
# Check if a cell is growing fast or not
################################################################################################
def checkFastGrowing(targetid, i):
    fastGrowing = [targetid]
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while (face != None):
            if face.getID() == targetid:
                    edges = qd.FaceEdgeIterator(face)
                    edge = edges.next()
                    while edge != None:
                            rightFace = edge.Right()
                            fastGrowing.append(rightFace.getID())
                            edge = edges.next()
                    #print "kappa for this face : ", face.getKappa()
            face = faces.next()
    if i in fastGrowing:
        return True
    return False
###############################################################################################
# Make a target cells and its neighbour fast growing
###############################################################################################
def fastGrowth(cell,faceid,fastkappa=2.0,fastgrowthvar=0.2):
        faces = qd.CellFaceIterator(cell)
        face = faces.next()
        while (face != None):
                if face.getID() == faceid:
                        #print "Setting Fast Kappa Growth for Face ID :", faceid
                        #print face.getKappa()
                        face.setKappa(fastkappa)
                        face.setGrowthVar(fastgrowthvar)
                        edges = qd.FaceEdgeIterator(face)
                        edge = edges.next()
                        while edge != None:
                                rightFace = edge.Right()
                                rightFace.setKappa(fastkappa)
                                rightFace.setGrowthVar(fastgrowthvar)
                                edge = edges.next()
                        #print "kappa for this face : ", face.getKappa()
                face = faces.next()
        return cell
def loadCellFromFile(step,numOfLayer = 8,eta = 1.,kappa = 0.2,fastkappa = 2.0, targetid = 135):
    ####################################
    loadedcell = qd.objReadCell("qdObject_step=%03d.obj"%step)
    loadedcell.setInitialParameters()
    # Flipping the Face ID after Loading to cell so that the face id before and after matches
    faces = qd.CellFaceIterator(loadedcell)
    facecount = loadedcell.countFaces()
    face= faces.next()
    while face != None:
        faceid = face.getID()
        face.setID(facecount-faceid+1)
        #print face.getID()
        face = faces.next()
    ######################################################
    #print "######################################################"
    #print "#"," "*10, "step %d"%step
    #print "######################################################"
    #settig target Form Matrix : 
    #TMF step corresponds to coordinate step
    ####
    sf.setTargetFormMatrix(loadedcell,step)
    #### CELL PARAMETERS ####
    loadedcell.setKappa(kappa)
    loadedcell.setEta(eta)
    #setTargetFormMatrix(loadedcell,step)
    loadedcell.setParameters()
    #calculating forces, stress-matrix and strain-matrix
    loadedcell.calculateStressStrain()
    loadedcell = fastGrowth(loadedcell,targetid,fastkappa = fastkappa)
    ######################################################
    #latdev.plotSurface(loadedcell,numOfLayer,name="dome_remake_%03d.png"%step,ids=False,alpha = 1.,elev = 7)
    #plotTargetFaceArrayGrowSurface(loadedcell,numOfLayer)
    return loadedcell

def plotTargetGrowSurface(cell, numOfLayer, step = None, targetid = 135,
                            ids = False,save = False, alpha = 0.8, Length=1.0,azim = -60.,elev=60,figsize = (20,16)):
    #calculating forces, stress-matrix and strain-matrix
    cell.calculateVertexForce()
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
    fig = plt.figure(frameon=False,figsize=figsize)
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
    faceids = []
    while (face != None):
        if face.getID() ==1:
            face  = faces.next()
            continue
        faceid = face.getID()#grabbing face id
        xlist = []
        ylist = []
        zlist = []
        faceids.append(faceid)
        #print "== Face ID : ", faceid, "=="
        xmean = face.getXCentralised()
        ymean = face.getYCentralised()
        zmean = face.getZCentralised()
        #print faceid
        ###############################################
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
        #ax.plot(xlist,ylist,zlist,'k',lw = 2)
        verts = [zip(np.array(xlist),np.array(ylist),np.array(zlist))]
        #pc = Poly3DCollection(verts,alpha = 0.1,linewidths = 4,zorder = 1, facecolor = 'w')
        #pc.set_edgecolor('k')
        #ax.add_collection3d(pc)
        ax.plot(xlist,ylist,zlist,'k',lw = 3)
        if ids:
            ax.text(xmean,ymean,zmean, face.getID(),fontsize = 20)
        face = faces.next()
        #if face.getID() == 1: break
    #plt.clf()
    ##PLOTTING ELLIPSE
    #numofface = cell.countFaces()
    #faceids = getNeighbouringFaces(cell,targetid)
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face != None:
        #face id : 
        i = face.getID()
        if i == 1:
            face = faces.next()
            continue
        # Strain 
        strainellipsepoints = getStrainEllipsePoints(cell,i,factor=4)
        verts = [zip(np.array(strainellipsepoints[0])[0],np.array(strainellipsepoints[1])[0], np.array(strainellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = 1.,linewidth=0,zorder = 100, facecolor = 'k')
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # Target Form
        TFellipsepoints = getTargetFormEllipsePoints(cell,i)
        verts = [zip(np.array(TFellipsepoints[0])[0],np.array(TFellipsepoints[1])[0], np.array(TFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = .7,linewidth=0, facecolor = 'r',zorder = 10,label="TFM" if i == 0 else "")
        pc.set_edgecolor('k')
        #ax.add_collection3d(pc)
        # Current Form
        CFellipsepoints = getCurrentFormEllipsePoints(cell,i)
        verts = [zip(np.array(CFellipsepoints[0])[0],np.array(CFellipsepoints[1])[0], np.array(CFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = 0.5,linewidth=0,zorder = 5, facecolor = 'g')
        pc.set_edgecolor('k')
        #ax.add_collection3d(pc)
        #######################
        face = faces.next()
    ### Now do growth on the target faces
    cell = sf.feedbackStrainGrowFaces(cell)
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face != None:
        #face id : 
        i = face.getID()
        if i == 1:
            face = faces.next()
            continue
        if checkFastGrowing(targetid, i):
            facecolor = 'g'
        else:
            facecolor = 'm'
        # Target Form
        TFellipsepoints = getTargetFormEllipsePoints(cell,i)
        verts = [zip(np.array(TFellipsepoints[0])[0],np.array(TFellipsepoints[1])[0], np.array(TFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = .5,linewidth=0, facecolor = facecolor,zorder = 10,label="TFM" if i == 0 else "")
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        face = faces.next()
    ax.view_init(azim = azim,elev=elev)
    #scatter1_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='r', marker = 's',ms=50)
    scatter2_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='g', marker = 's',ms=50)
    scatter3_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='k', marker = 's',ms=50)
    scatter4_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='m', marker = 's',ms=50)
    ax.set_title("Feedback Growth and strain visualised step %d"%step)
    ax.legend([ scatter2_proxy,scatter3_proxy,scatter4_proxy], ['Target Shape', 'Fast Growing',"Strain","Normal Growing"],fontsize = 30, numpoints = 1,loc = 0)
    if save:
        saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/ellipsePlot"
        import os
        if not os.path.exists(saveDirectory):
            os.makedirs(saveDirectory)
        plt.savefig(saveDirectory+r"/ellipsePlot_Time=%03d.png"%(step),transparent=True)
        plt.close()
    #ax.text(0,0,0,'P',fontsize = 50,zorder=-1)
    #plt.suptitle("Step =%03d"%step,fontsize = 30)
    #plt.savefig('initial_TFM_layer8.png', transparent=True)
    return




################################################################################################
#           Plotting Part 
################################################################################################
import argparse #argument parser, handles the arguments passed by command line

#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--start", help="Start of simulation step",default =1, type = int)
parser.add_argument('-e',"--end", help="End of simulation step", type = int)
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument('-l', "--layer", help = "The number of layers in the quadedge cell",type=int)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
                                   default = 1., type = float)
parser.add_argument("-k","--kappa",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
                                   default = .2, type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
                                   default = 1., type = float)
parser.add_argument("-q","--fastkappa",help = "Fast kappa, to grow a patch of cells fast, default = 0.0",
                                                                     default = 0.0, type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
                                   default = 0.0, type = float)
parser.add_argument("-n","--eta",help = "feedback parameter, default = 0",
                                   default = 0., type = float)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)

parser.add_argument("-t","--target", help = "Target face for faster growth", default = 10, type = float)
parser.add_argument("-u","--azimuthal", help = "azimuthal angle for display", default = -60, type = float)
parser.add_argument("-v","--elevation", help = "elevation angle for display", default = 60, type = float)

## Getting the arguments 
args = parser.parse_args()
#location = args.location
endStep = args.end
startStep = args.start
cylinder = args.cylinder
alpha = args.alpha
kappa = args.kappa
fastkappa = args.fastkappa
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
eta = args.eta
targetface = args.target
azim = args.azimuthal
elev = args.elevation

#################################################################################
import sys
import os
for step in range(startStep,endStep+1):
    percentStep = int((step-startStep)/float(endStep-startStep)*100)
    sys.stdout.write('\r'+"step : "+ str(step) +" "+"#"*percentStep+' '*(100-percentStep)+"%d%%"%percentStep)
    sys.stdout.flush()
    cell = loadCellFromFile(step,kappa = kappa,eta = eta, fastkappa = fastkappa)
    #print " Energy: ", cell.getEnergyCartesianVolume()
    ###########################################################
    #plotting and Savnig
    ############################################################
    #plotStrainMagnitude(cell,numOfLayer,step = step, targetface = targetface,save = True,azim = azim, elev = elev)
    #plotStrainSurface(cell,numOfLayer,step = step, save = True,azim = azim, elev = elev)
    #plotFaceArea(cell,numOfLayer,step = step, save = True,azim = azim, elev = elev)
    # Directory to save the figures
    saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/ellipsePlot"
    if not os.path.exists(saveDirectory):
        os.makedirs(saveDirectory)
    plotTargetGrowSurface(cell,numOfLayer,step = step, azim = azim , elev = elev, save = True)
    #Wlatdev.plotSurface(cell,numOfLayer,name=saveDirectory+r"/surface_%03d.png"%step,alpha = 0.8, azim = azim, elev = elev)
    plt.close("all")
################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

