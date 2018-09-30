####################################################################################
#                       CLUSTER CODE 
#            Task :
#                  To check time taken by all algorithms to perform relaxation
#
#            Initial condition :
#                        Uniform Dome
#            Bounds : 
#                   Tight bounds for all vertices to move, only allowed to move in a cube of 0.5 Length units
#                   basecells are fixed in their position
#            Growth : 
#                   No growth 
#            Termination Conditon : 
#                       Run the relaxation until the tolerance is reached or Maxtime is acheived
#            Algorithm :
#                       All of them are tested                            #
####################################################################################
#importing all the libraries
#importing all the libraries
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
#pd.set_option('precision',10)
#pd.set_option("display.max_columns",200)
#pd.set_option("display.max_rows",2000)
import sys 
sys.path.append("/home/jkhadka/plantdev")
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
import quadedge as qd
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
#from IPython.display import display, HTML
import nlopt #non linear optimizer
import argparse #argument parser, handles the arguments passed by command line
import os # to make directory
import time
#importing functions that the simulations rely on -- make sure to copy simulation_functions.py to the same directory as this file
from simulation_functions import *
from datetime import datetime
#timing 
##################################################################################################################################################
##################################################################################################################################################
programstart = time.time()
print "This is a test for Growth !!"
algorithm = "SBPLX"
print "The Used Algorithm :: ", algorithm
################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#adding arguments
parser.add_argument("-f","--flat", help = "with or without dome, makes it flat if argument is used", action= "store_true")
parser.add_argument("-u","--unlimited", help = "for unlimited time simulation until tolerance is reached", action= "store_true")

parser.add_argument("-e","--error",help = "value for relative error tolerance, default = 1e-4",
                                   default = 10**(-5), type = float)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
                                   default = 0., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
                                   default = 0.001, type = float)
parser.add_argument("-k","--kappa",help = "value to set for kappa growth Rate of Mo, default = 1.0",
                                   default = 1., type = float)
parser.add_argument("-t","--time",help = "total steps for relaxation to run, default = 3000",
                                   default = 1000, type = int)
parser.add_argument("-i","--interval",help = "time interval for saving the plot of tissue (other datas if applicable)",
                                   default = 30, type = int)
parser.add_argument("-l","--layer",help = "number of layers of cells to simulate",
                                   default = 5, type = int)
parser.add_argument("-s","--shape",help="Initial target shape, if used makes target shape to be Identity Matrix, else by default it is the current initial shape", action = "store_true")
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = 0.1, type = float)
parser.add_argument("-v","--varypressure",help="vary the pressure as multiple of 0.00001*(10)**(passedvalue)",default = 0, type = int)
parser.add_argument("-j","--jobname",help="jobname if given, the new folder to be created will be named jobname", default = None, type = str)
#parser.add_argument("-n","--numberofimage",help = "number of images or data points to store during simulation", default = 10, type = int)
##################################################################################################################################################
##################################################################################################################################################
#Setting Parameters
args = parser.parse_args()#arguments passed
#gathering the values passed by argument
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
totaltime = args.time
interval = args.interval
layer = args.layer
gamma = args.gamma
kappa = args.kappa
tolerance = 10.**(-4)
initialstep = 10.**(-2)
growthSteps = 1000 #performing 1000 growth steps
jobname = args.jobname
#chaning pressure if the Vary pressure is used
if args.varypressure != 0:
    gamma = 0.00001*10**(args.varypressure)
#making of hexagonal lattice
numOfLayer = layer
#geonetrical parameters of the tissue
Length = 1.# length of the sides of the faces
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on

numofface = 3*numOfLayer*(numOfLayer-1)+1
#number of vertices in the perimeter of lattice
perimeterVertexNum = 6+12*(numOfLayer-1)
perimeterHexNum = (numOfLayer!=1)*(6*(numOfLayer-1)-1) + 1 # number of perimeter hexgons
#perimeterVertexNum = 6+12*(numOfLayer-1)#number of vertices in the perimeter
#######################################################################################################################################
#               Making Directory
#######################################################################################################################################
###Making the directory to store the image files ####
#making Main directory
#making a directory to save coordinates after each relaxation step before growth kicks in 
#######################################################################################################################################
#               Making Directory
#######################################################################################################################################
###Making the directory to store the image files ####
#making Main directory
if not os.path.exists("m0=m"):
     os.makedirs("m0=m")
os.chdir("m0=m")
"""
#with_dome directory
if args.flat: #if flat option is used : the tissue will be grown as flat ; NOT dome
    withoutdome = "without_dome"
    if not os.path.exists(withoutdome):
         os.makedirs(withoutdome)
    os.chdir(withoutdome)
else:
    withdome = "with_dome"
    if not os.path.exists(withdome):
        os.makedirs(withdome)
    os.chdir(withdome)
"""
############################################################################################################
#if the argument -u is used then the loop is infinite until the convergence is reached with given tolerance
##############################################################################################################
if args.unlimited:
    totaltime = sys.maxint
    newpath = r'a=%(alpha).4f_b=%(beta).4f_g=%(gamma).4f_p=%(pressure)f_z=%(zeta)f_t=%(time)s_layer=%(layer).1d'%{"alpha":alpha,"beta":beta,"gamma":gamma,"pressure":pressure,"zeta":zeta,"time":"unlimited","layer":numOfLayer}
else:
    #now making path for simulation result
    newpath = r'a=%(alpha).4f_b=%(beta).4f_g=%(gamma).4f_p=%(pressure).4f_z=%(zeta).2f_t=%(time).4f_layer=%(layer).1d'%{"alpha":alpha,"beta":beta,"gamma":gamma,"zeta":zeta,"pressure":pressure,"time":totaltime,"layer":numOfLayer}
    #newpath = r'M0=identity_a=%(alpha)f_b=%(beta)f_p=%(pressure)f_t=%(time)f_layer=%(layer)d'%{"alpha":alpha,"beta":beta,"pressure":pressure,"time":totaltime,"layer":numOfLayer}

if not os.path.exists(newpath):
    os.makedirs(newpath)
#change directory to newpath
os.chdir(newpath)
#making a directory to save coordinates after each relaxation step before growth kicks in 
growthdirectory = 'growth_images'
if not os.path.exists(growthdirectory):
    os.makedirs(growthdirectory)
################################################################################################
##########################################################################################################
#                                SETTING INITIAL CONDITION 
##########################################################################################################

#makeCenteredHexagonLattice(numOfLayer, printcondition = False, projection = True) 
# Projection condition : it is to have the external vertex projected on to an external circle
cell, newedge = latdev.makeCenteredHexagonLattice(numOfLayer, True, False) 
###### put down the paramters to cell
cell.setAlpha(alpha)
cell.setBeta(beta)
cell.setPressure(pressure)
cell.setZeta(zeta)
cell.setGamma(gamma)
cell.setKappa(kappa)
#########    Printing Parameters    ########################
print "Layer Number : ", numOfLayer
print "Number of Faces : ", cell.countFaces()
print "Calculated : Radius = ", radius
print "Parameters"
print "------------------------"
print "alpha : ", cell.getAlpha()
print "beta : ", cell.getBeta()
print "pressure : ", cell.getPressure()
print "zeta : ", cell.getZeta()
print "Kappa : ", cell.getKappa()
print "tolerance : ", tolerance
print "gamma :",  gamma
###################################################################
#latdev.plot3DCell(cell)#plotting the funciton
#print "number of vertices in the boundary", perimeterVertexNum
############################################################################################
##################################################################
######## Making the DOME
##################################################################
#print "radius : ", radius
#zfactor = radius
#hemisphereradius = radius/np.arccos(0)
hemisphereradius = radius/1.5
#zfactor = 0.5*np.pi*radius - radius
#print hemisphereradius#, zfactor
zmin = 10
#now converting each vertex to spherical surface
vertices = qd.CellVertexIterator(cell)
vertex = vertices.next()
counter = 0#counter to change vertex positions
while vertex != None:
    x = vertex.getXcoordinate()
    y = vertex.getYcoordinate()
    z = vertex.getZcoordinate()
    #####scaling the coordinates to the spherical surface
    l = np.sqrt(x*x+y*y)
    a=l/hemisphereradius
    z=hemisphereradius*np.cos(a)
    if z>= 0.0:
        a = hemisphereradius*np.sin(a)/l
    else:
        z = 0.5*np.pi*hemisphereradius-l
        a = hemisphereradius/l
    x*=a
    y*=a
    vertex.setXcoordinate(x)#setting x coordinate
    vertex.setYcoordinate(y)#setting y coordinate
    vertex.setZcoordinate(z)#setting z coordinate
    if zmin > z:
        zmin = z
    counter += 1
    vertex = vertices.next()
#print zmin
###resaling vertices To make it positive
vertices = qd.CellVertexIterator(cell)
vertex = vertices.next()
counter = 0#counter to change vertex positions
while vertex != None:
    z = vertex.getZcoordinate()
    vertex.setZcoordinate(z-zmin)
    vertex = vertices.next()
#####################################################
cell = settingFirstParameters(cell)#setting all the parameters of the cell
#####################################################
latdev.plot3DCell(cell,name = r'initial_dome.png')
latdev.plotSurface(cell, numOfLayer,name = r'initial_dome_surface.png')
#########################################################################################
####                    Putting in the relaxed shape of the dome                      ###
#########################################################################################
relaxed_coordinates = np.load("/home/jkhadka/relaxed_coordinates/equidistant_hexagon/layer8/NxFrB8_S_z10g="+str(gamma)+"/relaxed_coordinates.npy")
relaxed_coordinates = relaxed_coordinates.reshape((3,int(len(relaxed_coordinates)/3)))
print "loading vertices !!!! "
########################################################################################################################
vertices = qd.CellVertexIterator(cell)
vertex = vertices.next()
counter = 0#counter to change vertex positions
while vertex != None:
    vertex.setXcoordinate(relaxed_coordinates[0,counter])#setting x coordinate
    vertex.setYcoordinate(relaxed_coordinates[1,counter])#setting y coordinate
    vertex.setZcoordinate(relaxed_coordinates[2,counter])#setting z coordinate
    counter += 1
    vertex = vertices.next()
print " vertices loaded !!!"
#########################################################################################
cell = settingParameters(cell)
#########################################################################################
#######################################################################################
####        IMPORTANT FUCNTIONS FOR SIMULATION                                     ####
#######################################################################################
############################################################################
###                 objective function for NLOPT                         ###
############################################################################
def energyobjective(inputcoordinates,grad):
    global maincounter, numOfLayer, alpha, beta, pressure, cell
    #making of hexagonal lattic
    #Reshape the tempcoordinates, to x-y-z arrays
    tempcoordinates = inputcoordinates.reshape((3,numberofvertices))
    ####iterating the vertices to feed the new coordinates in
    vertices = qd.CellVertexIterator(cell)
    vertex = vertices.next()
    counter = 0#counter to change vertex positions
    while vertex != None:
        vertex.setXcoordinate(tempcoordinates[0,counter])#setting x coordinate
        vertex.setYcoordinate(tempcoordinates[1,counter])#setting y coordinate
        vertex.setZcoordinate(tempcoordinates[2,counter])#setting z coordinate
        counter += 1
        vertex = vertices.next()
    #####################################################
    cell = settingParameters(cell)
    ######################################################
    #returning the total energy
    energyvalue = cell.getEnergyCartesianVolume()
    #release the current cell
    #printing the counter and if the passed value is different than the previous values passed
    #print maincounter, energyvalue
    #print np.subtract(tempcoordinates, tempcoordstore)
    return energyvalue
###################################################################################################################################
############################################################################
###            Function to print area of faces                           ###
############################################################################
def printareaofface():
    global cell
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face != None:
        print "Face ID : ", face.getID(), "Area : ", face.getAreaOfFace()
        face = faces.next()
    return

###################################################################################################################################
############################################################################
###                    Return Energy Values                              ###
############################################################################
def termsofenergy(inputcoordinates):
    global maincounter, numOfLayer, alpha, beta, pressure, cell
    #making of hexagonal lattice
    maincounter += 1 
    #Reshape the tempcoordinates, to x-y-z arrays
    tempcoordinates = inputcoordinates.reshape((3,numberofvertices))
    #to store the deformation on the cells
    faceids = []
    deformations = []
    ####iterating the vertices to feed the new coordinates in
    vertices = qd.CellVertexIterator(cell)
    vertex = vertices.next()
    counter = 0#counter to change vertex positions
    while vertex != None:
        vertex.setXcoordinate(tempcoordinates[0,counter])#setting x coordinate
        vertex.setYcoordinate(tempcoordinates[1,counter])#setting y coordinate
        vertex.setZcoordinate(tempcoordinates[2,counter])#setting z coordinate
        counter += 1
        vertex = vertices.next()
    #####################################################
    cell = settingParameters(cell)
    #####################################################
    #calculating the deformation of Area 
    faces = qd.CellFaceIterator(cell)
    face1 = faces.next()
    while face1 != None:
        if face1.getID() == 1: 
            face1 = faces.next()
            continue
        #print face1.getID()
        targetarea = face1.getTargetArea()
        currentarea = face1.getAreaOfFace()
        #print "targetarea :", targetarea, "current area: ",currentarea, "difference :",targetarea-currentarea
        #change of area (strain on area %)
        faceids.append(face1.getID())
        deformations.append(100*(currentarea - targetarea)/targetarea)
        face1 = faces.next()
    ######################################################
    #returning the total energy
    first= cell.getFirstTerm()
    second = cell.getSecondTerm()
    third  = cell.getThirdTerm()
    fourth = cell.getFourthTerm()
    volume = cell.getVolume()
    #release the current cell
    #printing the counter and if the passed value is different than the previous values passed
    #print maincounter, energyvalue
    #print np.subtract(tempcoordinates, tempcoordstore)
    return [first,second, third,volume, [faceids, deformations],fourth]

#####################################################
############################################################################
############################################################################
###                    Current Mean Determinant Target Area Plot         ###
############################################################################
def plotMeanTargetArea(cell,meandeterminantarray,timearray):
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    sumtargetdeterminant = 0.
    numofface = 0.
    while face != None:
        if face.getID() == 1:
            face = faces.next()
            continue
        sumtargetdeterminant += face.getTargetFormMatrixDeterminant()
        numofface += 1.
        face = faces.next()
    meandeterminant = sumtargetdeterminant/numofface
    meandeterminantarray.append(meandeterminant)
    ###now plotting
    plt.figure(10)
    plt.title("Mean determinant of Target Form Matrix ")
    plt.ylabel("Mean determinant")
    plt.xlabel("time")
    plt.plot(timearray, meandeterminantarray,'-x')
    plt.savefig('mean_target_form_determinant.png', transparent = True)
    plt.clf()
    plt.close()
    ###############
    return meandeterminantarray
"""
#########################################################################################################
#Correcting the outer layer initial bending
#########################################################################################################
faces = qd.CellFaceIterator(cell)
face = faces.next()
while face != None:
    if face.getID()<=(perimeterHexNum+1):
        if face.getID() == 1:
            face = faces.next()
            continue
        edges = qd.FaceEdgeIterator(face)
        edge = edges.next()
        while edge != None:
            vert = edge.Dest()
            x = vert.getNonCentralisedProjectedXcoordinate(face.getID())
            y = vert.getNonCentralisedProjectedYcoordinate(face.getID())
            z = vert.getNonCentralisedProjectedZcoordinate(face.getID())
            vert.setXcoordinate(x)
            vert.setYcoordinate(y)
            vert.setZcoordinate(z)
            #print x, y, z
            edge = edges.next()
    face = faces.next()
"""
########################################################################################################################
##intial parameters for optimization, the coordiantes of vertices
vertices = qd.CellVertexIterator(cell)
vertex = vertices.next()
vertexid = np.array([])
xcoord = np.array([])
ycoord = np.array([])
zcoord = np.array([])
while vertex != None: 
    #saving the ids
    vertexid = np.append(vertexid, vertex.getID())
    xcoord = np.append(xcoord, vertex.getXcoordinate())
    ycoord = np.append(ycoord, vertex.getYcoordinate())
    zcoord = np.append(zcoord, vertex.getZcoordinate())
    vertex = vertices.next()
coordinates = np.concatenate((xcoord,ycoord,zcoord))
numberofvertices = len(xcoord)
########################################################################################################################
##################################################################
###                bounds for optimization                    ####
##################################################################
##################################################################
###                BOUNDS for optimization                    ####
##################################################################
def getBoundsForOptimization(cell,radius, perimeterVertexNum):
    #getting the coordinates of the cell
    ##intial parameters for optimization, the coordiantes of vertices
    vertices = qd.CellVertexIterator(cell)
    vertex = vertices.next()
    #vertexid = np.array([])
    xcoord = np.array([])
    ycoord = np.array([])
    zcoord = np.array([])
    while vertex != None: 
        #saving the ids
        #vertexid = np.append(vertexid, vertex.getID())
        xcoord = np.append(xcoord, vertex.getXcoordinate())
        ycoord = np.append(ycoord, vertex.getYcoordinate())
        zcoord = np.append(zcoord, vertex.getZcoordinate())
        vertex = vertices.next()
    coordinates = np.concatenate((xcoord,ycoord,zcoord))
    numberofvertices = cell.countVertices()
    #print "lenght of xcoord : ", len(xcoord), "totalvertex", totalvertex
    basefreedom = 2*(radius+1)#in this case: we just let the vertices in the boundary move freely in their own X-Y plane
    epsbasefreedom = (1./5)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
    print "small freedome for base vertices : 20%"
    print "basefreedom", basefreedom
    print "epsbasefreedom", epsbasefreedom
    fac = 6
    #boundary set so that outer (periphery) vertices do not move
    uboundsx = np.zeros((numberofvertices))
    uboundsy = np.zeros((numberofvertices))
    uboundsz = np.zeros((numberofvertices))
    lboundsx = np.zeros((numberofvertices))
    lboundsy = np.zeros((numberofvertices))
    lboundsz = np.zeros((numberofvertices))
    ###############################
    vertices = qd.CellVertexIterator(cell)
    vertex = vertices.next()
    vertcounter = 0
    while vertex != None: 
        currentvertid = vertex.getID()
        #print vertcounter, " -> ",currentvertid
        ## DOME VERTEX
        if currentvertid > perimeterVertexNum:
            uboundsx[vertcounter] = xcoord[vertcounter]+0.5
            uboundsy[vertcounter] = ycoord[vertcounter]+0.5
            uboundsz[vertcounter] = zcoord[vertcounter]+0.5
            lboundsx[vertcounter] = xcoord[vertcounter]-0.5
            lboundsy[vertcounter] = ycoord[vertcounter]-0.5
            lboundsz[vertcounter] = zcoord[vertcounter]-0.5
        else:## CYLINDER VERTEX
            uboundsx[vertcounter] = xcoord[vertcounter]
            uboundsy[vertcounter] = ycoord[vertcounter]
            uboundsz[vertcounter] = zcoord[vertcounter]
            lboundsx[vertcounter] = xcoord[vertcounter]
            lboundsy[vertcounter] = ycoord[vertcounter]
            lboundsz[vertcounter] = zcoord[vertcounter]
        vertcounter += 1
        vertex = vertices.next()
    ### Upper bounds ###
    upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
    ##Lower bounds
    lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
    return (upperbounds, lowerbounds)
########################################################################################################################
#######################################################################################
#starting the simulation
print "Sequence started: "
for i in range(10,0,-1):
    print i, "..."
##############################################################################################################################################################################
############################################################################
###                    NLOPT configuration                               ###
############################################################################
#optimizing object
#SBPLX : Constraint optmization by Linear Approximations
# Local derivative-free optimization
global maincounter
#global tempcoordstore#array to store calculation parameter
maincounter = 0
stepcounter = 0
tempcoordstore = np.copy(coordinates)
###opt optimizer
#optimizer for nlopt
opt = nlopt.opt(nlopt.LN_SBPLX,numberofvertices*3)
###########################################################################
##########         Conditions on the relaxation                  ##########
###########################################################################
upperbounds, lowerbounds = getBoundsForOptimization(cell,radius,perimeterVertexNum)
opt.set_lower_bounds(lowerbounds)
opt.set_upper_bounds(upperbounds)
#setting minimization condition
opt.set_min_objective(energyobjective)
#setting the initial step for intiial guess
opt.set_initial_step(initialstep)
##stop criteria  of relative tolerance value
#opt.set_ftol_rel(tolerance)
##stop criteria of absolute tolerance value
#opt.set_ftol_abs(tolerance/100.)#two order of magnitude less than the relative
#setting tolerance for relative change in xvalues
opt.set_xtol_abs(tolerance)
print "xtol_abs : : ", tolerance
print "initial step ::", initialstep
###########################################################################
#opt.set_maxtime(interval)#time interval for optimization time limit
############################################################################
#running the code for number of times : totaltime/interval
#making an array to regularly store the optimized parameters
xopt = np.copy(coordinates)
#saving the initial cell
latdev.plot3DCell(cell,name = r'time=%(time)d.png'%{"time":stepcounter})
latdev.plotSurface(cell, numOfLayer,name = r'surface_time={0:d}.png'.format(stepcounter))
##########################################################################################
#Initialising the counters
##########################################################################################
#number of relaxation time steps
#Growth counter
growthtimearray = np.array([0.])
#array for which stepcounter did the growth took place
stepcounterGrowthArray = []
##########################################################################################
initialenergy = energyobjective(coordinates,0.)
tempterms = termsofenergy(coordinates)
#array to store the energy and time 
energyarray = np.array([initialenergy])
firsttermarray = np.array([tempterms[0]])
secondtermarray = np.array([tempterms[1]])
thirdtermarray = np.array([tempterms[2]])
timearray = np.array([stepcounter])
volumearray = np.array([tempterms[3]])
deformations = tempterms[4]
meanDeformation = np.array([np.mean(np.array(deformations[1]))])
stdDeformation = np.array([np.std(np.array(deformations[1]))])
fourthtermarray = np.array([tempterms[5]])
meandeterminantarray = []
#plotting deformation plot
plt.figure(11)
plt.title("Area strain Percentage")
plt.ylabel("Areastrain%")
plt.xlabel("faceids")
plt.plot(deformations[0], deformations[1], 'x')
plt.savefig('deformation_plot_0.png', transparent = True)
plt.clf()
#
saveGrowthFigures(cell,numOfLayer, 0)

#############################################################################################
####Plotting the CurrentFormMatrix, TargetFormMatrix and StrainMatrix
##getting matrices
##############################################################################################
def plotMatrices(cell, targetfaceid,matrixarray,growthtimearray, stepgrowtharray):
    #getting the matrix array
    mudeterminant = matrixarray[0]
    currentFormdeterminant = matrixarray[1]
    targetFormdeterminant = matrixarray[2]
    #calculating the matrices
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face!=None:
        faceid = face.getID()
        if faceid != targetfaceid:
            face= faces.next()
            continue
        mudeterminant.append((face.getMu1()*face.getMu4()-face.getMu2()*face.getMu3()))
        currentFormdeterminant.append((qd.getCurrentFormMatrix(face,0,0)*(qd.getCurrentFormMatrix(face,1,1))-qd.getCurrentFormMatrix(face,1,0)*(qd.getCurrentFormMatrix(face,0,1))))
        targetFormdeterminant.append((qd.getTargetFormMatrix(face,0,0)*(qd.getTargetFormMatrix(face,1,1))-qd.getTargetFormMatrix(face,1,0)*(qd.getTargetFormMatrix(face,0,1))))
        break
        #face = faces.next()
    #starting the plot
    plt.figure(21)
    plt.title("Determinant of Matrices for Face : %d"%targetfaceid)
    plt.plot(growthtimearray, mudeterminant,'-s',markeredgewidth=0.0,markeredgecolor=None,color='b', label="Mu")
    plt.plot(growthtimearray, currentFormdeterminant,'-o',markeredgewidth=0.0,markeredgecolor=None,color='g', label="CF")
    plt.plot(growthtimearray, targetFormdeterminant,'-v',markeredgewidth=0.0,markeredgecolor=None, color='r',label="TF")
    #making verticle lines for growthsteps
    for vline in stepgrowtharray:
        plt.axvline(x=vline, c = 'r')
    plt.legend(loc='best')
    plt.savefig('matrix_plot_faceid_%d.png'%targetfaceid,transparent=True)
    plt.clf()
    plt.close('all')
    ######
    return [mudeterminant,currentFormdeterminant,targetFormdeterminant]
###########################################################################################
###########################################################################################
#Faces ID : to plot Matrices 
faceidarray = []#array to store faceids
matrixdictionary = {}#dictionary storing the matrices determinant [[mu],[CF],[TF]]
faces =qd.CellFaceIterator(cell)
face = faces.next()
while face != None:
    if face.getID()==1:
        face = faces.next()
        continue
    if face.getID()%10 == 0:
        faceidarray.append(face.getID())
        matrixdictionary[face.getID()] = [[],[],[]]
    face = faces.next()
print "Matrix Dictionary \n", matrixdictionary
###########################################################################################
#now plotting for all the faces
for fid in faceidarray:
    matrixdictionary[fid] = plotMatrices(cell,fid,matrixdictionary[fid],growthtimearray,stepcounterGrowthArray)
###########################################################################################
#plotting the determinant of Target Form Matrix
meandeterminantarray = plotMeanTargetArea(cell,meandeterminantarray,growthtimearray)
###########################################################################################
###########################################################################################
#start of simulation ::
for growthcounter in range(growthSteps):
    np.save(r'coordinates_after_step={0:d}'.format(int(stepcounter)),xopt)#saving initial coordinates
    currentcoordinates = np.array([])
    if np.isnan(initialenergy):
        print "energy NaN ! :", energyarray[-1]
        printmatrixenergy(cell)
    else:
        print "start -----> ", "initial energy : ", energyarray[-1]
        ###################################################################################################
        ###                            PERFORMING GROWTH                                               ####
        ###################################################################################################
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "           Starting growth       STEP : : ", growthcounter+1
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        cell = growFaces(cell)
        ###############################
        #adding the current stepcounter to stepcounterGrowthArray
        stepcounterGrowthArray.append(stepcounter)
        ###############################
        saveGrowthFigures(cell,numOfLayer,growthcounter+1)
        growthtimearray = np.append(growthtimearray, growthcounter+1)
        np.save(r'growth_images/growth_coordinates_before_step={0:d}'.format(int(growthcounter+1)),xopt)#saving initial coordinates
        meandeterminantarray = plotMeanTargetArea(cell, meandeterminantarray, growthtimearray)
        ##################################################################################################
        ##################################################################################################
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "Growth Step : ", growthcounter+1, " Completed ! Next step to relaxation again!" 
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ###################################################################################################
        ###                            START OF THE RELAXATION                                         ####
        ###################################################################################################
        """
        print "###########################################################################################"
        print "Printing Matrix of Cells "
        print "###########################################################################################"
        printmatrixenergy(cell)
        """
        print "##############################################################################################################################################"
        print "now Optimization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! step  : ", stepcounter+1
        print "##############################################################################################################################################"
        #print "Intervalspassed : ", intervalpassed, "time passed : ", (-1+intervalpassed)*interval
        #printareaofface()
        #printmatrixenergy(cell)
        #############################################################################
        ##optimizing the returning the optimized parameters
        #############################################################################
        #comparing the currentarray with the xopt array to see if last output is the one being used for this optimization
        #print "Is the last output being used for optimization ? : ", np.array_equal(currentcoordinates,xopt)
        #print xopt
        #print nlopt.nlopt_get_initial_step(opt, xopt)
        #initialstep = opt.get_initial_step(xopt)
        #print initialstep
        #np.savetxt('coordinates_before_step_%(step).d'%{"step":intervalpassed},xopt)
        relaxationstarttime = time.time()
        upperbounds, lowerbounds = getBoundsForOptimization(cell,radius,perimeterVertexNum)
        opt.set_lower_bounds(lowerbounds)
        opt.set_upper_bounds(upperbounds)
        print "relaxation start time  :: ", str(datetime.now())
        try:
            xopt = opt.optimize(xopt)
            returncode = opt.last_optimize_result()#the result from the optimization 
            stepcounter += 1
        except nlopt.RoundoffLimited:
            #printmatrixenergy(cell)
            print "round off error"
            optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
            returncode = opt.last_optimize_result()#the result from the optimization 
            print "last optimized energy : ", optimizedenergy
            print "return code : ", returncode
            stepcounter += 1
            endofprogram(xopt, energyarray, firsttermarray, secondtermarray, thirdtermarray, volumearray, error = True)
            break
        except ValueError:
            print "Value error : ", xopt
            optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
            returncode = opt.last_optimize_result()#the result from the optimization 
            print "last optimized energy : ", optimizedenergy
            print "return code : ", returncode
            stepcounter += 1
            endofprogram(xopt, energyarray, firsttermarray, secondtermarray, thirdtermarray, volumearray, error = True)
            #quit()
        relaxationendtime = time.time()
        print "Time taken for relaxation step ", stepcounter, ": ", relaxationendtime - relaxationstarttime, "   in minutes : ", (relaxationendtime - relaxationstarttime)/60., "   in Hours : ", (relaxationendtime - relaxationstarttime)/3600.
        #getting initial step size
        np.save(r'coordinates_after_step={0:d}'.format(int(stepcounter)),xopt)
        currentcoordinates = np.copy(xopt)
        ############################################################################
        # Taking output from optimizer and plotting                                #
        ############################################################################
        optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
        optimizationTermination(returncode)
        #currenttime = interval*(1+intervalpassed)#calculating current time
        #printing the values
        print "Current Step : ", stepcounter, "returncode : ", returncode
        print "Optimized energyvalue : ", optimizedenergy
        #adding to the energy array
        energyarray  = np.append(energyarray, optimizedenergy)
        timearray = np.append(timearray, stepcounter)
        #grabbing the terms of energy
        tempterms = termsofenergy(xopt)
        #appending to the terms array
        firsttermarray = np.append(firsttermarray, tempterms[0])
        secondtermarray = np.append(secondtermarray, tempterms[1])
        thirdtermarray = np.append(thirdtermarray, tempterms[2])
        volumearray = np.append(volumearray, tempterms[3])
        deformations = tempterms[4]#taking the deformation datas
        meanDeformation = np.append(meanDeformation, np.mean(np.array(deformations[1])))
        stdDeformation = np.append(stdDeformation, np.std(np.array(deformations[1])))
        fourthtermarray = np.append(fourthtermarray,tempterms[5])
        #rearraging the optimizer 1D output to 3xN array
        newvertexpositions = xopt.reshape((3,numberofvertices))
        #setting the vertices of cell with new coordinates
        vertices = qd.CellVertexIterator(cell)
        vertex = vertices.next()
        counter = 0#counter to change vertex positions
        while vertex != None:
            vertex.setXcoordinate(newvertexpositions[0,counter])#setting x coordinate
            vertex.setYcoordinate(newvertexpositions[1,counter])#setting y coordinate
            vertex.setZcoordinate(newvertexpositions[2,counter])#setting z coordinate
            counter += 1
            vertex = vertices.next()
        ###Plotting all the plots
        #saving the plot
        latdev.plot3DCell(cell,name=r'time=%(time)d.png'%{"time":stepcounter})
        #plotting the surface only of cell
        latdev.plotSurface(cell, numOfLayer,name = r'surface_time={0:d}.png'.format(int(stepcounter)))
        #plotting the energy plot 
        plt.figure(2)
        plt.title("total energy")
        plt.ylabel("optimized energy values")
        plt.xlabel("time")
        plt.plot(timearray, energyarray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('energy_plot.png', transparent = True)
        plt.clf()
        #plotting the first term energy
        plt.figure(3)
        plt.title("First term of Energy")
        plt.ylabel("optimized First Term of energy")
        plt.xlabel("time")
        plt.plot(timearray, firsttermarray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('firstterm_plot.png', transparent = True)
        plt.clf()
        #plotting the first term energy
        plt.figure(4)
        plt.title("Second term of Energy")
        plt.ylabel("optimized Second Term of energy")
        plt.xlabel("time")
        plt.plot(timearray, secondtermarray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('secondterm_plot.png', transparent = True)
        plt.clf()
        #plotting the first term energy
        plt.figure(5)
        plt.title("Third term of Energy")
        plt.ylabel("optimized Third Term of energy")
        plt.xlabel("time")
        plt.plot(timearray, thirdtermarray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('thirdterm_plot.png', transparent = True)
        plt.clf()
        #plotting the first term energy
        plt.figure(6)
        plt.title("Volume of Cell")
        plt.ylabel("Volume")
        plt.xlabel("time")
        plt.plot(timearray, volumearray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('volume_plot.png', transparent = True)
        plt.clf()
        #plotting the energy plot 
        plt.figure(7)
        plt.title("Area Percent strain")
        plt.ylabel("Areastrain%")
        plt.xlabel("faceids")
        plt.plot(deformations[0], deformations[1], 'x')
        plt.savefig('deformation_plot_{0:d}.png'.format(int(stepcounter)), transparent = True)
        plt.clf()
        #plotting mean deformation
        plt.figure(8)
        plt.title("Mean relative deformation ")
        plt.ylabel("relative deformation")
        plt.xlabel("time")
        plt.errorbar(timearray, meanDeformation, stdDeformation, marker='x')
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('mean_deformation_plot.png', transparent = True)
        plt.clf()
        #plotting z-projection summation - Fourth term
        plt.figure(9)
        plt.title("Fourth term : Z-projection summation ")
        plt.ylabel("zprojection^2")
        plt.xlabel("time")
        plt.plot(timearray, fourthtermarray)
        for vline in stepcounterGrowthArray:
            plt.axvline(x=vline, c = 'r')
        plt.savefig('fourthterm_plot.png', transparent = True)
        plt.clf()
        #closing all figures
        plt.close('all')
        ################################################################################################
        #now plotting determinant of matrices for all the faces
        for fid in faceidarray:
            matrixdictionary[fid] = plotMatrices(cell,fid,matrixdictionary[fid],timearray,stepcounterGrowthArray)
print "done"
#######function to print energy and Mo and Mc
endofprogram(xopt, energyarray, firsttermarray, secondtermarray, thirdtermarray, volumearray, fourthtermarray)#saving the last relaxed coordinates
#printmatrixenergy(cell)
###################################################################################################################################

