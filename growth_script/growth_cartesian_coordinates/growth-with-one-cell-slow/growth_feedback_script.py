####################################################################################################################################
#                       CLUSTER CODE 
#           Description : 
#           ** The code slows one cell to grow slower than the rest of the cells 
#           ** 
#                 This code uses Cartesian Coordinate
#                 This script is made to be general propose for simulating Quadedge Tissue with options to run either
#                         a. Cylindrical + Dome top initial condition
#                         b. Only Dome top intial condition
#            Task :
#         Growth of Surface , can be
#							a. Cylinder 
#							b. Hemi-sphere
#					with Radial degree of freedom
#			##############################################
#            Constraints (BOUNDS) : 
#					a. WITH Radial CONSTRAINT on Cylindrical vertex in case of Cylinder
#					b. The external vertices are bounded to their position
#                   c. Infinite Penalty for crossng bending threshold of initial fourth term value
#			##############################################
#            Initial condition :
#                   a. Hemi-sphere is uniform Dome
#					b. Cylinder is made of "Uniform" hexagons
#			##############################################
#            Growth : 
#                   a. Growth is present
#					b. Is definined by the function in use
#							i. InflatedGrowth() -> growth is random with gaussian distribution proportional to current shape
#							ii. randomGrowth()  -> uniformly distributed random growth proportional to current shape 
#			##############################################
#            Termination Conditon : 
#                    a. SImulation runs until given number of growth step is completed
#            ##############################################
#			 Algorithm :
#                    a. SBPLX Algorithm is used
####################################################################################################################################
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
##################################################################################################################################################
programstart = time.time()
print "This is a test for Growth !!"
algorithm = "SBPLX"
print "The Used Algorithm :: ", algorithm
################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#adding arguments
parser.add_argument("-c","--cylinder", help = "if option is used, the initial condition is Cylinder, else by default it is Dome", action= "store_true")
parser.add_argument("-x","--exp", help = "If used, the growth is  Exponential, if not Growth is strain driven.", action= "store_true")
parser.add_argument("-e","--error",help = "value for relative error tolerance, default = 1e-5",
                                   default = 10**(-5), type = float)
parser.add_argument("-a","--alpha",help = "value to set for Alpha (weigth of first term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-n","--eta",help = "value to set for eta,  (weigth of Feedback in growth), default = 0",
                                   default = 0., type = float)
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
                                   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
                                   default = 0., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
                                   default = 0.001, type = float)
parser.add_argument("-k","--kappa",help = "value to set for kappa growth Rate of Mo, default = 1.0",
                                   default = 1., type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 180",
                                   default = 180., type = float)
parser.add_argument("-t","--time",help = "total steps for relaxation to run, default = 3000",
                                   default = 1000, type = int)
parser.add_argument("-i","--interval",help = "time interval for saving the plot of tissue (other datas if applicable)",
                                   default = 30, type = int)
parser.add_argument("-l","--layer",help = "number of layers of cells to simulate",
                                   default = 5, type = int)
parser.add_argument("-s","--cylinderalpha",help="alpha for cylinder to be different from dome", default = 0.)
parser.add_argument("-g","--gamma", help = "Gamme is the pressure from underneath the epidermis, comming from lower level cells. acting as the volume maximizing agent", default = .1, type = float)
parser.add_argument("-v","--varyfeedback",help="vary the Feedback as multiple of 0.00001*(10)**(passedvalue)",default = 0, type = int)
parser.add_argument("-j","--jobname",help="jobname if given, the new folder to be created will be named jobname", default = None, type = str)
#parser.add_argument("-n","--numberofimage",help = "number of images or data points to store during simulation", default = 10, type = int)
############################################################################
###         Function to Slow growth of a cell and its neighbourhood      ###
############################################################################
def slowGrowth(cell,faceid):
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while (face != None):
        if face.getID() == faceid:
            print "Slowed Growth ID :", faceid
            #print face.getKappa()
            face.setKappa(1.5)
            face.setGrowthVar(0.5)
            edges = qd.FaceEdgeIterator(face)
            edge = edges.next()
            while edge != None:
                rightFace = edge.Right()
                rightFace.setKappa(1.5)
                rightFace.setGrowthVar(0.5)
                edge = edges.next()
        face = faces.next()
    return
###################################################################################################
#								Setting Parameters
##################################################################################################
args = parser.parse_args()#arguments passed
#gathering the values passed by argument
alpha = args.alpha
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
totaltime = args.time
interval = args.interval
numOfLayer = args.layer
gamma = args.gamma
kappa = args.kappa
eta = args.eta
anglethreshold = args.angle
tolerance = args.error
initialstep = 10.**(-2)
growthSteps = 500 #performing 100 growth steps
initialStrain = 0.
jobname = args.jobname
### FACE TO SLOW THE GROWTH ###
if numOfLayer == 4:
    slowgrowthID = 31
elif numOfLayer == 6: 
    slowgrowthID = 57
elif numOfLayer == 8:
    slowgrowthID = 132
else: 
    slowgrowthID = 20
##################################################################################################
#						chaning pressure if the Vary pressure is used
#       varyfeedback is expected to be from the following set
#       varyfeedback = {1,2,3,4,5,...} 
#       for varyfeedback = 1 -> eta = 0, hence, No Feedback
#       for varyfeedback = {2,3,4,5,...} -> eta = 10^(-5)*10^(varyfeedback)
##################################################################################################
#########################
"""
#Supply args.varyfeedback in range (0, 25)#
# Division by 5 : -> Quotient is pressure & Remainder is feedback strength
varyVariable = args.varyfeedback - 1
varypressure = varyVariable // 5 
varyfeedback = varyVariable % 5 
gamma_array = np.insert(np.linspace(0.001,0.1,4),0,0)
#gamma = (10.**(-3.))*10.**(varypressure)
gamma = gamma_array[varypressure] #picking one gamma from the above array - gamma varies [0.005,0.1]
if varyfeedback == 0:
    eta = 0.# for varyfeedback = 1 -> 0 Feedback
else:
    eta = (10.**(-3.))*10.**(varyfeedback) #varrying feedback here for 4 orders of magnitude
"""
#Supply args.varyfeedback in range (0, 25)#
# Division by 5 : -> Quotient is pressure & Remainder is feedback strength
varyVariable = args.varyfeedback - 1
varypressure = varyVariable // 5
#varyfeedback = varyVariable % 5 
varyfeedback = 0
gamma_array = np.linspace(0.0001,0.1,25)
#gamma = (10.**(-3.))*10.**(varypressure)
gamma = gamma_array[varyVariable] #picking one gamma from the above array - gamma varies [0.005,0.1]
if varyfeedback == 0:
    eta = 0.# for varyfeedback = 1 -> 0 Feedback
else:
    eta = (10.**(-3.))*10.**(varyfeedback) #varrying feedback here for 4 orders of magnitude
##################################################################################################
#						Creating the output folder
##################################################################################################
if not os.path.exists("output"):
     os.makedirs("output")
os.chdir("output")
#now making path for simulation result
newpath = r'a=%(alpha).4f_g=%(gamma).4f_n=%(eta).4f_is=%(strain).2f_layer=%(layer).1d'%{"strain":initialStrain,"alpha":alpha,"beta":beta,"gamma":gamma,"zeta":zeta,"eta":eta,"time":totaltime,"layer":numOfLayer}
#newpath = r'M0=identity_a=%(alpha)f_b=%(beta)f_p=%(pressure)f_t=%(time)f_layer=%(layer)d'%{"alpha":alpha,"beta":beta,"pressure":pressure,"time":totaltime,"layer":numOfLayer}
if not os.path.exists(newpath):
    os.makedirs(newpath)
#change directory to newpath
os.chdir(newpath)
#making a directory to save coordinates after each relaxation+growth step before growth kicks in 
growthdirectory = 'growth_images'
if not os.path.exists(growthdirectory):
    os.makedirs(growthdirectory)
##########################################################################################################
#                                SETTING INITIAL CONDITION 
##########################################################################################################
#######  Geometrical parameters of the tissue  #########
Length = 1.# length of the sides of the faces
radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
numofface = 3*numOfLayer*(numOfLayer-1)+1
perimeterVertexNum = 6+12*(numOfLayer-1) #number of vertices in the perimeter of lattice
perimeterHexNum = (numOfLayer!=1)*(6*(numOfLayer-1)-1) + 1 # number of perimeter hexgons
#####################################################################################################################
#####                 MAKING THE INITIAL QUADEDGE TISSUE
#####################################################################################################################
#    Function to make the tissue
#		a. makeCenteredHexagonLattice(numOfLayer, printcondition = False, projection = True) 
#	 Projection condition : it is to have the external vertex projected on to an external circle
##############################################################################################
cell, newedge = latdev.makeCenteredHexagonLattice(numOfLayer, True, False) 
####################################################################################
#put down the paramters to cell
####################################################################################
cell.setAlpha(alpha)
cell.setBeta(beta)
cell.setPressure(pressure)
cell.setZeta(zeta)
cell.setGamma(gamma)
cell.setKappa(kappa)
cell.setEta(eta)
cell.setConvexAngleThreshold(anglethreshold)
cell.setInitialStrain(initialStrain)
cell.setGrowthVar(kappa/2.)
if args.cylinderalpha != 0.:#in case this option is used
    cell = setCylinderAlpha(cell,args.cylinderalpha)
    print "Alpha for Cylinder : ", args.cylinderalpha
####################################################################################
# Making Initial Condition 
#			a. args.cylinder == False -> Make Dome 
#			b. args.cylinder == True -> Make Cylinder
####################################################################################
if args.cylinder: 
	cell = makeCylinder(cell,numOfLayer)
	cellradius = radius/1.5
else:
	cell = makeDome(cell,numOfLayer)
	cellradius = radius
############################################################################################
# Slowing the growth of a target face 
############################################################################################
#slowGrowth(cell,slowgrowthID) ##Face id to slow the face growth
vmax=0 #change vmax if exceptionally high kappa for some cells are used
#print "All slow have small << 1 growth except neighbourhood of Face : ", slowgrowthID
print "UNIFORM GROWTH OF ALL CELLS"
####################################################################################
# Setting the initial parameters for the tissue now
####################################################################################
#cell = settingFirstParameters(cell)#setting all the parameters of the cell
#cell.setCylindrical()#SETTING Cylindrical coordinates for all the vertices
cell.setInitialParameters()
numberofvertices = cell.countVertices()
bendingThreshold = cell.getBendingThreshold()
#########    Printing Parameters    ########################
printCellParameters(cell,numOfLayer)
print "---------------------------------------------------"
print " CelL Growth : Gaussian Width : ", cell.getGaussianWidth()
print "---------------------------------------------------"
print "tolerance : ", tolerance
print "gamma :",  gamma
print "Cell Kappa : ", cell.getKappa()
print "Cell Growth Var ", cell.getGrowthVar()
print "The Bending Energy Bound : ", bendingThreshold
print "Angle Threshold : ", cell.getConvexAngleThreshold()
print "Initial Strain on Faces : ",cell.getInitialStrain()
print "------------------------------------------------"
############################################################################################################
#####################################################
#		Saving Initial Figures
#####################################################
latdev.plot3DCell(cell,name = r'initial_dome.png')
latdev.plotSurface(cell, numOfLayer,name = r'initial_dome_surface.png')
#######################################################################################
####        IMPORTANT FUCNTIONS FOR SIMULATION                                     ####
#######################################################################################
############################################################################
###                 objective function for NLOPT                         ###
############################################################################
def cartesianEnergyObjective(inputcoordinates,grad):
    global numOfLayer, alpha, beta, pressure, cell, functionCallCounter,bendingThreshold
    functionCallCounter += 1
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
    ######################################################
    #cell = settingParameters(cell)
    cell.setParameters()
    ######################################################
    #Checking the Bans : fourtherm < bending energy & Cell is Convex ! 
    """
    fourthterm = cell.getFourthTerm()
    bendingThreshold = cell.getBendingThreshold()
    if (fourthterm > bendingThreshold):
      return float('inf')
    if not (cell.isConvex()): #if cell is not convex
      return float('inf')
    """
    #returning the total energy
    energyvalue = cell.getEnergyCartesianVolume()
    return energyvalue

############################################################################
############################################################################
###                    NLOPT configuration                               
###			-------------------------------------------------
###			Optimizing object : cylinderEnergyObjective
### 		Algorithm : SBPLX
###						Constraint optmization by Linear Approximations
###						Local derivative-free optimization
############################################################################
coordinates = getCartesianCoordinate(cell)
tempcoordstore = np.copy(coordinates)
####################################################
#opt optimizer
opt = nlopt.opt(nlopt.LN_SBPLX,numberofvertices*3)#
###########################################################################
##########         Conditions on the relaxation                  ##########
###########################################################################
#setting minimization condition
opt.set_min_objective(cartesianEnergyObjective)#
#setting the initial step for intiial guess
opt.set_initial_step(initialstep)#
opt.set_xtol_abs(tolerance)#tolerance for convergence
print " Dimension of the Problem (NLOPT) : ", numberofvertices*3
print "xtol_abs : : ", tolerance
print "initial step ::", initialstep
########################################################################################
#starting the final intialization of variables needed for optimization & saving figures
########################################################################################
functionCallCounter = 0
########################################################################################
xopt = np.copy(coordinates)
#saving the initial cell
latdev.plot3DCell(cell,name = r'time=%(time)d.png'%{"time":0})
latdev.plotSurface(cell, numOfLayer,name = 'surface_time=%03.d.png'%(0))
########################################################################################
initialenergy = cartesianEnergyObjective(coordinates,0.)
tempterms = cartesianTermsofenergy(cell,coordinates)
#array to store the energy and time 
energyarray = np.array([initialenergy])
firsttermarray = np.array([tempterms[0]])
secondtermarray = np.array([tempterms[1]])
thirdtermarray = np.array([tempterms[2]])
timearray = np.array([0])
volumearray = np.array([tempterms[3]])
deformations = tempterms[4]
meanDeformation = np.array([np.mean(np.array(deformations[1]))])
stdDeformation = np.array([np.std(np.array(deformations[1]))])
fourthtermarray = np.array([tempterms[5]])
meandeterminantarray = []
relaxation_time_length=[]
functionCallCounterArray = []
#plotting deformation plot
plt.figure(11)
plt.title("Area strain Percentage")
plt.ylabel("Areastrain%")
plt.xlabel("faceids")
plt.plot(deformations[0], deformations[1], 'x')
plt.savefig('deformation_plot_0.png', transparent = True)
plt.clf()
###########################################################################################
###		Getting matrixdictionary : storage for Mu-,CF- and TF-Matrix
###########################################################################################
matrixdictionary, faceidarray = getInitialMatrixDictionary(cell)
###########################################################################################
#Adding to MatrixDictionary and Plotting for all the faces : Both are done through this 
#function
###########################################################################################
matrixdictionary = plotMatrixDictionary(cell,faceidarray, matrixdictionary)
###############################################################################################
#plotting the determinant of Target Form Matrix
meandeterminantarray = plotMeanTargetArea(cell,meandeterminantarray)#initial growth step = 0
###############################################################################################
## Saving Coordinates and Target Form Matrix
np.save(r'coordinates_step={0:d}'.format(int(0)),xopt)
#saving the TargetFormMatrix of Cell
saveTargetFormMatrix(cell,0)
##########################
plotGrowthRateSurface(cell,numOfLayer,name='growthRateSurface=%03.d.png'%(0),vmax=vmax)
###################################################################################################################
###							 			START OF SIMULATION 									
###################################################################################################################
for growthcounter in range(growthSteps):
    ####################################################################################
    ##	Check if Energy is NAN or not
    ##		If NAN : QUIT SIMULATION
    ##		IF NOT NAN : Continue
    ####################################################################################
    if np.isnan(initialenergy):
        print "energy NaN ! :", energyarray[-1]
        printmatrixenergy(cell)
        break
    else:
        print "Start Step :", growthcounter ," Initial energy : ", energyarray[-1]
        print "##############################################################################################################################################"
        print "now Optimization  : Step  : ", growthcounter
        print "##############################################################################################################################################"
        relaxationstarttime = time.time() #time before start of optimization
        upperbounds, lowerbounds = getBoundsForOptimization(cell,cellradius,perimeterVertexNum)
        #upperbounds, lowerbounds = getStrongBottomBoundsForOptimization(cell,cellradius,perimeterVertexNum)
        print " Upper Bounds greater than Lower Bounds : ", np.all(np.greater_equal(upperbounds,lowerbounds))
        opt.set_lower_bounds(lowerbounds)
        opt.set_upper_bounds(upperbounds)
        print "relaxation start time  :: ", str(datetime.now())
        functionCallCounter = 0
        #########################################################################################
        #		Start of relaxation
        #########################################################################################
        try:
            xopt = opt.optimize(xopt)
            returncode = opt.last_optimize_result()#the result from the optimization 
            #stepcounter += 1
        except nlopt.RoundoffLimited:
            #printmatrixenergy(cell)
            print "round off error"
            optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
            returncode = opt.last_optimize_result()#the result from the optimization 
            print "Termination Condition : ", optimizationTermination(returncode)
            print "last optimized energy : ", optimizedenergy
            print "return code : ", returncode
            #stepcounter += 1
            endofprogram(xopt, energyarray, firsttermarray, secondtermarray, thirdtermarray, volumearray, error = True)
            break
        except ValueError as caughterror:
            print "Value error : ", xopt
            print "Is any value Nan : ", np.all(np.isnan(xopt))
            print "Caught Error : ", caughterror
            print " IS Cell Convex ? : ", cell.isConvex()
            print " Fourth Term Value : ", cell.getFourthTerm()
            print " Bending Threshold : ", cell.getBendingThreshold()
            optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
            returncode = opt.last_optimize_result()#the result from the optimization 
            print "Termination Condition : ", optimizationTermination(returncode)
            print "last optimized energy : ", optimizedenergy
            print "return code : ", returncode
            print "Function Call Counter : ", functionCallCounter
            #stepcounter += 1
            print " TERMINATING THE PROGRAM"
            break
        ##########################################################################################
        ## End Of Relaxation
        ##########################################################################################
        relaxationendtime = time.time()
        relaxation_time_length.append((relaxationendtime-relaxationstarttime)/60)#recording the time taken for relaxation
        functionCallCounterArray.append(functionCallCounter)#saving the number of function call done in each growth step
        print "Time taken for relaxation step ", growthcounter, ": ", relaxationendtime - relaxationstarttime, "   in minutes : ", (relaxationendtime - relaxationstarttime)/60., "   in Hours : ", (relaxationendtime - relaxationstarttime)/3600.
        print "Number of Function Calls done : ", functionCallCounter
        # Getting Coordinates after relaxation
        #np.save(r'coordinates_after_step={0:d}'.format(int(growthcounter)),xopt)
        currentcoordinates = np.copy(xopt)
        ############################################################################
        # Taking output from optimizer and plotting                                #
        ############################################################################
        optimizedenergy = opt.last_optimum_value()#the most optimized value of energy
        optimizationTermination(returncode)
        #Printing the Values
        print "Current Step : ",growthcounter, "returncode : ", returncode
        print "Optimized energyvalue : ", optimizedenergy
        print "##########################################################################################"
        print 	"Calculation of Parameters and Plotting"
        print "##########################################################################################"
        energyarray  = np.append(energyarray, optimizedenergy)
        #grabbing the terms of energy
        tempterms = cartesianTermsofenergy(cell,xopt)
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
        ###########################################################################################
        #   Putting New Vertex Coordinates In The Vertices On The CELL 
        ###########################################################################################
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
        cell.setParameters()
        ###################################################################################################
        ###                            PERFORMING GROWTH                                               ####
        ###################################################################################################
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "     Performing growth       STEP : : ", growthcounter
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ## INFLATED GROWTH :: Growth proportional to current shape and with guassian randomness
        #cell = feedbackInflatedGrowFaces(cell)
        ## STRAIN GROWTH :: Growth proportional to current STRAIN and GAUSSIAN Noise
        if args.exp: #if this argument is true, then Exponential growth
            print " Exponential Growth with Feedback = %f"%(eta)
            cell = feedbackInflatedGrowFaces(cell)
        else:# IF false Strain Growth
            print " Strain Growth with Feedback = %f"%(eta)
            cell = feedbackStrainGrowFaces(cell)
        #Plotting the Growth Rate Surface
        plotGrowthRateSurface(cell,numOfLayer,name='growthRateSurface=%03.d.png'%(growthcounter+1),vmax=vmax)
        ###############################
        #saving figures after growth
        ###############################
        saveGrowthFigures(cell,numOfLayer,growthcounter+1)
        ## Saving coordinates before relaxation of growth ##
        np.save('growth_images/growth_coordinates_before_step=%03.d.png'%growthcounter,xopt)#saving initial coordinates
        meandeterminantarray = plotMeanTargetArea(cell, meandeterminantarray)
        ###########################################################################
        ## Saving Coordinates & Target Form Matrix
        ###########################################################################
        #saving coordinates
        np.save(r'coordinates_step={0:d}'.format(int(growthcounter+1)),xopt)
        #saving the TargetFormMatrix of Cell
        saveTargetFormMatrix(cell,growthcounter+1)
        ###########################################################################
        ## Plot parameters after the relaxation above
        ###########################################################################
        plotParametersAfterRelaxation(cell,numOfLayer,energyarray, firsttermarray, secondtermarray, 
                                      thirdtermarray, volumearray, deformations, meanDeformation, stdDeformation, 
                                      fourthtermarray,relaxation_time_length, functionCallCounterArray)
        #######################################################################################################################
	matrixdictionary = plotMatrixDictionary(cell,faceidarray, matrixdictionary)
	#######################################################################################################################
	### ~~~~~~~~END OF LOOP~~~~~~~
	#######################################################################################################################        
print "=============================================================================="
print "The Growth with steps ",growthSteps, "done"
print "=============================================================================="
print "				The END"
print "=============================================================================="