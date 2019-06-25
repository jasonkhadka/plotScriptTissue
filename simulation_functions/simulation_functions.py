###################################################################################################################################
# Simulation Functions
# This library has functions needed to run the simulations of tissue growth using Quadedge
# @ 2017 April Jason Khadka
#
###################################################################################################################################
import numpy as np
#import pandas as 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#pd.set_option('precision',10)
#pd.set_option("display.max_columns",200)
#pd.set_option("display.max_rows",2000)
import sys 
sys.path.append("/home/jkhadka/plantdev")
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
sys.path.append('/home/jkhadka/transferdata/scripts/strain_plots/')
sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scripts/strain_plots/')
sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scripts/simulation_functions/')
import quadedge as qd
import ellipse as ep
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
import os
########################################################################################
# Formatter for ticks
########################################################################################
import matplotlib.ticker
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
	def __init__(self, order=0, fformat="%1.2f", offset=True, mathText=True):
		self.oom = order
		self.fformat = fformat
		matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
	def _set_orderOfMagnitude(self, nothing):
		self.orderOfMagnitude = self.oom
	def _set_format(self, vmin, vmax):
		self.format = self.fformat
		if self._useMathText:
			self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)


######################################################################
#          Function to make initial Dome                            ##
######################################################################
def makeDome(cell, numOfLayer, Length = 1.):
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
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
		a=l/radius
		z=radius*np.cos(a)
		if z>= 0.0:
			a = radius*np.sin(a)/l
		else:
			z = 0.5*np.pi*radius-l
			a = radius/l
		x*=a
		y*=a
		vertex.setXcoordinate(x)#setting x coordinate
		vertex.setYcoordinate(y)#setting y coordinate
		vertex.setZcoordinate(z)#setting z coordinate
		if zmin > z:
			zmin = z
		counter += 1
		vertex = vertices.next()
	#################################
	###resaling vertices To make it positive
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	counter = 0#counter to change vertex positions
	while vertex != None:
		z = vertex.getZcoordinate()
		vertex.setZcoordinate(z-zmin)
		vertex = vertices.next()
	###################################
	cell.setInitialParameters()
	return cell
######################################################################
#          Function to initial Cylinder                             ##
######################################################################
def makeCylinder(cell,numOfLayer, Length = 1.,laterCylinder=None):
	if laterCylinder==None:
		laterCylinder = numOfLayer
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
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
	#########################################################
	###Making the Cylidrical Flanks of the Dome 
	#########################################################
	cell = latdev.makeEdgeHexagonsPentagons(cell,numOfLayer)
	cell = latdev.makeFirstCylindricalLayer(cell,numOfLayer)
	for _ in range(laterCylinder):cell = latdev.makeLaterCylindricalLayer(cell,numOfLayer)
	########Moving the whole thing in postive z direction
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	zmin = 10
	#counter = 0#counter to change vertex positions
	while vertex != None:
		z = vertex.getZcoordinate()
		if zmin > z:
			zmin = z
		counter += 1
		vertex = vertices.next()
	#cell = settingParameters(cell)
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	while vertex != None:
		z = vertex.getZcoordinate()
		vertex.setZcoordinate(z-zmin)
		vertex = vertices.next()
	############################
	cell.setInitialParameters()
	return cell
############################################################################
###        Function to load cell from file
############################################################################
def loadCellFromFile(step,numOfLayer = 8,resetids = False, location = None):
	####################################
	if location:
		if not ((os.path.isfile(location+"qdObject_step=%03d.obj"%step)) or (os.path.isfile(location+"TargetFormMatrix_step=%d.npy"%step))):
			return 
		loadedcell = qd.objReadCell(location+"qdObject_step=%03d.obj"%step)
	else:
		if not ((os.path.isfile("qdObject_step=%03d.obj"%step)) or (os.path.isfile("TargetFormMatrix_step=%d.npy"%step))):
			return 
		loadedcell = qd.objReadCell("qdObject_step=%03d.obj"%step)
	if resetids:
		# Flipping the Face ID after Loading to cell so that the face id before and after matches
		faces = qd.CellFaceIterator(loadedcell)
		facecount = loadedcell.countFaces()
		face= faces.next()
		while face != None:
			faceid = face.getID()
			face.setID(facecount-faceid+1)
			#print face.getID()
			face = faces.next()
	loadedcell.setInitialParameters()
	######################################################
	#print "######################################################"
	#print "#"," "*10, "step %d"%step
	#print "######################################################"
	#settig target Form Matrix : 
	#TMF step corresponds to coordinate step
	####
	setTargetFormMatrix(loadedcell,step, location = location)
	return loadedcell
############################################################################
###        Function to change the stiffness around the perifery of 
###        target cell and its 1 next neighbour
############################################################################
def changePerimeterPrimodialAlpha(cell,faceid,changedAlpha):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
			if face.getID() == faceid:
					#print "faceid"
					face.setAlpha(changedAlpha)
					edges = qd.FaceEdgeIterator(face)
					edge = edges.next()
					while edge != None:
							rightFace = edge.Right()
							rightFace.setAlpha(changedAlpha)
							# iterating the face again
							rightedges = qd.FaceEdgeIterator(rightFace)
							rightedge = rightedges.next()
							while rightedge != None:
								nextFace = rightedge.Right()
								#print nextFace.getID()
								nextFace.setAlpha(changedAlpha)
								rightedge = rightedges.next()
							edge = edges.next()
					#print "kappa for this face : ", face.getKappa()
			face = faces.next()
	#resetting alpha for all inner cells
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	changedAlpha = cell.getAlpha()
	while (face != None):
			if face.getID() == faceid:
					#print "faceid"
					face.setAlpha(changedAlpha)
					edges = qd.FaceEdgeIterator(face)
					edge = edges.next()
					while edge != None:
							rightFace = edge.Right()
							rightFace.setAlpha(changedAlpha)
							edge = edges.next()
					#print "kappa for this face : ", face.getKappa()
			face = faces.next()
	return cell
############################################################################
###        Function to change the feedback around the  
###        target cell and its 1 next neighbour
############################################################################
def changePrimodialEta(cell,faceid,changedAlpha):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
			if face.getID() == faceid:
					#print "faceid"
					face.setEta(changedAlpha)
					edges = qd.FaceEdgeIterator(face)
					edge = edges.next()
					while edge != None:
							rightFace = edge.Right()
							rightFace.setEta(changedAlpha)
							# iterating the face again
							rightedges = qd.FaceEdgeIterator(rightFace)
							rightedge = rightedges.next()
							while rightedge != None:
								nextFace = rightedge.Right()
								#print nextFace.getID()
								nextFace.setEta(changedAlpha)
								rightedge = rightedges.next()
							edge = edges.next()
					#print "kappa for this face : ", face.getKappa()
			face = faces.next()
	return cell
############################################################################
###        Function to change threshold on boundary cells                 ##
############################################################################
def setBoundaryDivisionFactor(cell, targetid, factor):
	faces = getPrimordiaBoundaryFaceList(cell, targetid, large=False)
	for face in faces:
		face.setDivisionFactor(factor)
		face.setDivisionThreshold()
	return
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def randomDivideCell(cell):
	#global cell
	print "##################################################"
	print "Perfomring Random Cell Division"
	print "##################################################"
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		face.divideRandom()
		face = faces.next()
	cell.setParameters()
	return cell
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def oneRandomDivideCell(cell):
	#global cell
	print "##################################################"
	print "Performing Random Cell Division"
	print "##################################################"
	import random as rn
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	facedivisionarray = []
	while face != None:
		if face.getAreaOfFace() > face.getDivisionThreshold():#large enough to divide
			facedivisionarray.append(face.getID())
		face = faces.next()
	#################################################
	#now choose one face randomly to divide
	#################################################
	if len(facedivisionarray) == 0:#if no face crossing threshold
		return cell
	else:
		facetodivide = rn.choice(facedivisionarray)
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == facetodivide:#large enough to divide
				face.divideRandom()
				break
			face = faces.next()
		cell.setParameters()
	return cell
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def oneRadialDivideCell(cell, targetface):
	#global cell
	print "##################################################"
	print "Performing Radial Cell Division"
	print "##################################################"
	import random as rn
	##########################################################
	# calculating radial/orthoradial vec before division
	##########################################################
	cell.setRadialOrthoradialVector(targetface)
	##########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	facedivisionarray = []
	while face != None:
		if face.getAreaOfFace() > face.getDivisionThreshold():#large enough to divide
			facedivisionarray.append(face.getID())
		face = faces.next()
	#################################################
	#now choose one face randomly to divide
	#################################################
	if len(facedivisionarray) == 0:#if no face crossing threshold
		return cell
	else:
		facetodivide = rn.choice(facedivisionarray)
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == facetodivide:#large enough to divide
				face.divideRadial()
				break
			face = faces.next()
		cell.setParameters()
	return cell
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def oneOrthoradialDivideCell(cell, targetface):
	#global cell
	print "##################################################"
	print "Performing Orthoradial Cell Division"
	print "##################################################"
	import random as rn
	##########################################################
	# calculating radial/orthoradial vec before division
	##########################################################
	cell.setRadialOrthoradialVector(targetface)
	##########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	facedivisionarray = []
	while face != None:
		if face.getAreaOfFace() > face.getDivisionThreshold():#large enough to divide
			facedivisionarray.append(face.getID())
		face = faces.next()
	#################################################
	#now choose one face randomly to divide
	#################################################
	if len(facedivisionarray) == 0:#if no face crossing threshold
		return cell
	else:
		facetodivide = rn.choice(facedivisionarray)
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == facetodivide:#large enough to divide
				face.divideOrthoradial()
				break
			face = faces.next()
		cell.setParameters()
	return cell
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def oneMaximalStressDivideCell(cell):
	#global cell
	print "##################################################"
	print "Performing Maximal Stress Cell Division"
	print "##################################################"
	import random as rn
	################################################
	# calculating stress/strain before division
	################################################
	cell.calculateStressStrain()
	################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	facedivisionarray = []
	while face != None:
		if face.getAreaOfFace() > face.getDivisionThreshold():#large enough to divide
			facedivisionarray.append(face.getID())
		face = faces.next()
	#################################################
	#now choose one face randomly to divide
	#################################################
	if len(facedivisionarray) == 0:#if no face crossing threshold
		return cell
	else:
		facetodivide = rn.choice(facedivisionarray)
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == facetodivide:#large enough to divide
				face.divideMaximalStress()
				break
			face = faces.next()
		cell.setParameters()
	return cell
############################################################################
###        Function to Divide Cells of this Tissue                       ###
############################################################################
def oneShortAxisDivideCell(cell):
	#global cell
	print "##################################################"
	print "Performing Short Axis Cell Division"
	print "##################################################"
	import random as rn
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	facedivisionarray = []
	while face != None:
		if face.getAreaOfFace() > face.getDivisionThreshold():#large enough to divide
			facedivisionarray.append(face.getID())
		face = faces.next()
	#################################################
	#now choose one face randomly to divide
	#################################################
	if len(facedivisionarray) == 0:#if no face crossing threshold
		return cell
	else:
		facetodivide = rn.choice(facedivisionarray)
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while face != None:
			if face.getID() == facetodivide:#large enough to divide
				face.divideShortAxis()
				break
			face = faces.next()
		cell.setParameters()
	return cell
######################################################################
#           Get Coordinates of all the vertices                     ##
######################################################################
## Cartesian Coordinates ##
def getCartesianCoordinate(cell):
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
	return np.concatenate((xcoord,ycoord,zcoord))
## Cylindrical Coordinates ##
def getCylindricalCoordinate(cell):
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertexid = np.array([])
	xcoord = np.array([])
	ycoord = np.array([])
	zcoord = np.array([])
	while vertex != None: 
		#saving the ids
		vertexid = np.append(vertexid, vertex.getID())
		xcoord = np.append(xcoord, vertex.getRcoordinate())
		ycoord = np.append(ycoord, vertex.getThetacoordinate())
		zcoord = np.append(zcoord, vertex.getZcoordinate())
		vertex = vertices.next()
	return np.concatenate((xcoord,ycoord,zcoord))
##############################################################################
###  Change Stiffness of cells that have their x-cordinate of their center####
###          -0.5<= x <=0.5                                               ####
##############################################################################
def changeStiffnessXBand(cell,factor = 10):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID() == 1:
			face = faces.next()
			continue
		#Changing stiffness for cells with x cordinate of center with 0.5 of origin
		if abs(face.getXCentralised())<=0.5:
			face.setAlpha(cell.getAlpha()*factor)
		face = faces.next()
	return cell
##############################################################################
###  Get Vertex of given ID in the given Cell                             ####
##############################################################################
def  getVertex(cell,vertexid):
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	while (vertex != None):
		if vertex.getID() == vertexid:
			return vertex
		vertex = vertices.next()
	return None
##############################################################################
# Get TFM determinant
##############################################################################
def getTargetFormMatrixDeterminant(face):
	return qd.getTargetFormMatrix(face,0,0)*qd.getTargetFormMatrix(face,1,1)-qd.getTargetFormMatrix(face,1,0)*qd.getTargetFormMatrix(face,0,1)
########################################################################
##  Function to set TargetFormMatrix at random
##  * step 1 : reduce the All faces' form matrix by 20% 0.8*M
##  * step 2 : random increase or decrease the matrix by 
##             uniform random number in range -25% to +25%
########################################################################
def setRandomizedTargetFormMatrix(cell):
	#############################
	import random 
	random.seed(1000)
	#############################
	initialStrain = 0.8
	lowerLim = 0.75
	upperLim = 1.25
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() ==1 :
			face = faces.next()
			continue
		#getting TargetForm Matrix from the faces and applying initial strain & 
		#increasing or decreasing it by random noise
		randomnoise = random.uniform(lowerLim, upperLim)
		targetForm = [randomnoise*initialStrain*qd.getTargetFormMatrix(face,0,0),
					  qd.getTargetFormMatrix(face,0,1),
					  qd.getTargetFormMatrix(face,1,0),
					  randomnoise*initialStrain*qd.getTargetFormMatrix(face,1,1)]
		#now assigning it back to the face
		qd.setTargetFormMatrix(face,0,0,targetForm[0])
		qd.setTargetFormMatrix(face,0,1,targetForm[1])
		qd.setTargetFormMatrix(face,1,0,targetForm[2])
		qd.setTargetFormMatrix(face,1,1,targetForm[3])
		####################
		face = faces.next()
	return cell
##################################################################
###         BOUNDS for optimization   CARTESIAN     
### Boundary vertices are fixed in their position 
### Cylinder is checked : If cylinder flank vertices are only allowed
###                        Z-FREEDOM            ####
##################################################################
def getBoundsForOptimization(cell,radius, perimeterVertexNum=0, Length = 1):
	print "------------------------------------------------------"
	print "|","Fixed Bounds on the external vertices         ","|"
	print "|","Cylinder flanks are x-y fixed & only z freedom","|"
	print "------------------------------------------------------"
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
	epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
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
		##       DOME VERTEX    X##########
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = fac*radius
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = -fac*radius
			lboundsy[vertcounter] = -fac*radius
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter] + Length#z-coordinates are allowed to move by 1 side lenght up or down
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter]
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = zcoord[vertcounter]
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
##################################################################
###         BOUNDS for optimization   CARTESIAN     
### Boundary vertices are fixed in their position 
### Cylinder is Not checked
### Flank vertices are completely allowed to move 
##################################################################
def getBoundsForOptimizationWithoutCylinderBound(cell,radius, perimeterVertexNum=0, Length = 1):
	print "----------------------------------------------"
	print "|","Fixed Bounds on the external vertices","|"
	print "----------------------------------------------"
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
	epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
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
		##  All vertices are allowed to vary in X-Y-Z within certain limit##########
		uboundsx[vertcounter] = fac*radius
		uboundsy[vertcounter] = fac*radius
		uboundsz[vertcounter] = fac*radius
		lboundsx[vertcounter] = -fac*radius
		lboundsy[vertcounter] = -fac*radius
		lboundsz[vertcounter] = 0.
		vertcounter += 1
		vertex = vertices.next()
	##################################################################
	###     base cells          ###
	##################################################################
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter]
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = zcoord[vertcounter]
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
##################################################################
###         BOUNDS for optimization :  CARTESIAN     
### Boundary vertices have 5% of average side length freedome ####
##################################################################
def getEpsBaseBoundsForOptimization(cell,radius, numOfLayer, perimeterVertexNum=0, Length = 1):
	print "-----------------------------------------------------"
	print "|",r"5 % Freedom on Bounds on the external vertices","|"
	print "-----------------------------------------------------"
	perimeterVertexNum = 6+12*(numOfLayer-1)  # number of Vertices in perimeter
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
	epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
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
		##       DOME VERTEX    X##########
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = fac*radius
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = -fac*radius
			lboundsy[vertcounter] = -fac*radius
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			print "cylinder"
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter] + Length#z-coordinates are allowed to move by 1 side lenght up or down
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	### First getting initial boundary vertex positions for this layer of cells ###
	##########################################
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	"""
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			thisvertex = getVertex(newcell,vertex.getID())#getting vertex with same vertex ID
			uboundsx[vertcounter] = thisvertex.getXcoordinate() + epsbasefreedom
			uboundsy[vertcounter] = thisvertex.getYcoordinate() + epsbasefreedom
			uboundsz[vertcounter] = thisvertex.getZcoordinate() + epsbasefreedom
			lboundsx[vertcounter] = thisvertex.getXcoordinate() - epsbasefreedom
			lboundsy[vertcounter] = thisvertex.getYcoordinate() - epsbasefreedom
			lboundsz[vertcounter] = thisvertex.getZcoordinate() - epsbasefreedom
		vertex = vertices.next()
		vertcounter += 1
	"""
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter] + epsbasefreedom
			uboundsy[vertcounter] = ycoord[vertcounter] + epsbasefreedom
			uboundsz[vertcounter] = zcoord[vertcounter] + epsbasefreedom
			lboundsx[vertcounter] = xcoord[vertcounter] - epsbasefreedom
			lboundsy[vertcounter] = ycoord[vertcounter] - epsbasefreedom
			lboundsz[vertcounter] = zcoord[vertcounter] - epsbasefreedom
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
##################################################################
###         BOUNDS for optimization   CARTESIAN     
### Boundary vertices are fixed in their position             ####
##################################################################
def getBoxBoundsForOptimization(cell,radius, perimeterVertexNum=0, Length = 1.):
	print "----------------------------------------------"
	print "|","Box Bounds on the external vertices","|"
	print "----------------------------------------------"
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
	epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
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
		##       DOME VERTEX    X##########
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = fac*radius
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = -fac*radius
			lboundsy[vertcounter] = -fac*radius
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter] + Length#z-coordinates are allowed to move by 1 side lenght up or down
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
		vertcounter += 1
		vertex = vertices.next()
	###############################
	###     base cells          ###
	###############################
	# First : Getting the z-max for the boudary cells
	zmax = 0.
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			if zcoord[vertcounter]>zmax:
				zmax = zcoord[vertcounter]
		vertex = vertices.next()
		vertcounter += 1
	##### Now making the bounds for Base Cells #####
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter] + np.sqrt(3.)*Length
			uboundsy[vertcounter] = ycoord[vertcounter] + np.sqrt(3.)*Length
			uboundsz[vertcounter] = zmax#getting maximal z value of the boundary as z upper bound
			lboundsx[vertcounter] = xcoord[vertcounter] - np.sqrt(3.)*Length
			lboundsy[vertcounter] = ycoord[vertcounter] - np.sqrt(3.)*Length
			lboundsz[vertcounter] = 0.#Low bound of Hexagon
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
######################################################################################
###         Plane BOUNDS for optimization   CARTESIAN     
### Boundary vertices are fixed in their plane     
###   That means -- Only X-Y values are allowed to cahnge for boundary vertex
######################################################################################
def getPlaneBoundsForOptimization(cell,radius, perimeterVertexNum=0, Length = 1.):
	print "Plane Bounds on the external vertices"
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
	epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
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
		##       DOME VERTEX    X##########
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = fac*radius
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = -fac*radius
			lboundsy[vertcounter] = -fac*radius
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter] + Length#z-coordinates are allowed to move by 1 side lenght up or down
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter] + np.sqrt(3)*Length#x-y plane freedom
			uboundsy[vertcounter] = ycoord[vertcounter] + np.sqrt(3)*Length#x-y plane freedom
			uboundsz[vertcounter] = zcoord[vertcounter]# z-bounds are fixed
			lboundsx[vertcounter] = xcoord[vertcounter] - np.sqrt(3)*Length#x-y plane freedom
			lboundsy[vertcounter] = ycoord[vertcounter] - np.sqrt(3)*Length#x-y plane freedom
			lboundsz[vertcounter] = zcoord[vertcounter]
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
####################################################################
### Strong Bottom BOUNDS for optimization : CYLINDRICAL Coord   ####
###  Only the bottom vertices are bound. 
###  The bound : 
###             Z : Bounded
###             R : Bounded
###             Theta : Bounded   
####################################################################
def getCylindricalStrongBottomBoundsForOptimization(cell,radius, perimeterVertexNum, Length = 1.):
	#getting the coordinates of the cell
	##intial parameters for optimization, the coordiantes of vertices
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	xcoord = np.array([])
	ycoord = np.array([])
	zcoord = np.array([])
	while vertex != None: 
		#saving the ids
		#vertexid = np.append(vertexid, vertex.getID())
		xcoord = np.append(xcoord, vertex.getRcoordinate())
		ycoord = np.append(ycoord, vertex.getThetacoordinate())
		zcoord = np.append(zcoord, vertex.getZcoordinate())
		vertex = vertices.next()
	coordinates = np.concatenate((xcoord,ycoord,zcoord))
	numberofvertices = cell.countVertices()
	#print "lenght of xcoord : ", len(xcoord), "totalvertex", totalvertex
	lowRBound = radius-0.01*radius#in this case: we just let the vertices in the boundary move freely in their own X-Y plane
	upRBound = radius+0.01*radius
	#epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
	print "radius : ", radius, "Low bound : ",lowRBound, "Up Bound: ", upRBound
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
		##       DOME VERTEX    X##########
		## X => R  -> [0,inf]
		## Y => Theta -> [-pi,pi]
		## Z => Z  -> [-inf,inf]
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = np.pi#theta
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = 0.
			lboundsy[vertcounter] = -1.*np.pi#theta
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = fac*radius#R
			uboundsy[vertcounter] = np.pi#THETA
			uboundsz[vertcounter] = fac*radius#z-coordinates are allowed to move by 1 side lenght up or down
			####LOWER BOUNDS###
			lboundsx[vertcounter] = 0.#R
			lboundsy[vertcounter] = -1.*np.pi#THETA
			#lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
			lboundsz[vertcounter] = 0.#z minimum = 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external (BOTTOM) vertex then new bounds applied
			uboundsx[vertcounter] = xcoord[vertcounter]#R - Fixed
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter]#Z - Fixed
			lboundsx[vertcounter] = xcoord[vertcounter]#R - Fixed
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = zcoord[vertcounter]#Z - Fixed
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
##################################################################
###      Bottom BOUNDS for optimization : CYLINDRICAL Coord   ####
###  Only the bottom vertices are bound. 
###  The bound is in X-Y plane, that is only z-coordinate is fixed
##################################################################
def getCylindricalBottomBoundsForOptimization(cell,radius, perimeterVertexNum, Length = 1.):
	#getting the coordinates of the cell
	##intial parameters for optimization, the coordiantes of vertices
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	xcoord = np.array([])
	ycoord = np.array([])
	zcoord = np.array([])
	while vertex != None: 
		#saving the ids
		#vertexid = np.append(vertexid, vertex.getID())
		xcoord = np.append(xcoord, vertex.getRcoordinate())
		ycoord = np.append(ycoord, vertex.getThetacoordinate())
		zcoord = np.append(zcoord, vertex.getZcoordinate())
		vertex = vertices.next()
	coordinates = np.concatenate((xcoord,ycoord,zcoord))
	numberofvertices = cell.countVertices()
	#print "lenght of xcoord : ", len(xcoord), "totalvertex", totalvertex
	lowRBound = radius-0.01*radius#in this case: we just let the vertices in the boundary move freely in their own X-Y plane
	upRBound = radius+0.01*radius
	#epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
	print "radius : ", radius, "Low bound : ",lowRBound, "Up Bound: ", upRBound
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
		##       DOME VERTEX    X##########
		## X => R  -> [0,inf]
		## Y => Theta -> [-pi,pi]
		## Z => Z  -> [-inf,inf]
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = np.pi#theta
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = 0.
			lboundsy[vertcounter] = -1.*np.pi#theta
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = fac*radius#R
			uboundsy[vertcounter] = np.pi#THETA
			uboundsz[vertcounter] = fac*radius#z-coordinates are allowed to move by 1 side lenght up or down
			####LOWER BOUNDS###
			lboundsx[vertcounter] = 0.#R
			lboundsy[vertcounter] = -1.*np.pi#THETA
			#lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
			lboundsz[vertcounter] = 0.#z minimum = 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external (BOTTOM) vertex then new bounds applied
			uboundsx[vertcounter] = fac*radius#R
			uboundsy[vertcounter] = np.pi#THETA
			uboundsz[vertcounter] = zcoord[vertcounter]
			lboundsx[vertcounter] = 0.#R
			lboundsy[vertcounter] = -1.*np.pi#THETA
			lboundsz[vertcounter] = zcoord[vertcounter]
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
##################################################################
###         BOUNDS for optimization   CYLINDRICAL             ####
##################################################################
def getCylindricalBoundsForOptimization(cell,radius, perimeterVertexNum, Length = 1.):
	#getting the coordinates of the cell
	##intial parameters for optimization, the coordiantes of vertices
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	xcoord = np.array([])
	ycoord = np.array([])
	zcoord = np.array([])
	while vertex != None: 
		#saving the ids
		#vertexid = np.append(vertexid, vertex.getID())
		xcoord = np.append(xcoord, vertex.getRcoordinate())
		ycoord = np.append(ycoord, vertex.getThetacoordinate())
		zcoord = np.append(zcoord, vertex.getZcoordinate())
		vertex = vertices.next()
	coordinates = np.concatenate((xcoord,ycoord,zcoord))
	numberofvertices = cell.countVertices()
	#print "lenght of xcoord : ", len(xcoord), "totalvertex", totalvertex
	lowRBound = radius-0.01*radius#in this case: we just let the vertices in the boundary move freely in their own X-Y plane
	upRBound = radius+0.01*radius
	#epsbasefreedom = (1./20)*(2.*np.pi*radius)/perimeterVertexNum#vertex is allowed to move 5% of average perimeter side length
	#print "small freedome for base vertices : 20%"
	#print "basefreedom", basefreedom
	#print "epsbasefreedom", epsbasefreedom
	print "radius : ", radius, "Low bound : ",lowRBound, "Up Bound: ", upRBound
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
		##       DOME VERTEX    X##########
		## X => R  -> [0,inf]
		## Y => Theta -> [-pi,pi]
		## Z => Z  -> [-inf,inf]
		if vertex.getDomePosition():#If domePosition == True -> then Vertex is in Dome
			uboundsx[vertcounter] = fac*radius
			uboundsy[vertcounter] = np.pi#theta
			uboundsz[vertcounter] = fac*radius
			lboundsx[vertcounter] = 0.
			lboundsy[vertcounter] = -1.*np.pi#theta
			lboundsz[vertcounter] = 0.
		else:## #If domePosition == False -> then Vertex is in Cylindrical Flanks
			##     Cylinder VERTEX    X##########
			uboundsx[vertcounter] = upRBound#R
			uboundsy[vertcounter] = np.pi#THETA
			#uboundsz[vertcounter] = zcoord[vertcounter] + Length#z-coordinates are allowed to move by 1 side lenght up or down
			uboundsz[vertcounter] = fac*radius#z-coordinates are allowed to move by 1 side lenght up or down
			####LOWER BOUNDS###
			lboundsx[vertcounter] = lowRBound#R
			lboundsy[vertcounter] = -1.*np.pi#THETA
			#lboundsz[vertcounter] = ((zcoord[vertcounter] - Length)>0)*(zcoord[vertcounter] - Length)#z-coordinates are allowed to move by 1 side length up or down or until z >= 0
			lboundsz[vertcounter] = 0.#z minimum = 0
		vertcounter += 1
		vertex = vertices.next()
	###base cells###
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	vertcounter = 0
	while vertex != None: 
		currentvertid = vertex.getID()
		if checkExternalVertex(vertex):#if the vertex is external vertex then new bounds applied
			"""uboundsx[vertcounter] = xcoord[vertcounter]
			uboundsy[vertcounter] = ycoord[vertcounter]
			uboundsz[vertcounter] = zcoord[vertcounter]
			lboundsx[vertcounter] = xcoord[vertcounter]
			lboundsy[vertcounter] = ycoord[vertcounter]
			lboundsz[vertcounter] = zcoord[vertcounter]"""
			uboundsx[vertcounter] = xcoord[vertcounter]#R
			uboundsy[vertcounter] = ycoord[vertcounter]#theta
			uboundsz[vertcounter] = zcoord[vertcounter]#Z
			lboundsx[vertcounter] = xcoord[vertcounter]#R
			lboundsy[vertcounter] = ycoord[vertcounter]#Theta
			lboundsz[vertcounter] = zcoord[vertcounter]#Z
		vertex = vertices.next()
		vertcounter += 1
	### Upper bounds ###
	upperbounds = np.concatenate((uboundsx,uboundsy,uboundsz))
	##Lower bounds
	lowerbounds = np.concatenate((lboundsx,lboundsy,lboundsz))
	return (upperbounds, lowerbounds)
####################################################################
#~Function : settingParameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~Task: Sets all the parameters on the cell to run simulation ~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~>>INPUT<<~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~cell - Quadedge cell~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~>>>>RETURN<<<<~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~cell - returns the cell after updating all the parameters~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################################################
def settingFirstParameters(cell):
	#setting projected coordiantes and other parameters
	for _ in range(1):
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setProjectedCoordinate()
			face1 = faces.next()
		#####################################################
		#setting parameters for all vertices
		vertices = qd.CellVertexIterator(cell)
		vertex = vertices.next()
		while vertex != None:
			vertex.setparameters()#setting Alpha, beta, Gamma and F1,f2,f3 and Ak for vertices
			vertex = vertices.next()
		#####################################################
		#setting => this gives the CurrentFormMatrix and Mu-Matrix (strain Matrix) is only set as CurrentFormMatrix as TargetForm is 0 Matrix
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setMu()
			face1 = faces.next()
		#####################################################
		#setting  TargetFormMatrix = CurrentFormMatrix for initial condition
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			face1.setTempTargetFormMatrixCurrent()
			face1 = faces.next()
		#####################################################
		#setting => CurrentFormMatrix with comparision to TargetFormMatrix
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setMu()
			face1 = faces.next()
		#####################################################
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setAreaOfFace()
			face1 = faces.next()
		"""
		#setting parameters for all vertices
		vertices = qd.CellVertexIterator(cell)
		vertex = vertices.next()
		while vertex != None:
			vertex.setDerivatives()
			vertex = vertices.next()
		"""
		#####################################################
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setEnergyTerms()
			face1 = faces.next()
	return cell
###################################################################################################################################
def settingParameters(cell):
	#setting projected coordiantes and other parameters
	for _ in range(1):
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setProjectedCoordinate()
			face1 = faces.next()
		#####################################################
		#setting parameters for all vertices
		vertices = qd.CellVertexIterator(cell)
		vertex = vertices.next()
		while vertex != None:
			vertex.setparameters()#setting Alpha, beta, Gamma and F1,f2,f3 and Ak for vertices
			vertex = vertices.next()
		#####################################################
		#setting => this gives the CurrentFormMatrix and Mu-Matrix (strain Matrix) is only set as CurrentFormMatrix as TargetForm is 0 Matrix
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setMu()
			face1 = faces.next()
		#####################################################
		#####################################################
		#####################################################
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setAreaOfFace()
			face1 = faces.next()
		"""
		#setting parameters for all vertices
		vertices = qd.CellVertexIterator(cell)
		vertex = vertices.next()
		while vertex != None:
			vertex.setDerivatives()
			vertex = vertices.next()
		"""
		#####################################################
		faces = qd.CellFaceIterator(cell)
		face1 = faces.next()
		while face1 != None:
			#print face1.getID()
			face1.setEnergyTerms()
			face1 = faces.next()
	return cell
###################################################################################################################################

###################################################################################################################################
############################################################################
###                    Function to Print all the values of Cell          ###
############################################################################
#print all the energy and matrix for faces
def printmatrixenergy(cell):
	faces = qd.CellFaceIterator(cell)
	face1 = faces.next()
	while face1 != None:
		#print face1.getID()
		print "*****************************************"
		face1.printTargetFormMatrix()
		print "Face ID : ", face1.getID()
		print "first term : ", face1.getFirstTerm()
		print "second term : ", face1.getSecondTerm()
		print "third term : ", face1.getThirdTerm()
		print "energy of face :", face1.getEnergy()
		face1 = faces.next()
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print "total Energy : ", cell.getEnergy()
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	return
###################################################################################################################################
###################################################################################################################################

###################################################################################################################################
############################################################################
#           to call at the end of program to save the coordiantes         ##
############################################################################
def endofprogram(txopt, totalenergy, firsttermarray, secondtermarray, thirdtermarray, volumearray,fourthtermarray, error = False):
	#saving the values at the end of program
	tempcoordatearray = np.array(txopt)
	if not error:
		np.save("final_coordinates",tempcoordatearray)
		#energyterms = np.vstack((firsttermarray, secondtermarray, thirdtermarray))
		np.save("energyterms.txt", np.c_[totalenergy, firsttermarray, secondtermarray, thirdtermarray, volumearray, fourthtermarray])
	else:
		np.save("ERROR_final_coordinates",tempcoordatearray)
		#energyterms = np.vstack((firsttermarray, secondtermarray, thirdtermarray))
		np.save("ERROR_energyterms", np.c_[totalenergy, firsttermarray, secondtermarray, thirdtermarray, volumearray,fourthtermarray])
	return
###################################################################################################################################
###################################################################################################################################
### Set Cylinder Alpha
###################################################################################################################################
def setCylinderAlpha(cell,cyaplha):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if not face.getDomePosition():#If domePosition == True -> then Vertex is in Dome ELSE it is in Cylinder
			face.setAlpha(cyaplha)
		face = faces.next()
	return cell
###################################################################################################################################
###############################################################################################################
###            Save all the images to Growth Folder
###############################################################################################################
def growFaces(cell):
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.grow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
###############################################################################################################
###            Inflated growth
###############################################################################################################
def inflatedGrowFaces(cell):
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.inflatedGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
###############################################################################################################
###           Feedback growth
###############################################################################################################
#### GROWTH WITH INFLATION -> Proportinal to a Intrinsic form ####
def feedbackInflatedGrowFaces(cell):
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackInflatedGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
#### GROWTH WITH STRAIN -> Proportional to Strain on the cell ####
def feedbackStrainGrowFaces(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackStrainGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
#### GROWTH WITH STRAIN Limiting Feedback -> Proportional to Strain on the cell ####
def feedbackLimitingStrainGrowFaces(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackLimitingStrainGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
#### GROWTH WITH STRAIN and Target matrix -> Proportional to Strain and on Target Matrix on the cell ####
def feedbackStrainProportionalGrowFaces(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackStrainProportionalGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
#### GROWTH WITH STRAIN, RANDOM DIRECRION and Target matrix -> Proportional to Strain and on Target Matrix on the cell ####
def feedbackRandomStrainProportionalGrowFaces(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackRandomStrainProportionalGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
###############################################################################################################
###           Feedback growth : Except the outer cells
###############################################################################################################
### GROWTH WITH STRAIN -> Proportional to Strain on the cell ####
def feedbackStrainGrowFacesWithoutOuter(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	print "Growing Face : "
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		# Skipping the faces that are outer
		if checkExternalFace(face):
			face = faces.next()
			continue
		print "Face ID :", face.getID()
		face.feedbackStrainGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
#### GROWTH WITH STRAIN and Target matrix -> Proportional to Strain and on Target Matrix on the cell ####
def feedbackStrainProportionalGrowFacesWithoutOuter(cell): 
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		# Skipping the faces that are outer
		if checkExternalFace(face):
			face = faces.next()
			continue
		face.feedbackStrainProportionalGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
###############################################################################################################
###           Feedback growth : Except the outer cells
###############################################################################################################
#### GROWTH WITH Constant -> Proportional to a constant intrinsic rate on the cell ####
def feedbackConstantGrowFaces(cell):
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		face.feedbackConstantGrow()
		face = faces.next()
	#now setting new parameters on the cell
	#cell = settingParameters(cell)
	cell.setParameters()
	return cell
###################################################################################################################################
###############################################################################################################
###            Save All the TargetFormMatrix to a folder
###             Matrix is saved as a dictionary {faceid: matrix}
###############################################################################################################
def saveTargetFormMatrix(cell,growthcounter):
	#getting all the target form matrix
	targetformmatrixDictionary = {}#dictionary to save all the targetformmatrix
	## Getting and saving matrix in array
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	matrixArray = np.zeros((2,2))
	#starting Iteration
	while face != None:
		faceid = face.getID()
		matrixArray[0,0] = qd.getTargetFormMatrix(face, 0,0)
		matrixArray[0,1] = qd.getTargetFormMatrix(face, 0,1)
		matrixArray[1,0] = qd.getTargetFormMatrix(face, 1,0)
		matrixArray[1,1] = qd.getTargetFormMatrix(face, 1,1)
		#saving this in dictionary
		#print matrixArray
		targetformmatrixDictionary[faceid] = np.copy(matrixArray)
		face = faces.next()
	#saving the dictionary
	np.save("TargetFormMatrix_step=%d"%growthcounter,targetformmatrixDictionary)
	return
###############################################################################################################
###         Loading the TargetFormMatrix from file and setting it to the faces for given growthstep
###############################################################################################################
def setTargetFormMatrix(cell, growthcounter, location = None):
	### Loading the TargetFormMatrix
	if location:
		loadedDictionary = np.load(location+"TargetFormMatrix_step=%d.npy"%growthcounter).item()
	else:
		loadedDictionary = np.load("TargetFormMatrix_step=%d.npy"%growthcounter).item()
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		faceid = face.getID()
		matrixArray = loadedDictionary[faceid]
		for i in range(2):
			for j in range(2):  
				qd.setTargetFormMatrix(face, i,j, matrixArray[i][j])
		face = faces.next()
	cell.setParameters()
	######
	return cell
################################################################################################
# Get Strain Matrix of target face
################################################################################################
def getStrainMatrix(cell,targetid):
	cell.calculateStressStrain()# Stress-Strain calculation
	faces = qd.CellFaceIterator(cell)
	#performing growth in the face
	face = faces.next()
	while face != None:
		if face.getID() != targetid:
			face = faces.next()
			continue
		strainmatrix = np.array([[face.getStrainValue(0,0),face.getStrainValue(0,1)],
									[face.getStrainValue(1,0),face.getStrainValue(1,1)]])
		break
	return strainmatrix
###############################################################################################################
###         Loading the Cylindrical coordinate from file and setting it to the faces for given growthstep
###############################################################################################################
def setCylindricalCoordinate(cell, growthcounter):
	### Loadign the coordinates
	numberofvertices = cell.countVertices()
	coordinates=np.load("coordinates_step=%d.npy"%(growthcounter))
	newvertexpositions = coordinates.reshape((3,numberofvertices))
	#setting the vertices of cell with new coordinates
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	counter = 0#counter to change vertex positions
	while vertex != None:
		vertex.setRcoordinate(newvertexpositions[0,counter])#setting x coordinate
		vertex.setThetacoordinate(newvertexpositions[1,counter])#setting y coordinate
		vertex.setZcoordinate(newvertexpositions[2,counter])#setting z coordinate
		counter += 1
		vertex = vertices.next()
	cell.setCartesian()
	cell.setParameters()
	return cell
###############################################################################################################
###         Loading the Cartesian coordinate from file and setting it to the faces for given growthstep
###############################################################################################################
def setCartesianCoordinate(cell, growthcounter):
	### Loadign the coordinates
	numberofvertices = cell.countVertices()
	coordinates=np.load("coordinates_step=%d.npy"%(growthcounter))
	newvertexpositions = coordinates.reshape((3,numberofvertices))
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
	return cell
###############################################################################################################
###         Loading the Cylindrical coordinate from file and setting it to the faces for given growthstep
###############################################################################################################
def setTFMMatchedCartesianCoordinate(cell, numOfLayer,gamma, z):
	### Loadign the coordinates
	CFFileName = r"/home/jkhadka/transferdata/coordinates_storage/matching-TFM/relaxed/z=%d/layer=%d/g=%.4fCF.npy"%(z,numOfLayer,gamma)
	TFFileName = r"/home/jkhadka/transferdata/coordinates_storage/matching-TFM/relaxed/z=%d/layer=%d/g=%.4fTF.npy"%(z,numOfLayer,gamma)
	#Loading the Coordinates
	numberofvertices = cell.countVertices()
	coordinates=np.load(CFFileName)
	newvertexpositions = coordinates.reshape((3,numberofvertices))
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
	### Loading the TargetFormMatrix
	cell.setInitialParameters()
	loadedDictionary = np.load(TFFileName).item()
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		faceid = face.getID()
		matrixArray = loadedDictionary[faceid]
		for i in range(2):
			for j in range(2):  
				qd.setTargetFormMatrix(face, i,j, matrixArray[i][j])
		face = faces.next()
	cell.setParameters()
	return cell
###############################################################################################################
###                 Function to return the Initial shape for  Matrix Dictionary of Cell
###         Matrix Dictionary : the dictionary that stores DETERMINANT for each Face of following : 
###                                 a.  Mu-Matrix
###                                 b.  Current Form-Matrix
###                                 c.  Target Form-Matrix
###############################################################################################################
def getInitialMatrixDictionary(cell,targetfaceid=0):
	faceidarray = []#array to store faceids
	matrixdictionary = {}#dictionary storing the matrices determinant [[mu],[CF],[TF]]
	faces =qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID()==1:
			face = faces.next()
			continue
		if face.getID()%10 == 0 or face.getID() == targetfaceid:
			faceidarray.append(face.getID())
			matrixdictionary[face.getID()] = [[],[],[]]
		face = faces.next()
	return matrixdictionary, faceidarray
###################################################################################################################################
###############################################################################################################
###            Save all the images to Growth Folder
###############################################################################################################
def saveGrowthFigures(cell,numOfLayer, growthcounter):
	#saving the images of cell before the after growthcounter steps
	latdev.plot3DCell(cell,name='growth_images/Tissue_before_step=%03.d.png'%growthcounter)
	#plotting the surface only of cell
	latdev.plotSurface(cell, numOfLayer,name = 'growth_images/Tissue_surface_before_step=%03.d.png'%growthcounter)
	return
###############################################################################################################
###            Print the reason for Termination of Optmizaion
###############################################################################################################
def optimizationTermination(returnvalue):
	returnValueDictionary = {1 : "NLOPT_SUCCESS : Generic success return value", 
							 2 : "NLOPT_STOPVAL_REACHED : Optimization stopped because stopval was reached.",
							 3 : "NLOPT_FTOL_REACHED : Stopped because ftol_rel or ftol_abs was reached.",
							 4 : "NLOPT_XTOL_REACHED : Optimization stopped because xtol_rel or xtol_abs was reached.",
							 5 : "NLOPT_MAXEVAL_REACHED : Optimization stopped because maxeval was reached.",
							 6 : "NLOPT_MAXTIME_REACHED : Optimization stopped because maxtime was reached.",
							'-1' : "NLOPT_FAILURE :Generic failure code. ",
							'-2' : "NLOPT_INVALID_ARGS :Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera)",
							'-3' : "NLOPT_OUT_OF_MEMORY : Ran out of memory. ",
							'-4' : "NLOPT_ROUNDOFF_LIMITED : Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.) ", 
							'-5' : "NLOPT_FORCED_STOP : Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimizations nlopt_opt object opt from the users objective function or constraints"
							}
	if returnvalue in returnValueDictionary:
		print returnValueDictionary[returnvalue]
	else:
		print "Return Value : ",returnvalue," and associated Termination Code Not Found!"
	return
##########################################################################################################################################################
############################################################################
###                    Return Energy Values                              ###
############################################################################
def cartesianTermsofenergy(cell, inputcoordinates):
	#global maincounter, numOfLayer, alpha, beta, pressure, cell
	#making of hexagonal lattice
	#maincounter += 1 
	#Reshape the tempcoordinates, to x-y-z arrays
	numberofvertices = cell.countVertices()
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
	#cell = settingParameters(cell)
	cell.setParameters()
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
############################################################################
###       Return Energy Values : Takes CYLINDRICAL COORDINATES           ###
############################################################################
def cylinderTermsofenergy(cell, inputcoordinates):
	#global maincounter, numOfLayer, alpha, beta, pressure, cell
	#making of hexagonal lattice
	#maincounter += 1 
	#Reshape the tempcoordinates, to x-y-z arrays
	numberofvertices = cell.countVertices()
	tempcoordinates = inputcoordinates.reshape((3,numberofvertices))
	#to store the deformation on the cells
	faceids = []
	deformations = []
	####iterating the vertices to feed the new coordinates in
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	counter = 0#counter to change vertex positions
	while vertex != None:
		vertex.setRcoordinate(tempcoordinates[0,counter])#setting x coordinate
		vertex.setThetacoordinate(tempcoordinates[1,counter])#setting y coordinate
		vertex.setZcoordinate(tempcoordinates[2,counter])#setting z coordinate
		counter += 1
		vertex = vertices.next()
	#####################################################
	#cell = settingParameters(cell)
	cell.setCartesian()#setting the cartesian coordinates from the input Cylindrical
	cell.setParameters()
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
##############################################################################################################################
############################################################################
###                    Current Mean Determinant Target Area Plot         ###
############################################################################
def plotMeanTargetArea(cell,meandeterminantarray):
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
	timearray = range(len(meandeterminantarray))
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

################################################################################################################################################
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
################################################################################################################################################

#############################################################################################
### Function to take Plot the determinants Current Form Matrix, Target Form Matrix and Strain Matrix
### of all given faceid found in matrixdictionary (storange for determinants for all matrices)
###             growthStep : Current growthStep
##############################################################################################
def plotMatrixDictionary(cell, faceidarray, matrixdictionary):
	for fid in faceidarray:
			matrixdictionary[fid] = plotMatrices(cell,fid,matrixdictionary[fid])
	return matrixdictionary
#############################################################################################
### Plotting the CurrentFormMatrix, TargetFormMatrix and StrainMatrix
### Plot is made determinant of matrix vs GrowthtimeArray
##############################################################################################
def plotMatrices(cell, targetfaceid,matrixarray):
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
	growthStep = len(mudeterminant)## Each growth step adds 1 element to the determinant array, so total growthstep = total length
	plt.figure(21)
	plt.title("Determinant of Matrices for Face : %d"%targetfaceid)
	plt.plot(range(growthStep), mudeterminant,'-s',markeredgewidth=0.0,markeredgecolor=None,color='b', label="Mu")
	plt.plot(range(growthStep), currentFormdeterminant,'-o',markeredgewidth=0.0,markeredgecolor=None,color='g', label="CF")
	plt.plot(range(growthStep), targetFormdeterminant,'-v',markeredgewidth=0.0,markeredgecolor=None, color='r',label="TF")
	plt.legend(loc='best')
	plt.savefig('matrix_plot_faceid_%d.png'%targetfaceid,transparent=True)
	plt.clf()
	plt.close('all')
	######
	return [mudeterminant,currentFormdeterminant,targetFormdeterminant]
##########################################################################################################################################################
###############################################################################################################
###            Function to check if a given vertex is external vertex(vertex on boundary or not)
###############################################################################################################
def checkExternalVertex(vertex):
	edges = qd.VertexEdgeIterator(vertex)
	edge = edges.next()
	while edge != None:
		if edge.Left().getID() == 1 :
			return True
		edge = edges.next()
	return False
##########################################################################################
#       Function to check if given face is adjacent to outer face
##########################################################################################
def checkExternalFace(face):
	edges = qd.FaceEdgeIterator(face)
	edge = edges.next()
	while edge != None:
		adface = edge.Right()
		if adface.getID() == 1:
			return True
		edge = edges.next()
	return False
##########################################################################################################################################################
###############################################################################################################
###      function to print all the parameters
###############################################################################################################
def printCellParameters(cell,numOfLayer, Length = 1.):
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length
	print "Layer Number : ", numOfLayer
	print "Number of Faces : ", cell.countFaces()
	print "Calculated : Radius = ", radius
	print "------------------------------------------------"
	print "              Parameters"
	print "------------------------------------------------"
	print "alpha : ", cell.getAlpha()
	print "beta : ", cell.getBeta()
	print "pressure : ", cell.getPressure()
	print "zeta : ", cell.getZeta()
	print "gamma : ",cell.getGamma()
	print "Kappa : ", cell.getKappa()
	print "Feedback (Eta) :", cell.getEta()
	print "Initial Strain : ", cell.getInitialStrain()
	print "Sigma : ", cell.getSigma()
	return


##########################################################################################################################################################
###############################################################################################################
###      Plotting Function to plot all the plotting after relaxation has been done
###############################################################################################################
###########################################################################################
def plotParametersAfterRelaxation(cell,numOfLayer,energyarray, firsttermarray, secondtermarray, 
		thirdtermarray, volumearray, deformations, meanDeformation, stdDeformation, 
		fourthtermarray,relaxation_time_length, 
		functionCallCounterArray, divisionCounter, 
		sumWallLength, helfrichEnergy,meanfacearea):
		#time array for x-axis of plot
		import matplotlib.ticker as ticker
		print " STARTING TO PLOT THE PARAMETERS"
		growthstep = len(energyarray)
		timearray = range(growthstep)#each time step = new addition to energyarray
		#plotting the cell
		latdev.plot3DCell(cell,name='time=%03.d.png'%(growthstep-1))
		#plotting the surface only of cell
		latdev.plotSurface(cell, numOfLayer,name = 'surface_time=%03.d.png'%(growthstep-1))
		#plotting the energy plot 
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("total energy")
		plt.ylabel("optimized energy values")
		plt.xlabel("time")
		ax.plot(timearray, energyarray,'--o')
		plt.savefig('energy_plot.png', transparent = True)
		plt.clf()
		#plotting the first term energy
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("First term of Energy")
		plt.ylabel("optimized First Term of energy")
		plt.xlabel("time")
		plt.plot(timearray, firsttermarray,'--o')
		plt.savefig('firstterm_plot.png', transparent = True)
		plt.clf()
		#plotting the first term energy
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Second term of Energy")
		plt.ylabel("optimized Second Term of energy")
		plt.xlabel("time")
		plt.plot(timearray, secondtermarray,'--o')
		plt.savefig('secondterm_plot.png', transparent = True)
		plt.clf()
		#plotting the first term energy
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Third term of Energy")
		plt.ylabel("optimized Third Term of energy")
		plt.xlabel("time")
		plt.plot(timearray, thirdtermarray,'--o')
		plt.savefig('thirdterm_plot.png', transparent = True)
		plt.clf()
		#plotting the first term energy
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Volume of Cell")
		plt.ylabel("Volume")
		plt.xlabel("time")
		plt.plot(timearray, volumearray,'--o')
		plt.savefig('volume_plot.png', transparent = True)
		plt.clf()
		#plotting the energy plot 
		fig,ax = plt.subplots(1,1)
		ax.grid()
		plt.title("Area Percent strain")
		plt.ylabel("Areastrain%")
		plt.xlabel("faceids")
		plt.plot(deformations[0], deformations[1], 'x')
		plt.savefig('deformation_plot_{0:d}.png'.format(int(growthstep-1)), transparent = True)
		plt.clf()
		#plotting mean deformation
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Mean relative deformation ")
		plt.ylabel("relative deformation")
		plt.xlabel("time")
		plt.errorbar(timearray, meanDeformation, stdDeformation, marker='x')
		plt.savefig('mean_deformation_plot.png', transparent = True)
		plt.clf()
		#plotting z-projection summation - Fourth term
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Fourth term : Z-projection summation ")
		plt.ylabel("zprojection^2")
		plt.xlabel("time")
		plt.plot(timearray, fourthtermarray,'--o')
		plt.axhline(y=cell.getBendingThreshold(), color ='r',linestyle='-')
		plt.savefig('fourthterm_plot.png', transparent = True)
		plt.clf()
		#plotting time required for each relaxation time step
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Time required for each time step in minutes")
		plt.xlabel("relaxation time step")
		plt.ylabel("time required in minutes")
		plt.plot(relaxation_time_length,'--o')
		plt.savefig('relaxation_time_step.png', transparent=True)
		plt.clf()
		#plotting function call required for each relaxation step
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Function Call for each time step")
		plt.xlabel("growth step")
		plt.ylabel("Function call required")
		plt.plot(functionCallCounterArray,'--o')
		plt.savefig('functioncall.png', transparent=True)
		plt.clf()
		#plotting the number of cell division that happend  till now
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Number of Cell Divisions = %d, Time=%d"%(divisionCounter[-1],growthstep))
		plt.xlabel("growth step")
		plt.ylabel("Cumulative number of division")
		plt.plot(divisionCounter,'--o')
		plt.savefig('cumulativeDivisionNumber.png', transparent=True)
		plt.clf()
		#plotting the number of cell division that happend  till now
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Sum of All Wall Length")
		plt.xlabel("growth step")
		plt.ylabel("Sum of wall length")
		plt.plot(sumWallLength,'--o')
		plt.savefig('sumWallLength.png', transparent=True)
		plt.clf()
		#plotting the Helfrich bending energy that happend  till now
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Helfrich Energy omega = %d"%(cell.getOmega()))
		plt.xlabel("growth step")
		plt.ylabel("mean face area")
		plt.plot(meanfacearea,'--o')
		plt.savefig('meanfacearea.png', transparent=True)
		plt.clf()
		#plotting the mean face area
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Mean Face Area")
		plt.xlabel("growth step")
		plt.ylabel("Helfrich Energy")
		plt.plot(helfrichEnergy,'--o')
		plt.savefig('helfrichEnergy.png', transparent=True)
		plt.clf()
		#plotting the number of cell division that happend  in each step
		#subtracting division counter of ith step with i-1th step
		divisionCounterRolled = np.roll(divisionCounter,1)
		divisionCounterRolled[0] = 0
		eachStepCount = divisionCounter - divisionCounterRolled
		fig,ax = plt.subplots(1,1)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
		ax.grid()
		plt.title("Number of Cell Divisions in Each Step, Time = %d"%(growthstep))
		plt.xlabel("growth step")
		plt.ylabel("Number of Division")
		plt.plot(eachStepCount,'--o')
		plt.savefig('divisionNumber.png', transparent=True)
		plt.clf()
		#closing all figures
		plt.close('all')
		########################################################
		return
##########################################################################################
#       Function to Plot the Last Growth Rate for the cells
##########################################################################################
def plotGrowthRateSurface(cell, numOfLayer, name=None, alpha = 0.5, Length=1.0,vmax=0,azim = 0, elev = 90):
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-.7*radius,.7*radius))
	ax.set_ylim((-.7*radius,.7*radius))
	ax.set_zlim((0*radius,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	#fig = plt.figure()
	#ax = fig.gca(projection='3d')
	##iterating through the cell##
	##########################################################
	#Get the Max Of Kappa 
	##########################################################
	if vmax == 0:
		vmax = cell.getKappa() + cell.getGrowthVar()  
		vmin = cell.getKappa() - cell.getGrowthVar()
	else:
		vmin = 0.
	##########################################################
	#### Making the COLOR BAR #########################
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	scalarMap.set_array(np.linspace(vmin,vmax,10))
	#cbar = plt.colorbar(scalarMap,orientation="horizontal",shrink=.5, pad=.0, aspect=10)
	cbar = plt.colorbar(scalarMap,shrink=.5, pad=.0, aspect=10)
	cbar.set_label(r"Growth Rate $\kappa$")
	##########################################################
	##########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		faceid = face.getID()#grabbing face id
		if faceid == 1:
			face  = faces.next()
			continue
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		kappa = face.getLastGrowthRate()
		color = scalarMap.to_rgba(kappa)
		ax.add_collection3d(Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolors = color))
		face = faces.next()
				#if face.getID() == 1: break
	#ax.axis("off")
	ax.view_init(elev=elev, azim=azim)
	if name == None:#plot the figure
		plt.show()
	else:
		plt.savefig(name, transparent = True)
	plt.close()
	return
##########################################################################################
#       Function to Plot Stress on the Surface of the Tissue
##########################################################################################
def plotStressSurface(cell, numOfLayer, step = None, alpha = 0.8, 
	Length=1.0,azim = 0, elev = 7, ids= False, name = None,
	zaxisoffset = 0.3,vmax = None, vmin = None):
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
	fig = plt.figure(frameon=False,figsize=(6,6))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-.7*radius,.7*radius))
	ax.set_ylim((-.7*radius,.7*radius))
	ax.set_zlim((0*radius,(1.-zaxisoffset)*1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if (face.getID() == 1) or checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getStressEigenVector1()
		eigenvec2 = face.getStressEigenVector2()
		#################################################
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#normalised: 
		#print face.getID(), eigenvalue1, eigenvalue2
		eigenmax = max(eigenvalue1, eigenvalue2)
		eigenmin = min(eigenvalue1, eigenvalue2)
		eigenvalueratioarray.append(eigenmax-eigenmin)
		#just the difference :eigenvalueratioarray.append(abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2)))#/max(abs(eigenvalue1),abs(eigenvalue2)))
		#########~~~ EIGEN VEC 1 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec1,0))
		V.append(qd.doublearray_getitem(eigenvec1,1))
		W.append(qd.doublearray_getitem(eigenvec1,2))
		#########~~~ EIGEN VEC 2 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec2,0))
		V.append(qd.doublearray_getitem(eigenvec2,1))
		W.append(qd.doublearray_getitem(eigenvec2,2))
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	if vmax == None:
		maxEigenValueRatio = (max(eigenvalueratioarray))
		minEigenValueRatio = (min(eigenvalueratioarray))
	else:
		maxEigenValueRatio = vmax
		minEigenValueRatio = vmin
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	maxeigenvalue = max(np.abs(eigenvalue))
	eigenvalue = np.array(eigenvalue)/maxeigenvalue
	for i in range(len(X)):
		veclength = eigenvalue[i]*np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',length = veclength,pivot='middle',
			zorder=10, linewidths = 1)
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4, linewidths = 2)#quiver without length for older matplotlib
	########    ########    ########    ########    ########
	#                 Plotting the Cell                    #
	########    ########    ########    ########    ########
	######### Color Map
	jet = cm = plt.get_cmap('plasma') 
	cNorm  = colors.Normalize(vmin=minEigenValueRatio, vmax=maxEigenValueRatio)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		#normalised : 
		eigenmax = max(eigenvalue1, eigenvalue2)
		eigenmin = min(eigenvalue1, eigenvalue2)
		ratio = eigenmax-eigenmin
		#just difference : ratio = abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2))#/max(abs(eigenvalue1),abs(eigenvalue2))
		color = scalarMap.to_rgba(ratio)
		if checkExternalFace(face):
			color = 'w'
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		#ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		if ids:
			ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	scalarMap._A = []
	ax.view_init(azim = azim,elev=elev)
	cbar_ax2 = fig.add_axes([0.9, 0.25, 0.05, 0.5])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax2,format=OOMFormatter(-1, mathText=True))
	clrbar.set_label(r"$\sigma_2-\sigma_1$")
	if name:
		plt.savefig(name+'.png',format= 'png', transparent = True, bbox_inches="tight")
		plt.savefig(name+'.pdf',format = 'pdf', transparent = True, bbox_inches="tight")
	else:
		if step != None:
			plt.title("Time = %d"%step)
			plt.savefig('stress_Plot_%d.png'%step, transparent=True)
		else:
			plt.show()
	#plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return
##########################################################################################
#   Function to Plot radial/orthorad feedbackcorrection
#   on the Surface of the Tissue
##########################################################################################
def plotRadialOrthoradialFeedbackCorrectionSurface(cell,
	numOfLayer, step = None, alpha = 0.8, 
	Length=1.0,azim = 0, elev = 7, ids= False, name = None,zaxisoffset=0.3):
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,6))
	ax1 = fig.add_subplot(121,projection = '3d')
	ax2 = fig.add_subplot(122,projection = '3d')
	for ax in [ax1,ax2]:
		ax.set_xlim((-.7*radius,.7*radius))
		ax.set_ylim((-.7*radius,.7*radius))
		ax.set_zlim((0*radius,(1.-zaxisoffset)*1.4*radius))
		ax.axis('off')
		ax.xaxis.pane.set_edgecolor('black')
		ax.yaxis.pane.set_edgecolor('black')
		ax.xaxis.pane.fill = False
		ax.yaxis.pane.fill = False
		ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1 or checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getRadialVector()
		eigenvec2 = face.getOrthoradialVector()
		#################################################
		eigenvalue1 = face.getRadialFeedbackCorrection()
		eigenvalue2 = face.getOrthoradialFeedbackCorrection()
		#################################################
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2)
		#################################################
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#################################################
		#normalised: 
		#print face.getID(), eigenvalue1, eigenvalue2
		eigenvalueratioarray.append(abs(abs(eigenvalue1)- abs(eigenvalue2))/max(abs(eigenvalue1),abs(eigenvalue2)))
		#just the difference :eigenvalueratioarray.append(abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2)))#/max(abs(eigenvalue1),abs(eigenvalue2)))
		####################################################
		#getting the centroid coordinate 
		centroidvector = [face.getXCentralised(),
						face.getYCentralised(),
						face.getZCentralised()]
		#eigenvectors
		vector1 = [qd.doublearray_getitem(eigenvec1,0),
					qd.doublearray_getitem(eigenvec1,1),
					qd.doublearray_getitem(eigenvec1,2)]
		vector2 = [qd.doublearray_getitem(eigenvec2,0),
					qd.doublearray_getitem(eigenvec2,1),
					qd.doublearray_getitem(eigenvec2,2)]
		#computing mid vectors
		#midvec1 = [(p+q)/2. for p,q in zip(centroidvector,vector1)]
		#midvec2 = [(p+q)/2. for p,q in zip(centroidvector,vector2)]

		#now adding for plotting

		#vec1
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector1[0])
		V.append(vector1[1])
		W.append(vector1[2])

		#vec2
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector2[0])
		V.append(vector2[1])
		W.append(vector2[2])
		####################################################
		face = faces.next()
	###
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	maxeigenvalue = max(np.abs(eigenvalue))
	mineigenvalue = min(np.abs(eigenvalue))
	eigenvalue = np.array(eigenvalue)/maxeigenvalue
	####################################################
	for i in range(0,len(X),2):
		color1 = 'g'
		color2 = 'g'
		veclength1 = (eigenvalue[i])*np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		veclength2 = (eigenvalue[i+1])*np.sqrt((U[i+1])**2+(V[i+1])**2+(W[i+1])**2)
		if eigenvalue[i]<0.:
			color1 = 'r'
		if eigenvalue[i+1]<0.:
			color2 = 'r'
		ax1.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=color1,length = veclength1,
			pivot='middle',zorder=10, linewidths = 1,)
		#        headwidth= 1e-8, headlength=1e-8)
		ax2.quiver(X[i+1], Y[i+1], Z[i+1], U[i+1], V[i+1], W[i+1],color=color2,length = veclength2,
			pivot='middle',zorder=10, linewidths = 1)
		#headwidth=0.01, headlength=0.01)
	#########################################################
	#                 Plotting the Cell                     #
	#########################################################
	mineigenvalue1 = min(eigenvalue1array)
	maxeigenvalue1 = max(eigenvalue1array)
	mineigenvalue2 = min(eigenvalue2array)
	maxeigenvalue2 = max(eigenvalue2array)
	######### Color Maps ###################
	jet1 = cm = plt.get_cmap('magma') 
	cNorm1  = colors.Normalize(vmin=mineigenvalue1, vmax=maxeigenvalue1)
	scalarMapRad = cmx.ScalarMappable(norm=cNorm1, cmap=jet1)
	###################
	jet2 = cm = plt.get_cmap('magma') 
	cNorm2  = colors.Normalize(vmin=mineigenvalue2, vmax=maxeigenvalue2)
	scalarMapOrth = cmx.ScalarMappable(norm=cNorm2, cmap=jet2)
	#########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		#################################################
		eigenvalue1 = face.getRadialFeedbackCorrection()
		eigenvalue2 = face.getOrthoradialFeedbackCorrection()
		#################################################
		color1 = scalarMapRad.to_rgba(eigenvalue1)
		color2 = scalarMapOrth.to_rgba(eigenvalue2)
		if checkExternalFace(face):
			color1 = 'w'
			color2 = 'w'
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color1)
		pc.set_edgecolor('k')
		ax1.add_collection3d(pc)
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color2)
		pc.set_edgecolor('k')
		ax2.add_collection3d(pc)
		#################################################
		#ax1.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#ax2.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#################################################
		if ids:
			ax2.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
			ax1.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		#################################################
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	scalarMapRad._A = []
	scalarMapOrth._A = []
	labels = [r'Radial reorganisation',r'Orthoradial reorganisation']
	
	formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
	for ax,scalarMap,label in zip([ax1,ax2],[scalarMapRad,scalarMapOrth],labels):
		ax.view_init(azim = azim,elev=elev)
		clrbar = plt.colorbar(scalarMap,ax = ax,shrink=0.5,aspect = 10,format=OOMFormatter(-3, mathText=True),
							  pad=-0.1)
		clrbar.set_label(label)
	#################################################
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1,wspace = 0.)
	#################################################
	if name:
		fig.savefig(name+'.png',format = 'png', transparent = True, bbox_inches="tight",dpi = 300)
		fig.savefig(name+'.pdf',format = 'pdf', transparent=True,bbox_inches='tight')
	else:
		if step != None:
			plt.title("Time = %d"%step)
			plt.savefig('stress_Plot_%d.png'%step, transparent=True,bbox_inches='tight')
		else:
			plt.show()
	#################################################
	return

##########################################################################################
#   Function to Plot principal direction feedbackcorrection
#   on the Surface of the Tissue
##########################################################################################

def plotPrincipalDeformationFeedbackCorrectionSurface(cell,
	numOfLayer, step = None, alpha = 0.8, 
	Length=1.0,azim = 0, elev = 7, ids= False, name = None,zaxisoffset=0.3):
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,6))
	ax1 = fig.add_subplot(121,projection = '3d')
	ax2 = fig.add_subplot(122,projection = '3d')
	for ax in [ax1,ax2]:
		ax.set_xlim((-.7*radius,.7*radius))
		ax.set_ylim((-.7*radius,.7*radius))
		ax.set_zlim((0*radius,(1.-zaxisoffset)*1.4*radius))
		ax.axis('off')
		ax.xaxis.pane.set_edgecolor('black')
		ax.yaxis.pane.set_edgecolor('black')
		ax.xaxis.pane.fill = False
		ax.yaxis.pane.fill = False
		ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1 or checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getDeformationEigenVector1()
		eigenvec2 = face.getDeformationEigenVector2()
		#################################################
		eigenvalue1 = face.getPrincipalDeformationDirection1FeedbackCorrection()
		eigenvalue2 = face.getPrincipalDeformationDirection2FeedbackCorrection()
		#################################################
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2)
		#################################################
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#################################################
		#normalised: 
		#print face.getID(), eigenvalue1, eigenvalue2
		eigenvalueratioarray.append(abs(abs(eigenvalue1)- abs(eigenvalue2))/max(abs(eigenvalue1),abs(eigenvalue2)))
		#just the difference :eigenvalueratioarray.append(abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2)))#/max(abs(eigenvalue1),abs(eigenvalue2)))
		####################################################
		#getting the centroid coordinate 
		centroidvector = [face.getXCentralised(),
						face.getYCentralised(),
						face.getZCentralised()]
		#eigenvectors
		vector1 = [qd.doublearray_getitem(eigenvec1,0),
					qd.doublearray_getitem(eigenvec1,1),
					qd.doublearray_getitem(eigenvec1,2)]
		vector2 = [qd.doublearray_getitem(eigenvec2,0),
					qd.doublearray_getitem(eigenvec2,1),
					qd.doublearray_getitem(eigenvec2,2)]
		#computing mid vectors
		#midvec1 = [(p+q)/2. for p,q in zip(centroidvector,vector1)]
		#midvec2 = [(p+q)/2. for p,q in zip(centroidvector,vector2)]

		#now adding for plotting

		#vec1
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector1[0])
		V.append(vector1[1])
		W.append(vector1[2])

		#vec2
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector2[0])
		V.append(vector2[1])
		W.append(vector2[2])
		####################################################
		face = faces.next()
	###
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	maxeigenvalue = max(np.abs(eigenvalue))
	mineigenvalue = min(np.abs(eigenvalue))
	eigenvalue = np.array(eigenvalue)/maxeigenvalue
	####################################################
	for i in range(0,len(X),2):
		color1 = 'g'
		color2 = 'g'
		veclength1 = (eigenvalue[i])*np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		veclength2 = (eigenvalue[i+1])*np.sqrt((U[i+1])**2+(V[i+1])**2+(W[i+1])**2)
		if eigenvalue[i]<0.:
			color1 = 'r'
		if eigenvalue[i+1]<0.:
			color2 = 'r'
		ax1.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=color1,length = veclength1,
			pivot='middle',zorder=10, linewidths = 1,)
		#        headwidth= 1e-8, headlength=1e-8)
		ax2.quiver(X[i+1], Y[i+1], Z[i+1], U[i+1], V[i+1], W[i+1],color=color2,length = veclength2,
			pivot='middle',zorder=10, linewidths = 1)
		#headwidth=0.01, headlength=0.01)
	#########################################################
	#                 Plotting the Cell                     #
	#########################################################
	mineigenvalue1 = min(eigenvalue1array)
	maxeigenvalue1 = max(eigenvalue1array)
	mineigenvalue2 = min(eigenvalue2array)
	maxeigenvalue2 = max(eigenvalue2array)
	######### Color Maps ###################
	jet1 = cm = plt.get_cmap('magma') 
	cNorm1  = colors.Normalize(vmin=mineigenvalue1, vmax=maxeigenvalue1)
	scalarMapRad = cmx.ScalarMappable(norm=cNorm1, cmap=jet1)
	###################
	jet2 = cm = plt.get_cmap('magma') 
	cNorm2  = colors.Normalize(vmin=mineigenvalue2, vmax=maxeigenvalue2)
	scalarMapOrth = cmx.ScalarMappable(norm=cNorm2, cmap=jet2)
	#########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		#################################################
		eigenvalue1 = face.getPrincipalDeformationDirection1FeedbackCorrection()
		eigenvalue2 = face.getPrincipalDeformationDirection2FeedbackCorrection()
		#################################################
		color1 = scalarMapRad.to_rgba(eigenvalue1)
		color2 = scalarMapOrth.to_rgba(eigenvalue2)
		if checkExternalFace(face):
			color1 = 'w'
			color2 = 'w'
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color1)
		pc.set_edgecolor('k')
		ax1.add_collection3d(pc)
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color2)
		pc.set_edgecolor('k')
		ax2.add_collection3d(pc)
		#################################################
		#ax1.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#ax2.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#################################################
		if ids:
			ax2.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
			ax1.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		#################################################
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	scalarMapRad._A = []
	scalarMapOrth._A = []
	labels = [r'Deformation direction 1',r'Deformation direction 2']
	
	formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
	for ax,scalarMap,label in zip([ax1,ax2],[scalarMapRad,scalarMapOrth],labels):
		ax.view_init(azim = azim,elev=elev)
		clrbar = plt.colorbar(scalarMap,ax = ax,shrink=0.5,aspect = 10,format=OOMFormatter(-3, mathText=True),
							  pad=-0.1)
		clrbar.set_label(label)
	#################################################
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1,wspace = 0.)
	#################################################
	if name:
		fig.savefig(name+'.png',format = 'png', transparent = True, bbox_inches="tight",dpi = 300)
		fig.savefig(name+'.pdf',format = 'pdf', transparent=True,bbox_inches='tight')
	else:
		if step != None:
			plt.title("Time = %d"%step)
			plt.savefig('stress_Plot_%d.png'%step, transparent=True,bbox_inches='tight')
		else:
			plt.show()
	#################################################
	return
##########################################################################################
#   Function to Plot principal direction feedbackcorrection
#   on the Surface of the Tissue
##########################################################################################

def plotTracePrincipalDeformationFeedbackCorrectionSurface(cell,
	numOfLayer, step = None, alpha = 0.8, 
	Length=1.0,azim = 0, elev = 7, ids= False, name = None,zaxisoffset=0.3):
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,6))
	ax1 = fig.add_subplot(121,projection = '3d')
	ax2 = fig.add_subplot(122,projection = '3d')
	for ax in [ax1,ax2]:
		ax.set_xlim((-.7*radius,.7*radius))
		ax.set_ylim((-.7*radius,.7*radius))
		ax.set_zlim((0*radius,(1.-zaxisoffset)*1.4*radius))
		ax.axis('off')
		ax.xaxis.pane.set_edgecolor('black')
		ax.yaxis.pane.set_edgecolor('black')
		ax.xaxis.pane.fill = False
		ax.yaxis.pane.fill = False
		ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1 or checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getDeformationEigenVector1()
		eigenvec2 = face.getDeformationEigenVector2()
		#################################################
		eigenvalue1 = face.getPrincipalDeformationDirection1FeedbackCorrection()
		eigenvalue2 = face.getPrincipalDeformationDirection2FeedbackCorrection()
		#################################################
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2)
		#################################################
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#################################################
		#normalised: 
		#print face.getID(), eigenvalue1, eigenvalue2
		eigenvalueratioarray.append(eigenvalue1+eigenvalue2)
		#just the difference :eigenvalueratioarray.append(abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2)))#/max(abs(eigenvalue1),abs(eigenvalue2)))
		####################################################
		#getting the centroid coordinate 
		centroidvector = [face.getXCentralised(),
						face.getYCentralised(),
						face.getZCentralised()]
		#eigenvectors
		vector1 = [qd.doublearray_getitem(eigenvec1,0),
					qd.doublearray_getitem(eigenvec1,1),
					qd.doublearray_getitem(eigenvec1,2)]
		vector2 = [qd.doublearray_getitem(eigenvec2,0),
					qd.doublearray_getitem(eigenvec2,1),
					qd.doublearray_getitem(eigenvec2,2)]
		#computing mid vectors
		#midvec1 = [(p+q)/2. for p,q in zip(centroidvector,vector1)]
		#midvec2 = [(p+q)/2. for p,q in zip(centroidvector,vector2)]

		#now adding for plotting

		#vec1
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector1[0])
		V.append(vector1[1])
		W.append(vector1[2])

		#vec2
		X.append(centroidvector[0])
		Y.append(centroidvector[1])
		Z.append(centroidvector[2])
		#getting the vector headings
		U.append(vector2[0])
		V.append(vector2[1])
		W.append(vector2[2])
		####################################################
		face = faces.next()
	###
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	maxeigenvalue = max(np.abs(eigenvalue))
	mineigenvalue = min(np.abs(eigenvalue))
	eigenvalue = np.array(eigenvalue)/maxeigenvalue
	####################################################
	for i in range(0,len(X),2):
		color1 = 'g'
		color2 = 'g'
		veclength1 = (eigenvalue[i])*np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		veclength2 = (eigenvalue[i+1])*np.sqrt((U[i+1])**2+(V[i+1])**2+(W[i+1])**2)
		if eigenvalue[i]<0.:
			color1 = 'r'
		if eigenvalue[i+1]<0.:
			color2 = 'r'
		ax1.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=color1,length = veclength1,
			pivot='middle',zorder=10, linewidths = 1,)
		#        headwidth= 1e-8, headlength=1e-8)
		#ax2.quiver(X[i+1], Y[i+1], Z[i+1], U[i+1], V[i+1], W[i+1],color=color2,length = veclength2,
		#    pivot='middle',zorder=10, linewidths = 1)
		#headwidth=0.01, headlength=0.01)
	#########################################################
	#                 Plotting the Cell                     #
	#########################################################
	mineigenvalueratio = min(eigenvalueratioarray)
	maxeigenvalueratio = max(eigenvalueratioarray)
	#mineigenvalue2 = min(eigenvalue2array)
	#maxeigenvalue2 = max(eigenvalue2array)
	######### Color Maps ###################
	jet1 = cm = plt.get_cmap('magma') 
	cNorm1  = colors.Normalize(vmin=mineigenvalueratio, vmax=maxeigenvalueratio)
	scalarMapRad = cmx.ScalarMappable(norm=cNorm1, cmap=jet1)
	###################
	#jet2 = cm = plt.get_cmap('magma') 
	#cNorm2  = colors.Normalize(vmin=mineigenvalue2, vmax=maxeigenvalue2)
	#scalarMapOrth = cmx.ScalarMappable(norm=cNorm2, cmap=jet2)
	#########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		#################################################
		eigenvalue1 = face.getPrincipalDeformationDirection1FeedbackCorrection()
		eigenvalue2 = face.getPrincipalDeformationDirection2FeedbackCorrection()
		#################################################
		color1 = scalarMapRad.to_rgba(eigenvalue1+eigenvalue2)
		#print eigenvalue2+eigenvalue1, face.getRadialFeedbackCorrection()+face.getOrthoradialFeedbackCorrection()
		color2 = scalarMapRad.to_rgba((face.getRadialFeedbackCorrection()+face.getOrthoradialFeedbackCorrection()))
		if checkExternalFace(face):
			color1 = 'w'
			color2 = 'w'
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color1)
		pc.set_edgecolor('k')
		ax1.add_collection3d(pc)
		#################################################
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color2)
		pc.set_edgecolor('k')
		ax2.add_collection3d(pc)
		#################################################
		#ax1.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#ax2.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='m')
		#################################################
		if ids:
			ax2.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
			ax1.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		#################################################
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	scalarMapRad._A = []
	#scalarMapOrth._A = []
	labels = [r'Tr(Deform feedback)',r'Tr(Rad/Orthorad Feedback)']
	formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
	for ax,scalarMap,label in zip([ax1,ax2],[scalarMapRad,scalarMapRad],labels):
		ax.view_init(azim = azim,elev=elev)
		clrbar = plt.colorbar(scalarMap,ax = ax,shrink=0.5,aspect = 10,format=OOMFormatter(-3, mathText=True),
							  pad=-0.1)
		clrbar.set_label(label)
	#################################################
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1,wspace = 0.)
	#################################################
	if name:
		fig.savefig(name+'.png',format = 'png', transparent = True, bbox_inches="tight",dpi = 300)
		fig.savefig(name+'.pdf',format = 'pdf', transparent=True,bbox_inches='tight')
	else:
		if step != None:
			plt.title("Time = %d"%step)
			plt.savefig('stress_Plot_%d.png'%step, transparent=True,bbox_inches='tight')
		else:
			plt.show()
	#################################################
	return
##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotPrimaryStrainSurface(cell, numOfLayer,targetface = 0,step = None, alpha = 0.8, Length=1.0):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	#import matplotlib
	#matplotlib.use('agg')
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	################################################################################
	###########################################
	### Making the directory to save the figure
	###########################################
	import os
	directory = 'strain-plots'
	if not os.path.exists(directory):
		os.makedirs(directory)
	################################################################################
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-.5*radius,.5*radius))
	ax.set_ylim((-.5*radius,.5*radius))
	ax.set_zlim((-0.,1.*radius))
	ax.set_xlabel("x-axis")
	ax.set_ylabel("y-axis")
	ax.set_zlabel("z-axis")
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	#ax.xaxis.pane.fill = False
	#ax.yaxis.pane.fill = False
	#ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Strain vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face = faces.next()
			continue
		"""if face.getID() <= (perimeterHexNum+1):
			face  = faces.next()
			continue"""
		eigenvec1 = face.getStrainEigenVector1()
		eigenvec2 = face.getStrainEigenVector2()
		eigenvalue1 = face.getStrainEigenValue1()
		eigenvalue2 = face.getStrainEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#########~~~ EIGEN VEC 1 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec1,0))
		V.append(qd.doublearray_getitem(eigenvec1,1))
		W.append(qd.doublearray_getitem(eigenvec1,2))
		 #########~~~ EIGEN VEC 2 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec2,0))
		V.append(qd.doublearray_getitem(eigenvec2,1))
		W.append(qd.doublearray_getitem(eigenvec2,2))
		#ax.scatter(X[-1],Y[-1],Z[-1])
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	eigenvalue = np.array(eigenvalue)
	mineigenvalue = np.min(np.abs(eigenvalue))
	maxeigenvalue = np.max(np.abs(eigenvalue))
	print "------------------------------------------"
	print "---     Max & Min Eigenvalue           ---"
	print "------------------------------------------"
	print "Max Strain eigen Value : ", maxeigenvalue
	print "Min Strain eigen Value : ", mineigenvalue
	print "------------------------------------------"
	print "---     Rescaling Eigen value          ---"
	print "------------------------------------------"
	rescalevalue = maxeigenvalue
	#rescalevalue = 0.0021906004466#MAX strain
	eigenvalue= eigenvalue/rescalevalue
	print "Max Strain eigen Value : ", np.max(eigenvalue)
	print "Min Strain eigen Value : ", np.min(eigenvalue)
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	for i in range(len(X)):
		veclength = eigenvalue[i]*np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		#print eigenvalue[i]
		if eigenvalue[i]>=0:
			ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',length = veclength,pivot='tail',zorder=4)
		else:
			ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',length = veclength,zorder=4)
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
	########    ########    ########    ########    ########
	#                 Plotting the Cell                    #
	########    ########    ########    ########    ########
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=0, vmax=1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		if False: #face.getID() <= (perimeterHexNum+1):
			color='w'
		else:
			#print face.getZCentralised()
			eigenvalue1 = face.getStrainEigenValue1()
			eigenvalue2 = face.getStrainEigenValue2()
			ratio = abs(abs(eigenvalue1)-abs(eigenvalue2))/(max(abs(eigenvalue1),abs(eigenvalue2)))
			#print "face ID : ", face.getID(), " ratio : ", ratio
			if eigenvalue1 < 0 or eigenvalue2 < 0:
				#color = scalarMap.to_rgba(ratio)
				color = 'r'
			else:
				color = 'w'
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		if faceid == targetface:
			ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	scalarMap._A = []
	#clrbar = plt.colorbar(scalarMap)
	#clrbar.set_label("Magnitude of Strain Anisotrophy")
	ax.view_init(elev=90,azim=0)
	if  step != None:
			plt.title("Time = %d"%step)
			plt.savefig(directory+r"/"+'strain_Plot_%d.png'%step, transparent=True)
	else :
		plt.show()
	#print xcenarray-0.5
	#plt.close("all")
	return


##########################################################################################
#       Function to Plot Magnitude of Normal Forces on the Faces of the Cell
##########################################################################################
def plotNormalForce(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0):
	###########################################
	### Making the directory to save the figure
	###########################################
	import os
	forceDirectory = 'normalForcePlot'
	if not os.path.exists(forceDirectory):
		os.makedirs(forceDirectory)
	###########################################
	### FORMATING FOR COLORBAR SCALE 
	###########################################
	def fmt(x, pos):
		a, b = '{:.2e}'.format(x).split('e')
		b = int(b)
		return r"${} \times 10^{{{}}}$".format(a, b)
	###########################################
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	import matplotlib.ticker as ticker
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	magnitudeNormalForce =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() ==1:
			face  = faces.next()
			continue
		#if face.getID() <= (perimeterHexNum+1):
		#    face  = faces.next()
		#    continue
		#Getting the Normal Force
		normalForce = face.getNormalForce()
		#################################################
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(normalForce,0))
		V.append(qd.doublearray_getitem(normalForce,1))
		W.append(qd.doublearray_getitem(normalForce,2))
		#################################################
		magnitudeNormalForce.append(np.sqrt(U[-1]**2+V[-1]**2+W[-1]**2))
		#################################################
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxMagnitude = (max(magnitudeNormalForce))
	minMagnitude = (min(magnitudeNormalForce))
	#print "Min & Max Magnitude", minMagnitude, maxMagnitude
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	for i in range(len(X)):
		veclength = magnitudeNormalForce[i]
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',length = veclength,pivot='tail',zorder=4, linewidths = 2)
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4, linewidths = 2)#quiver without length for older matplotlib
	###############################################################
	#                 Plotting the Cell                           #
	###############################################################
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	faceCounter = 0
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() ==1:
			face  = faces.next()
			continue
		#if face.getID() <= (perimeterHexNum+1):
		#    face  = faces.next()
		#    continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		########################################################################################
		color = scalarMap.to_rgba(magnitudeNormalForce[faceCounter])
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		ax.add_collection3d(Poly3DCollection(verts,alpha = alpha,facecolors = color,linewidths=1,zorder=0))
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
		faceCounter+= 1
		#if face.getID() == 1: break
	#plt.clf()
	scalarMap._A = []
	clrbar = plt.colorbar(scalarMap, format=ticker.FuncFormatter(fmt))
	clrbar.set_label("Magnitude of Normal Force")
	if step != None:
		plt.title("Time = %d"%step)
		plt.savefig(forceDirectory+r"/"+'Normal_Force_Plot_%d.png'%step, transparent=True)
	else:
		plt.show()
	plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return

##########################################################################################
#       Function to Plot STRAIN Magnitude (Trace of Strain) on the Surface of the Tissue
##########################################################################################
def plotStrainMagnitude(cell, numOfLayer,step = None,targetface =10, alpha = 0.8, Length=1.0, save=False,azim = -70, elev=50):
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
	 fig = plt.figure(frameon=False,figsize=(10,8))
	 #fig = plt.figure(frameon=False)
	 fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	 ax = Axes3D(fig)
	 ax.set_xlim((-0.7*radius,0.7*radius))
	 ax.set_ylim((-0.7*radius,0.7*radius))
	 ax.set_zlim((-0.,1.4*radius))
	 ax.axis('off')
	 ax.xaxis.pane.set_edgecolor('black')
	 ax.yaxis.pane.set_edgecolor('black')
	 #ax.xaxis.pane.fill = False
	 #ax.yaxis.pane.fill = False
	 #ax.zaxis.pane.fill = False
	 ########################################################################
	 #                    Plotting the Strain vectors                       #
	 ########################################################################
	 X = []
	 Y = []
	 Z = []
	 U = []
	 V = []
	 W = []
	 eigenvalue =[]
	 eigenvalue1array = []
	 eigenvalue2array = []
	 eigenvalueratioarray = []
	 tracearray = []
	 faces = qd.CellFaceIterator(cell)
	 face = faces.next()
	 while face != None:
		  if face.getID() == 1:
				face  = faces.next()
				continue
		  if checkExternalFace(face):
				face  = faces.next()
				continue
		  eigenvec1 = face.getStrainEigenVector1()
		  eigenvec2 = face.getStrainEigenVector2()
		  eigenvalue1 = face.getStrainEigenValue1()
		  eigenvalue2 = face.getStrainEigenValue2()
		  eigenvalue.append(eigenvalue1)
		  eigenvalue.append(eigenvalue2) 
		  eigenvalue1array.append(eigenvalue1)
		  eigenvalue2array.append(eigenvalue2)
		  #straindeterminantarray.append(face.getStrainDeterminant())
		  tracearray.append(face.getStrainTrace())
		  #print face.getID(), straindeterminantarray[-1]
		  #########~~~ EIGEN VEC 1 ~~~#########
		  #getting the centralised coordinate of centroid
		  X.append(face.getXCentralised())
		  Y.append(face.getYCentralised())
		  Z.append(face.getZCentralised())
		  #getting the vector headings
		  U.append(qd.doublearray_getitem(eigenvec1,0))
		  V.append(qd.doublearray_getitem(eigenvec1,1))
		  W.append(qd.doublearray_getitem(eigenvec1,2))
			#########~~~ EIGEN VEC 2 ~~~#########
		  #getting the centralised coordinate of centroid
		  X.append(face.getXCentralised())
		  Y.append(face.getYCentralised())
		  Z.append(face.getZCentralised())
		  #getting the vector headings
		  U.append(qd.doublearray_getitem(eigenvec2,0))
		  V.append(qd.doublearray_getitem(eigenvec2,1))
		  W.append(qd.doublearray_getitem(eigenvec2,2))
		  #ax.scatter(X[-1],Y[-1],Z[-1])
		  face = faces.next()
	 ###getting Maximum Eigenvalue ratio
	 maxTrace = (max(tracearray))
	 minTrace = (min(tracearray))
	 #print "Max Strain magnitude :", maxTrace
	 #print "Min Strain magnitude :", minTrace
	 #print " rescalling by max Value"
	 maxEigenValue = max(map(abs,eigenvalue))
	 for i in range(len(X)):
		  #print " veclength : ", veclength, (eigenvalue[i]/maxEigenValue)
		  veclength = (eigenvalue[i]/maxEigenValue)
		  if veclength <= 0:
				colorvec = 'r'
		  else:
				colorvec = 'k'
		  ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=colorvec,length = veclength,pivot='tail',zorder=4)
		  #ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
	 ########    ########    ########    ########    ########
	 #                 Plotting the Cell                    #
	 ########    ########    ########    ########    ########
	 ######### Color Map
	 jet = cm = plt.get_cmap('plasma') 
	 maxvalue = maxTrace
	 minvalue = minTrace
	 #print "Max value", maxvalue, " minvalue", minvalue
	 cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	 scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	 faces = qd.CellFaceIterator(cell)
	 ###################
	 face = faces.next()
	 xcenarray = []
	 ycenarray = []
	 zcenarray = []
	 while (face != None):
		  if face.getID() == 1:
				face  = faces.next()
				continue
		  faceid = face.getID()#grabbing face id
		  xlist = []
		  ylist = []
		  zlist = []
		  xproj = []
		  yproj = []
		  zproj = []
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
		  verts = [zip(xlist, ylist,zlist)]
		  #adding to 3d plot
		  xcenarray.append(face.getXCentralised())
		  ycenarray.append(face.getYCentralised())
		  zcenarray.append(face.getZCentralised())
		  #print "face ID : ", face.getID(), " ratio : ", ratio
		  #ratio = (face.getStrainTrace()+minTrace)/(minTrace+maxTrace)
		  ratio = face.getStrainTrace()
		  if checkExternalFace(face):
				ratio = 0.
		  color = scalarMap.to_rgba(ratio)
		  #print face.getZCentralised(), alpha_fac
		  #ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		  pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		  pc.set_edgecolor('k')
		  ax.add_collection3d(pc)
		  ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
		  if faceid == targetface:
				ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1],faceid,
						 fontsize=20,color='black')
		  face = faces.next()
		  #if face.getID() == 1: break
	 #plt.clf()
	 ax.view_init(azim = azim, elev = elev)
	 scalarMap._A = []
	 clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	 clrbar.set_label("Rescaled Magnitude of Strain", fontsize = 20)
	 ax.set_title("Time %d : Magnitude of Strain : Max %.4f   ; Min %.4f"%(step,maxTrace,minTrace), fontsize = 20)
	 #print xcenarray-0.5
	 #plt.close("all")
	 if save:
		  saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressMagnitude"
		  import os
		  if not os.path.exists(saveDirectory):
				os.makedirs(saveDirectory)
		  plt.savefig(saveDirectory+r"/strainMagnitude_Time=%03d.png"%(step),transparent=True)
		  plt.close()
	 return
##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotStrainDifferenceSurface(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False,azim = -70, elev=50, norm=True):
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
	fig = plt.figure(frameon=False,figsize=(10,8))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	#ax.xaxis.pane.fill = False
	#ax.yaxis.pane.fill = False
	#ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Strain vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		if checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getStrainEigenVector1()
		eigenvec2 = face.getStrainEigenVector2()
		eigenvalue1 = face.getStrainEigenValue1()
		eigenvalue2 = face.getStrainEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = abs(eigenvalue1)
		eigen2 = abs(eigenvalue2)
		ratio = abs(max(eigen1,eigen2)- min(eigen1,eigen2))/max(abs(eigenvalue1),abs(eigenvalue2))
		ratio = abs(max(eigen1,eigen2)- min(eigen1,eigen2))#/max(abs(eigenvalue1),abs(eigenvalue2))
		eigenvalueratioarray.append(ratio)
		#########~~~ EIGEN VEC 1 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec1,0))
		V.append(qd.doublearray_getitem(eigenvec1,1))
		W.append(qd.doublearray_getitem(eigenvec1,2))
		 #########~~~ EIGEN VEC 2 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec2,0))
		V.append(qd.doublearray_getitem(eigenvec2,1))
		W.append(qd.doublearray_getitem(eigenvec2,2))
		#ax.scatter(X[-1],Y[-1],Z[-1])
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxEigenValueRatio = (max(eigenvalueratioarray))
	minEigenValueRatio = (min(eigenvalueratioarray))
	maxEigenValue = max(map(abs,eigenvalue))
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	########    ########    ########    ########    ########
	#                 Plotting the Cell                    #
	########    ########    ########    ########    ########
	######### Color Map
	jet = cm = plt.get_cmap('plasma')
	if norm:
		maxvalue = 1.#maxEigenValueRatio
		minvalue = minEigenValueRatio
		normMax = maxEigenValueRatio # value to normalize by
	else:
		maxvalue = maxEigenValueRatio
		minvalue = minEigenValueRatio
		normMax = 1.#maxEigenValueRatio # value to normalize by
	######################################################
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
	  if face.getID() == 1:
			face  = faces.next()
			continue
	  faceid = face.getID()#grabbing face id
	  xlist = []
	  ylist = []
	  zlist = []
	  xproj = []
	  yproj = []
	  zproj = []
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
	  verts = [zip(xlist, ylist,zlist)]
	  #adding to 3d plot
	  xcenarray.append(face.getXCentralised())
	  ycenarray.append(face.getYCentralised())
	  zcenarray.append(face.getZCentralised())
	  #print face.getZCentralised()
	  eigenvalue1 = face.getStrainEigenValue1()
	  eigenvalue2 = face.getStrainEigenValue2()
	  eigen1 = abs(eigenvalue1)
	  eigen2 = abs(eigenvalue2)
	  ratio = abs(max(eigen1,eigen2)- min(eigen1,eigen2))/max(abs(eigenvalue1),abs(eigenvalue2))
	  ratio = abs(max(eigen1,eigen2)- min(eigen1,eigen2))/normMax#/max(abs(eigenvalue1),abs(eigenvalue2))
	  if checkExternalFace(face):
			ratio = 0.
	  #print "face ID : ", face.getID(), " ratio : ", ratio
	  color = scalarMap.to_rgba(ratio)
	  #print face.getID(), ratio
	  #print face.getZCentralised(), alpha_fac
	  #ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
	  pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
	  pc.set_edgecolor('k')
	  ax.add_collection3d(pc)
	  ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
	  #ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1],face.getID(),fontsize = 18)        
	  face = faces.next()
	  #if face.getID() == 1: break
	#plt.clf()
	for i in range(len(X)):
	  veclength = np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
	  #print " veclength : ", veclength, (eigenvalue[i]/maxEigenValue)
	  veclength *= (eigenvalue[i]/maxEigenValue)
	  if veclength <= 0:
			colorvec = 'r'
	  else:
			colorvec = 'k'
	  ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=colorvec,length = veclength,pivot='tail',zorder=-1)
	  #ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.873, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label("Magnitude of Stress Anisotrophy", fontsize = 18)
	clrbar.ax.tick_params(labelsize=14)
	#ax.set_title("Time %d : Anisotrophy of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	ax.view_init(azim = azim, elev= elev)
	if save:
	  saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressSurface"
	  import os
	  if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
	  plt.savefig(saveDirectory+r"/strainDifference_Time=%03d.png"%(step),transparent=True)
	  plt.close()
	return
##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotStrainSurface(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False,azim = -70, elev=50):
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
	fig = plt.figure(frameon=False,figsize=(10,8))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	#ax.xaxis.pane.fill = False
	#ax.yaxis.pane.fill = False
	#ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Strain vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		if checkExternalFace(face):
			face  = faces.next()
			continue
		eigenvec1 = face.getStrainEigenVector1()
		eigenvec2 = face.getStrainEigenVector2()
		eigenvalue1 = face.getStrainEigenValue1()
		eigenvalue2 = face.getStrainEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = abs(eigenvalue1)
		eigen2 = abs(eigenvalue2)
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		eigenvalueratioarray.append(ratio)
		#########~~~ EIGEN VEC 1 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec1,0))
		V.append(qd.doublearray_getitem(eigenvec1,1))
		W.append(qd.doublearray_getitem(eigenvec1,2))
		 #########~~~ EIGEN VEC 2 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec2,0))
		V.append(qd.doublearray_getitem(eigenvec2,1))
		W.append(qd.doublearray_getitem(eigenvec2,2))
		#ax.scatter(X[-1],Y[-1],Z[-1])
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxEigenValueRatio = (max(eigenvalueratioarray))
	minEigenValueRatio = (min(eigenvalueratioarray))
	maxEigenValue = max(map(abs,eigenvalue))
	#print "Max Eigen Value Ration :", maxEigenValueRatio
	#print "Min Eigen Value Ratio :", minEigenValueRatio
	########    ########    ########    ########    ########
	#                 Plotting the Cell                    #
	########    ########    ########    ########    ########
	######### Color Map
	jet = cm = plt.get_cmap('plasma') 
	maxvalue = 1#maxEigenValueRatio
	minvalue = 0#minEigenValueRatio
	normMax = 2 # value to normalize by
	cNorm  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		eigenvalue1 = face.getStrainEigenValue1()
		eigenvalue2 = face.getStrainEigenValue2()
		eigen1 = abs(eigenvalue1)
		eigen2 = abs(eigenvalue2)
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(max(eigen1,eigen2))
		ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		if checkExternalFace(face):
			ratio = 0.
		#print "face ID : ", face.getID(), " ratio : ", ratio
		color = scalarMap.to_rgba(ratio)
		#print face.getID(), ratio
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1],face.getID(),fontsize = 18)        
		face = faces.next()
		#if face.getID() == 1: break
	#plt.clf()
	for i in range(len(X)):
		veclength = np.sqrt((U[i])**2+(V[i])**2+(W[i])**2)
		#print " veclength : ", veclength, (eigenvalue[i]/maxEigenValue)
		veclength *= (eigenvalue[i]/maxEigenValue)
		if veclength <= 0:
			colorvec = 'r'
		else:
			colorvec = 'k'
		ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=colorvec,length = veclength,pivot='tail',zorder=-1)
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.873, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label("Magnitude of Stress Anisotrophy", fontsize = 18)
	clrbar.ax.tick_params(labelsize=14)
	#ax.set_title("Time %d : Anisotrophy of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	ax.view_init(azim = azim, elev= elev)
	if save:
		saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressSurface"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		plt.savefig(saveDirectory+r"/strainSurface_Time=%03d.png"%(step),transparent=True)
		plt.close()
	return


##########################################################################################
#       Function to Plot Magnitude of Normal Forces on the Faces of the Cell
##########################################################################################
def plotFaceArea(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0,save = False,azim = -70, elev=50):
	 ###########################################
	 ### Making the directory to save the figure
	 ###########################################   
	 import os
	 directory = 'faceAreaPlot'
	 if not os.path.exists(directory):
		  os.makedirs(directory)
	 ###########################################
	 ### FORMATING FOR COLORBAR SCALE 
	 ###########################################
	 def fmt(x, pos):
		  a, b = '{:.2e}'.format(x).split('e')
		  b = int(b)
		  return r'${} \times 10^{{{}}}$'.format(a, b)
	 ###########################################
	 #calculating forces, stress-matrix and strain-matrix
	 #cell.calculateVertexForce()
	 cell.calculateStrain()
	 import matplotlib.colors as colors
	 import matplotlib.cm as cmx
	 import matplotlib.ticker as ticker
	 #import the libraries
	 from mpl_toolkits.mplot3d import Axes3D
	 import matplotlib as mpl
	 from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	 import numpy as np
	 import matplotlib.pyplot as plt
	 #limits of the plot
	 radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	 #plotting part
	 fig = plt.figure(frameon=False,figsize=(10,8))
	 fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	 ax = Axes3D(fig)
	 ax.set_xlim((-0.7*radius,0.7*radius))
	 ax.set_ylim((-0.7*radius,0.7*radius))
	 ax.set_zlim((-0.,1.4*radius))
	 ax.axis('off')
	 ax.xaxis.pane.set_edgecolor('black')
	 ax.yaxis.pane.set_edgecolor('black')
	 ax.xaxis.pane.fill = False
	 ax.yaxis.pane.fill = False
	 ax.zaxis.pane.fill = False
	 ########################################################################
	 #                    Plotting the Stress vectors                       #
	 ########################################################################
	 X = []
	 Y = []
	 Z = []
	 U = []
	 V = []
	 W = []
	 areaArray =[]
	 faces = qd.CellFaceIterator(cell)
	 face = faces.next()
	 while face != None:
		  if face.getID() ==1 :
				face  = faces.next()
				continue
		  areaArray.append(face.getAreaOfFace())
		  #################################################
		  face = faces.next()
	 ###getting Maximum Eigenvalue ratio
	 maxMagnitude = (max(areaArray))
	 minMagnitude = (min(areaArray))
	 #print "Min & Max Magnitude", minMagnitude, maxMagnitude
	 ###############################################################
	 #                 Plotting the Cell                           #
	 ###############################################################
	 ######### Color Map
	 jet = cm = plt.get_cmap('viridis') 
	 cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	 scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	 faces = qd.CellFaceIterator(cell)
	 ###################
	 face = faces.next()
	 faceCounter = 0
	 xcenarray = []
	 ycenarray = []
	 zcenarray = []
	 while (face != None):
		  if face.getID()==1:
				face  = faces.next()
				continue
		  faceid = face.getID()#grabbing face id
		  xlist = []
		  ylist = []
		  zlist = []
		  xproj = []
		  yproj = []
		  zproj = []
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
		  verts = [zip(xlist, ylist,zlist)]
		  #adding to 3d plot
		  xcenarray.append(face.getXCentralised())
		  ycenarray.append(face.getYCentralised())
		  zcenarray.append(face.getZCentralised())
		  ########################################################################################
		  color = scalarMap.to_rgba(face.getAreaOfFace())
		  #print face.getZCentralised(), alpha_fac
		  #ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		  pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		  pc.set_edgecolor('k')
		  ax.add_collection3d(pc)
		  ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		  #ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		  face = faces.next()
		  faceCounter+= 1
		  #if face.getID() == 1: break
	 #plt.clf()
	 #ax.view_init(elev=90., azim=0.)#viewing angle from top
	 scalarMap._A = []
	 #clrbar = plt.colorbar(scalarMap, format=ticker.FuncFormatter(fmt))
	 #clrbar.set_label("Area of Cell")
	 clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	 clrbar.set_label("Area of Cell")
	 ax.set_title("Time %d : Area of Cell : Max %.4f Min %.4f"%(step,maxMagnitude,minMagnitude), fontsize = 20)
	 ax.view_init(azim = azim, elev= elev)
	 if save:
		  saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/faceArea"
		  import os
		  if not os.path.exists(saveDirectory):
				os.makedirs(saveDirectory)
		  plt.savefig(saveDirectory+r"/faceArea_Time=%03d.png"%(step),transparent=True)
		  plt.close()
	 #plt.close("all")
	 #return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	 return
##########################################################################################
#       Function to Plot stiffness of the cells
##########################################################################################
def plotStiffness(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0,azim = 0, elev = 90):
	###########################################
	### Making the directory to save the figure
	###########################################
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	#cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	import matplotlib.ticker as ticker
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	areaArray =[]
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() ==1 :
			face  = faces.next()
			continue
		thisAlpha = face.getAlpha()
		if thisAlpha == 0.:
			thisAlpha = cell.getAlpha()
		areaArray.append(thisAlpha)
		#################################################
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxMagnitude = (max(areaArray))
	minMagnitude = (min(areaArray))
	#print "Min & Max Magnitude", minMagnitude, maxMagnitude
	###############################################################
	#                 Plotting the Cell                           #
	###############################################################
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	faceCounter = 0
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID()==1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		########################################################################################
		facealpha = face.getAlpha()
		if facealpha == 0:
			facealpha = cell.getAlpha()
		color = scalarMap.to_rgba(facealpha)
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
		faceCounter+= 1
		#if face.getID() == 1: break
	#plt.clf()
	ax.view_init(elev=elev, azim=azim)#viewing angle from top
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.877, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	#clrbar = plt.colorbar(scalarMap)
	clrbar.set_label("Stiffness of Cell")
	if step != None:
		#plt.suptitle("Time = %d"%step)
		plt.savefig('stiffness_Plot_%03d.png'%step, transparent=True)
	else:
		plt.show()
	plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return

##########################################################################################
#       Function to Plot Eta of the cells
##########################################################################################
def plotEta(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0,azim = 0, elev = 7):
	###########################################
	### Making the directory to save the figure
	###########################################
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	#cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	import matplotlib.ticker as ticker
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	areaArray =[]
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() ==1 :
			face  = faces.next()
			continue
		areaArray.append(face.getEta())
		#################################################
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxMagnitude = (max(areaArray))
	minMagnitude = (min(areaArray))
	#print "Min & Max Magnitude", minMagnitude, maxMagnitude
	###############################################################
	#                 Plotting the Cell                           #
	###############################################################
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	faceCounter = 0
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID()==1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		########################################################################################
		facealpha = face.getEta()
		if facealpha == 0:
			facealpha = cell.getEta()
		color = scalarMap.to_rgba(facealpha)
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
		faceCounter+= 1
		#if face.getID() == 1: break
	#plt.clf()
	ax.view_init(elev=elev, azim=azim)#viewing angle from top
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.877, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	ax.view_init(azim = -60, elev = 60)
	#clrbar = plt.colorbar(scalarMap)
	clrbar.set_label("Eta of Tissue, Cell eta = %d"%cell.getEta())
	if step != None:
		#plt.suptitle("Time = %d"%step)
		plt.savefig('eta_Plot_%03d.png'%step, transparent=True)
	else:
		plt.show()
	plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return



##########################################################################################
#       Function to Plot Magnitude of Normal Forces on the Faces of the Cell
##########################################################################################
def plotForceAreaRelation(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0):
	###########################################
	### Making the directory to save the figure
	###########################################   
	import os
	directory = 'faceAreaRelationPlot'
	if not os.path.exists(directory):
		os.makedirs(directory)
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
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	perimeterHexNum = (numOfLayer!=1)*(6*(numOfLayer-1)-1) + 1 # number of perimeter hexgons
	print perimeterHexNum
	#plotting part
	U = []
	V = []
	W = []
	magnitudeNormalForce =[]
	areaofface=[]
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() <= (perimeterHexNum+1):
			face  = faces.next()
			continue
		#Getting the Normal Force
		normalForce = face.getNormalForce()
		U.append(qd.doublearray_getitem(normalForce,0))
		V.append(qd.doublearray_getitem(normalForce,1))
		W.append(qd.doublearray_getitem(normalForce,2))
		#################################################
		magnitudeNormalForce.append(np.sqrt(U[-1]**2+V[-1]**2+W[-1]**2))
		areaofface.append(face.getAreaOfFace())
		#################################################
		face = faces.next()
	########################################################
	##Fiting Function 
	func = np.polyfit(np.array(areaofface), np.array(magnitudeNormalForce), 1)
	plt.plot(areaofface,magnitudeNormalForce,'x')
	fitrange = np.array([min(areaofface),max(areaofface)])
	plt.plot(fitrange, func[1]+func[0]*fitrange,'r')
	plt.xlabel('Area of Face')
	plt.ylabel('magnitude of Normal Force')
	plt.annotate('y=%.4fx+%.4f'%(func[0],func[1]),xy=(0.7,0.2),xycoords='axes fraction',xytext=(0.2,0.8),textcoords='axes fraction')
	if step != None:
		plt.title("Time = %d"%step)
		plt.savefig(directory+r'/'+'Area_Force_Relation_Plot_%d.png'%step, transparent=True)
	else:
		plt.show()
	#plt.close("all")
	return


##########################################################################################
#       Function to Plot Area Distribution of Cells
##########################################################################################
def plotAreaDistribution(cell, numOfLayer, variancearray = [],step = None, alpha = 0.8, Length=1.0):
	###########################################
	### Making the directory to save the figure
	###########################################   
	import os
	directory = 'areaDistributionPlot'
	if not os.path.exists(directory):
		os.makedirs(directory)
	#calculating forces, stress-matrix and strain-matrix
	import matplotlib as mpl
	import numpy as np
	import matplotlib.pyplot as plt
	areaofface=[]
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		areaofface.append(face.getAreaOfFace())
		#################################################
		face = faces.next()
	########################################################
	## Calculating variance 
	areaofface = np.array(areaofface)
	variancearray.append(np.var(areaofface))
	########################################################
	## Plotting Figure
	########################################################
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(1,1,1)
	ax1.hist(areaofface,bins=np.arange(2.2,2.75,0.01),facecolor='green')
	ax1.set_xlabel('Area of Face')
	ax1.set_ylabel('Frequency')
	if step != None:
		ax1.set_title("Time = %d"%step)
		plt.savefig(directory+r'/'+'area_distribution_%d.png'%step, transparent=True)
	else:
		plt.show()
	########################################################
	fig2 = plt.figure(2)
	ax2 = fig2.add_subplot(1,1,1)
	ax2.plot(variancearray,'-x')
	ax2.set_title("Variance of Area")
	ax2.set_xlabel('Time step')
	ax2.set_ylabel('Variance of Area')
	if step != None:
		plt.savefig(directory+r'/'+'variance_of_area.png', transparent=True)
	else:
		plt.show()
	#plt.close("all")
	return variancearray

##########################################################################################
#       Function to Plot Stiffness of the cells
##########################################################################################
def plotFaceStiffness(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0):
	###########################################
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	cell.calculateStressStrain()
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	import matplotlib.ticker as ticker
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.7*radius,0.7*radius))
	ax.set_ylim((-0.7*radius,0.7*radius))
	ax.set_zlim((-0.,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	stiffnessArray =[]
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		facestiffness = face.getAlpha()
		if facestiffness == 0:
			facestiffness = cell.getAlpha()
		stiffnessArray.append(facestiffness)
		#################################################
		face = faces.next()
	###getting Maximum Eigenvalue ratio
	maxMagnitude = (max(stiffnessArray))
	minMagnitude = (min(stiffnessArray))
	#print "Min & Max Magnitude", minMagnitude, maxMagnitude
	###############################################################
	#                 Plotting the Cell                           #
	###############################################################
	######### Color Map
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=minMagnitude, vmax=maxMagnitude)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	faces = qd.CellFaceIterator(cell)
	###################
	face = faces.next()
	faceCounter = 0
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
		#print "== Face ID : ", faceid, "=="
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		########################################################################################
		facestiffness = face.getAlpha()
		if facestiffness == 0:
			facestiffness = cell.getAlpha()
		color = scalarMap.to_rgba(facestiffness)
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		ax.add_collection3d(Poly3DCollection(verts,alpha = alpha,facecolors = color,linewidths=1,zorder=0))
		#ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
		faceCounter+= 1
		#if face.getID() == 1: break
	#plt.clf()
	scalarMap._A = []
	clrbar = plt.colorbar(scalarMap)
	clrbar.set_label("Stiffness of Cell")
	if step != None:
		plt.title("Time = %d"%step)
	if step != None:
		plt.savefig('stiffness_Plot_%d.png'%step, transparent=True)
	else:
		plt.show()
	#plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return
########################################################################################
# Function to plot a target face from given quadedge cell
########################################################################################

def plotTargetFace(cell, numOfLayer,targetid, name = None, alpha = 0.8, Length=1.0,azim = None,elev= None):
	#calculating forces, stress-matrix and strain-matrix
	#cell.calculateVertexForce()
	#cell.calculateStressStrain()
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
	fig = plt.figure(frameon=False,figsize=(10,8))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	#fig.set_aspect(aspect='equal', adjustable='box')
	ax = Axes3D(fig)
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
		if face.getID() != targetid:
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
			ax.text(xCoord1,yCoord1,zCoord1, vertex.getID(),fontsize = 10)
			edge = edges.next()
		xlist.append(xlist[0])
		ylist.append(ylist[0])
		zlist.append(zlist[0])
		#ax.plot(xlist,ylist,zlist,'k',lw = 2)
		verts = [zip(np.array(xlist),np.array(ylist),np.array(zlist))]
		pc = Poly3DCollection(verts,alpha = 0.1,linewidths = 4,zorder = 1, facecolor = 'b')
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.plot(xlist,ylist,zlist,'k',lw = 2)
		ax.scatter(xlist,ylist,zlist,s = 8,c='k',marker = 'o',alpha = 1)
		ax.text(xmean,ymean,zmean, face.getID(),fontsize = 10)
		ax.scatter(xmean,ymean,zmean,s=8,c='k',marker = 'o',alpha = 1)
		#face = faces.next()
		break
		#if face.getID() == 1: break
	maxX = np.max(np.array(xlist))
	maxY = np.max(np.array(ylist))
	maxZ = np.max(np.array(zlist))
	minX = np.min(np.array(xlist))
	minY = np.min(np.array(ylist))
	minZ = np.min(np.array(zlist))
	###################################
	ax.set_xlim((minX,maxX))
	ax.set_ylim((minY,maxY))
	ax.set_zlim((minZ,maxZ))
	#plt.clf()
	if azim != None and elev != None:
		ax.view_init(azim = azim,elev=elev)
	if name == None:#plot the figure
		plt.show()
	else:
		plt.savefig(name, transparent = True)
	return
########################################################################################
# Function to plot Minimum gaussian curvature for each time step till given step
########################################################################################
def plotMinimumGaussianCurvature(step, numOfLayer=8.,save = False,resetids = True):
	###################
	fig = plt.figure(1)
	ax1 = fig.add_subplot(111)
	#ax2.axis('equal')
	###################
	minimumCurvature=[]
	for currentstep in range(step+1):
		cell = loadCellFromFile(currentstep, numOfLayer,resetids = resetids)
		curvatureArray = []
		faces = qd.CellFaceIterator(cell)
		face = faces.next()
		while (face != None):
			if face.getID() == 1:
				face  = faces.next()
				continue
			curvatureArray.append(face.getGaussianCurvature())
			face = faces.next()
		###################
		vertices = qd.CellVertexIterator(cell)
		vertex = vertices.next()
		while (vertex != None):
			#print "== Face ID : ", faceid, "=="
			if checkExternalVertex(vertex):
				vertex = vertices.next()
				continue
			curvatureArray.append(vertex.getGaussianCurvature())
			vertex = vertices.next()
		###################
		#print curvatureArray
		minimumCurvature.append(np.min(np.array(curvatureArray)))
		###################
	ax1.plot(minimumCurvature,'-x')
	ax1.set_ylabel("Minimum Gaussian Curvature",fontsize=15)
	ax1.set_xlabel("time Step",fontsize=15)
	if save:
		plt.savefig('minimum_gaussian_curvature.png', transparent = True)
	else:
		plt.show()
	return
########################################################################################
# Function to return face of given id
# if no face is found, None is returned
########################################################################################
def getFace(cell,targetid):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == targetid:
			break
		face = faces.next()
	return face
########################################################################################
# Function to return area of target cell and cells arround it
########################################################################################
def getAreaTargetRegion(cell, targetid):
	face = getFace(cell,targetid)
	edges = qd.FaceEdgeIterator(face)
	edge = edges.next()
	faceArea = face.getAreaOfFace()
	while edge != None:
		rightFace = edge.Right()
		faceArea += rightFace.getAreaOfFace()
		edge = edges.next()
	return faceArea
########################################################################################
# Function to return Surface area of total tissue
########################################################################################
def getSurfaceArea(cell):
	surfacearea= 0
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face =faces.next()
			continue
		surfacearea += face.getAreaOfFace()
		face =faces.next()
	return surfacearea
########################################################################################
# Function to plot Unit vectors
########################################################################################
def plotUnitVectors(cell, numOfLayer, step = None,save=False,directory=".", alpha = 0.8, Length=1.0):
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
	fig = plt.figure(frameon=False,figsize=(12,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	limfac = 0.7
	ax.set_xlim((-limfac*radius,limfac*radius))
	ax.set_ylim((-limfac*radius,limfac*radius))
	ax.set_zlim((-0.,2*limfac*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	########################################################################
	#                    Plotting the Stress vectors                       #
	########################################################################
	X = []
	Y = []
	Z = []
	U = []
	V = []
	W = []
	eigenvalue =[]
	eigenvalue1array = []
	eigenvalue2array = []
	eigenvalueratioarray = []
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getID() == 1:
			face  = faces.next()
			continue
		eigenvec1 = face.getUnitx()
		eigenvec2 = face.getUnity()
		eigenvec3 = face.getUnitz()
		#################################################
		eigenvalue1 = 1
		eigenvalue2 = 1
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		#########~~~ EIGEN VEC 1 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec1,0))
		V.append(qd.doublearray_getitem(eigenvec1,1))
		W.append(qd.doublearray_getitem(eigenvec1,2))
		#########~~~ EIGEN VEC 2 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec2,0))
		V.append(qd.doublearray_getitem(eigenvec2,1))
		W.append(qd.doublearray_getitem(eigenvec2,2))
		#########~~~ EIGEN VEC 3 ~~~#########
		#getting the centralised coordinate of centroid
		X.append(face.getXCentralised())
		Y.append(face.getYCentralised())
		Z.append(face.getZCentralised())
		#getting the vector headings
		U.append(qd.doublearray_getitem(eigenvec3,0))
		V.append(qd.doublearray_getitem(eigenvec3,1))
		W.append(qd.doublearray_getitem(eigenvec3,2))
		####################################################
		face = faces.next()
	from itertools import cycle
	color = cycle(['r','b','g'])
	for i in range(len(X)):
		ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=color.next(),length = 1,pivot='tail',zorder=4, linewidths = 2)
	########################################################
	#                 Plotting the Cell                    #
	########################################################
	###################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	xcenarray = []
	ycenarray = []
	zcenarray = []
	while (face != None):
		if face.getID() == 1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		xcenarray.append(face.getXCentralised())
		ycenarray.append(face.getYCentralised())
		zcenarray.append(face.getZCentralised())
		#print face.getZCentralised()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		#normalised : 
		color = 'w'
		#ax.add_collection3d(Poly3DCollection(verts,alpha = alpha,facecolors = color,linewidths=1,zorder=0))
		ax.plot(xlist,ylist, zlist, lw = 2, c = 'k')
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1], 'o',c='r')
		#ax.text(xcenarray[-1], ycenarray[-1],zcenarray[-1], face.getID())
		face = faces.next()
	if step != None:
		plt.title("x = r, y = b, z = g, Time = %d"%step)
	if save:
		plt.savefig(directory+r'/unitvec_plot_%03d.png'%step, transparent=True)
	else:
		plt.show()
	plt.close()
	return
########################################################################################
# Function to plot spontaneous mean curvature for a quadedge cell
########################################################################################
def plotSpontaneousMeanCurvatureTriangulation(cell, numOfLayer,save=True,threshold=-0.0, step = None, alpha = 0.8, Length=1.0,azim = -60.,elev=60):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(20,16))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	#fig.set_aspect(aspect='equal', adjustable='box')
	ax = Axes3D(fig)
	ax.set_xlim((-.7*radius,.7*radius))
	ax.set_ylim((-.7*radius,.7*radius))
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
	maxcurvature = cell.getMeanCurvatureWidth()
	mincurvature = -1.*cell.getMeanCurvatureWidth()
	###################################
	######### Color Map
	###################################
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=mincurvature, vmax=maxcurvature)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	####################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID() ==1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		#print "== Face ID : ", faceid, "=="
		xmean = face.getXCentralised()
		ymean = face.getYCentralised()
		zmean = face.getZCentralised()
		################################
		edges = qd.FaceEdgeIterator(face)
		ax.scatter(xmean,ymean,zmean,s = 50, c=scalarMap.to_rgba(face.getInitialMeanCurvature()))
		edge = edges.next()
		while edge != None:
			vertexOrg = edge.Org()
			vertexDest = edge.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			ax.plot([xmean, xCoord1,xCoord2],[ymean, yCoord1,yCoord2],[zmean, zCoord1,zCoord2], 'k', lw = 1)
			if checkExternalVertex(vertexOrg):
				ax.scatter(xCoord1,yCoord1,zCoord1,s = 50, c='k')
			else:
				ax.scatter(xCoord1,yCoord1,zCoord1,s = 50, c=scalarMap.to_rgba(vertexOrg.getInitialMeanCurvature()))
			edge = edges.next()
		################################
		face = faces.next()
	################################
	ax.view_init(azim = azim,elev=elev)
	###################################
	# color bar
	###################################
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.84, 0.2, 0.04, 0.55])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	clrbar.ax.tick_params(labelsize=25)
	clrbar.set_label("Spontaneous Curvature", fontsize = 30)
	if save:
		plt.savefig('spontaneous_curvature_layer=%d_elev=%d_azim=%d_step=%03d.png'%(numOfLayer,elev,azim,step), transparent = True)
	###################################
	plt.close()
	return
########################################################################################
# Function to calculate weighted gaussian curvature
########################################################################################
def getFaceWeightedGaussianCurvature(face):
	edges = qd.FaceEdgeIterator(face)
	edge = edges.next()
	curvature = [face.getGaussianCurvature()]
	areaMixed = [face.getAreaMixed()]
	while edge != None:
		vertexDest = edge.Dest()
		curvature.append(vertexDest.getGaussianCurvature())
		areaMixed.append(vertexDest.getAreaMixed())
		#####################
		edge = edges.next()
	return np.sum(np.multiply(curvature,areaMixed))/np.sum(areaMixed)#returning weighted curvature
########################################################################################
# Function to calculate weighted mean curvature
########################################################################################
def getFaceWeightedMeanCurvature(face):
	edges = qd.FaceEdgeIterator(face)
	edge = edges.next()
	curvature = [face.getMeanCurvature()]
	areaMixed = [face.getAreaMixed()]
	while edge != None:
		vertexDest = edge.Dest()
		curvature.append(vertexDest.getMeanCurvature())
		areaMixed.append(vertexDest.getAreaMixed())
		#####################
		edge = edges.next()
	return np.sum(np.multiply(curvature,areaMixed))/np.sum(areaMixed)#returning weighted curvature

########################################################################################
# Function to plot spontaneous mean curvature for a quadedge cell
########################################################################################
def plotMeanCurvatureTriangulation(cell, numOfLayer,save = True, threshold=-0.0, step = None, alpha = 0.8, Length=1.0,azim = -60.,elev=60):
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np
	import matplotlib.pyplot as plt
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(20,16))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	#fig.set_aspect(aspect='equal', adjustable='box')
	ax = Axes3D(fig)
	ax.set_xlim((-.7*radius,.7*radius))
	ax.set_ylim((-.7*radius,.7*radius))
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
	maxcurvature = cell.getMeanCurvatureWidth()
	mincurvature = -1.*cell.getMeanCurvatureWidth()
	###################################
	######### Color Map
	###################################
	jet = cm = plt.get_cmap('viridis') 
	cNorm  = colors.Normalize(vmin=mincurvature, vmax=maxcurvature)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	####################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID() ==1:
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		#print "== Face ID : ", faceid, "=="
		xmean = face.getXCentralised()
		ymean = face.getYCentralised()
		zmean = face.getZCentralised()
		################################
		edges = qd.FaceEdgeIterator(face)
		ax.scatter(xmean,ymean,zmean,s = 50, c=scalarMap.to_rgba(face.getMeanCurvature()))
		edge = edges.next()
		while edge != None:
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edge.Org()
			vertexDest = edge.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			ax.plot([xmean, xCoord1,xCoord2],[ymean, yCoord1,yCoord2],[zmean, zCoord1,zCoord2], 'k', lw = 1)
			if checkExternalVertex(vertexOrg):
				ax.scatter(xCoord1,yCoord1,zCoord1,s = 50, c='k')
			else:
				ax.scatter(xCoord1,yCoord1,zCoord1,s = 50, c=scalarMap.to_rgba(vertexOrg.getMeanCurvature()))
			edge = edges.next()
		################################
		face = faces.next()
	################################
	ax.view_init(azim = azim,elev=elev)
	###################################
	# color bar
	###################################
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.84, 0.2, 0.04, 0.55])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	clrbar.ax.tick_params(labelsize=25)
	clrbar.set_label("Mean Curvature", fontsize = 30)
	if save:
		plt.savefig('mean_curvature_layer=%d_elev=%d_azim=%d_step=%03d.png'%(numOfLayer,elev,azim,step), transparent = True)
	###################################
	plt.close()
	return
########################################################################################
# Function to plot gaussian curvature surface for a quadedge cell
########################################################################################
###################################################################
def plotGaussianCurvatureSurface(cell, numOfLayer,threshold=-10000.0,
				ids= False,step = None, save= False, name= None,fileformat = 'pdf',
				alpha = 0.8, Length=1.0,azim = -60.,elev=60,
				zaxisoffset=0.3,
				colormap = 'viridis'):
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
	#plotting part
	fig = plt.figure(frameon=False,figsize=(20,16))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	#fig.set_aspect(aspect='equal', adjustable='box')
	ax = Axes3D(fig)
	limfac =1.-zaxisoffset
	ax.set_xlim((-limfac*radius,limfac*radius))
	ax.set_ylim((-limfac*radius,limfac*radius))
	ax.set_zlim((-0.,0.7*2*limfac*radius))
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
	meanpointx = []
	meanpointy = []
	meanpointz = []
	#ax.scatter(np.mean(meanpointx),np.mean(meanpointy),np.mean(meanpointz),s = 60, c= 'm')
	#targetface = getFace(cell, 135)
	#targetx = targetface.getXCentralised()
	#targety = targetface.getYCentralised()
	#targetz = targetface.getZCentralised()
	#ax.scatter(targetx,targety,targetz,s = 60, c= 'g')
	#ax.plot([targetx,np.mean(meanpointx)],[targety,np.mean(meanpointy)],[targetz,np.mean(meanpointz)],lw = 4,c= 'r')
	
	###################
	curvaturedict = {}
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID() == 1 or checkExternalFace(face):
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		curvaturedict[faceid] = getFaceWeightedGaussianCurvature(face)#getFaceWeightedMeanCurvature(face)
		################################
		face = faces.next()
	###################
	curvatureArray = np.array(curvaturedict.values())
	maxcurvature = np.max(curvatureArray)
	mincurvature = -1*maxcurvature#np.min(curvatureArray)
	###################################
	######### Color Map
	###################################
	jet = cm = plt.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=mincurvature, vmax=maxcurvature)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	###################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		faceid = face.getID()#grabbing face id
		if face.getID() == 1:
			face = faces.next()
			continue
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		if checkExternalFace(face):
				#gausscurve = curvaturedict[faceid]
				facecolor = 'slategray'
		else:
				gausscurve = curvaturedict[faceid]
				if gausscurve<threshold:
					facecolor = 'r'
				else:
					facecolor =  scalarMap.to_rgba(gausscurve)
		#facecolor = 'b'
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = facecolor)
		ax.plot(xlist,ylist,zlist,'k', lw = 1)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		# adding ids for face
		if ids: 
			ax.text(xmean,ymean,zmean,faceid)
		face = faces.next()
	ax.view_init(azim = azim,elev=elev)
	###################################
	# color bar
	###################################
	import matplotlib.ticker
	print mincurvature, maxcurvature
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.55])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax,
						  ticks = np.linspace(mincurvature,maxcurvature,3),
						  format=OOMFormatter(-2, mathText=True))#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label(r"Gaussian curvature, $K$")
	#plt.show()
	#plt.suptitle("Step =%03d"%step,fontsize = 30)
	if save:
		if name:
			plt.savefig(name+'.'+fileformat,format = fileformat,bbox_inches='tight', transparent=True)
		else:
			plt.savefig('gaussianCurvatureSurface-weighted'+'.'+fileformat,format = fileformat,bbox_inches='tight', transparent=True)
	#print table.draw()
	return
###################################################################
# to plot surface with mean curvature
###################################################################
def plotMeanCurvatureSurface(cell, numOfLayer,threshold=-1000.0,
				ids= False,step = None, save= False, name=None,fileformat='pdf',
				alpha = 0.8, Length=1.0,azim = -60.,elev=60,
				zaxisoffset=0.3,
				colormap = 'viridis'):
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
	#plotting part
	fig = plt.figure(frameon=False,figsize=(20,16))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	#fig.set_aspect(aspect='equal', adjustable='box')
	ax = Axes3D(fig)
	limfac =1.-zaxisoffset
	ax.set_xlim((-limfac*radius,limfac*radius))
	ax.set_ylim((-limfac*radius,limfac*radius))
	ax.set_zlim((-0.,0.7*2*limfac*radius))
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
	meanpointx = []
	meanpointy = []
	meanpointz = []
	unneccesary = False
	if unneccesary:
		print "NOT"
		edge = getEdge(cell,211,210)
		edgenext = edge
		####grabbing the origin of edge####
		#centralised coordiante
		vertexOrg = edge.Org()
		vertexDest = edge.Dest()
		#print vertex.getID()
		xCoord1 = vertexOrg.getXcoordinate()
		yCoord1 = vertexOrg.getYcoordinate()
		zCoord1 = vertexOrg.getZcoordinate()
		xCoord2 = vertexDest.getXcoordinate()
		yCoord2 = vertexDest.getYcoordinate()
		zCoord2 = vertexDest.getZcoordinate()
		meanpointx.append(xCoord2)
		meanpointy.append(yCoord2)
		meanpointz.append(zCoord2)
		ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'m', lw = 4)
		while True:
			edgenext = edgenext.Rprev()
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edgenext.Org()
			vertexDest = edgenext.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			meanpointx.append(xCoord2)
			meanpointy.append(yCoord2)
			meanpointz.append(zCoord2)
			ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'r', lw = 4)
			########################################################
			edgenext = edgenext.Lnext()
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edgenext.Org()
			vertexDest = edgenext.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			meanpointx.append(xCoord2)
			meanpointy.append(yCoord2)
			meanpointz.append(zCoord2)
			ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'r', lw = 4)
			########################################################
			edgenext = edgenext.Lnext()
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edgenext.Org()
			vertexDest = edgenext.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			meanpointx.append(xCoord2)
			meanpointy.append(yCoord2)
			meanpointz.append(zCoord2)
			ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'r', lw = 4)
			###############################################################
			edgenext = edgenext.Rprev()
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edgenext.Org()
			vertexDest = edgenext.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			meanpointx.append(xCoord2)
			meanpointy.append(yCoord2)
			meanpointz.append(zCoord2)
			ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'r', lw = 4)
			########################################################
			edgenext = edgenext.Lnext()
			####grabbing the origin of edge####
			#centralised coordiante
			vertexOrg = edgenext.Org()
			vertexDest = edgenext.Dest()
			#print vertex.getID()
			xCoord1 = vertexOrg.getXcoordinate()
			yCoord1 = vertexOrg.getYcoordinate()
			zCoord1 = vertexOrg.getZcoordinate()
			xCoord2 = vertexDest.getXcoordinate()
			yCoord2 = vertexDest.getYcoordinate()
			zCoord2 = vertexDest.getZcoordinate()
			meanpointx.append(xCoord2)
			meanpointy.append(yCoord2)
			meanpointz.append(zCoord2)
			ax.plot([xCoord1,xCoord2],[yCoord1,yCoord2], [zCoord1,zCoord2],'r', lw = 4)
			#############################################################
			if edgenext.Org().getID() == edge.Org().getID():
				break
	#ax.scatter(np.mean(meanpointx),np.mean(meanpointy),np.mean(meanpointz),s = 60, c= 'm')
	#targetface = getFace(cell, 135)
	#targetx = targetface.getXCentralised()
	#targety = targetface.getYCentralised()
	#targetz = targetface.getZCentralised()
	#ax.scatter(targetx,targety,targetz,s = 60, c= 'g')
	#ax.plot([targetx,np.mean(meanpointx)],[targety,np.mean(meanpointy)],[targetz,np.mean(meanpointz)],lw = 4,c= 'r')
	
	###################
	curvaturedict = {}
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		if face.getID() == 1 or checkExternalFace(face):
			face  = faces.next()
			continue
		faceid = face.getID()#grabbing face id
		curvaturedict[faceid] = getFaceWeightedMeanCurvature(face)#getFaceWeightedMeanCurvature(face)
		################################
		face = faces.next()
	###################
	curvatureArray = np.array(curvaturedict.values())
	maxcurvature = np.max(curvatureArray)
	mincurvature = np.min(curvatureArray)
	###################################
	######### Color Map
	###################################
	jet = cm = plt.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=mincurvature, vmax=maxcurvature)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	###################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		faceid = face.getID()#grabbing face id
		if face.getID() == 1:
			face = faces.next()
			continue
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		if checkExternalFace(face):
				#gausscurve = curvaturedict[faceid]
				facecolor = 'slategray'
		else:
				gausscurve = curvaturedict[faceid]
				if gausscurve<threshold:
					facecolor = 'r'
				else:
					facecolor =  scalarMap.to_rgba(gausscurve)
		#facecolor = 'b'
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = facecolor)
		ax.plot(xlist,ylist,zlist,'k', lw = 1)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		# adding ids for face
		if ids: 
			ax.text(xmean,ymean,zmean,faceid)
		face = faces.next()
	ax.view_init(azim = azim,elev=elev)
	###################################
	# color bar
	###################################
	import matplotlib.ticker
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.55])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax,
						  ticks = np.linspace(mincurvature,maxcurvature,3),
						  format=OOMFormatter(-2, mathText=True))#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label(r"Mean Curvature, $H$")
	print mincurvature, maxcurvature
	#clrbar.ax.tick_params(labelsize=30)
	#clrbar.ax.set_yticklabels(['{:.3f}'.format(x) for x in np.linspace(mincurvature,maxcurvature,4)], fontsize=30)
	#formatter = matplotlib.ticker.ScalarFormatter(useMathText=False)
	#clrbar.ax.yaxis.set_major_formatter(formatter)
	#clrbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


	#plt.show()
	#plt.suptitle("Step =%03d"%step,fontsize = 30)
	if save:
		if name:
			plt.savefig(name+'.'+fileformat,format =fileformat,bbox_inches='tight', transparent=True)
		else:
			plt.savefig('meanCurvatureSurface-weighted'+'.'+fileformat,format = fileformat,bbox_inches='tight', transparent=True)
	#print table.draw()
	return
####################################################################################################
# Get array of all vertices on boundary of primordia
####################################################################################################
def getPrimordiaBoundaryVertexList(cell, targetid, large=False):
	########################################################
	"""
	Get a list of all vertices arround a primordia
	targetid : center of primordia
	large : if true, large primordia is calculated for (2 layers)
			if false, small primordia (1 layer)
	"""
	face = getFace(cell, targetid)
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
# Give list of all edge arround primiordia
###################################################################
def getPrimordiaBoundaryEdgeList(cell, targetid, large=False):
	########################################################
	"""
	Get a list of all edges arround a primordia
	targetid : center of primordia
	large : if true, large primordia is calculated for (2 layers)
			if false, small primordia (1 layer)
	"""
	face = getFace(cell, targetid)
	edge = face.getEdge()
	edgeList = []
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
			edgeList.append(edge)
			#####################
			for _ in range(2):
				edge = edge.Lnext()
				edgeList.append(edge)
			#####################
			edge = edge.Rprev()
			edgeList.append(edge)
			#####################
			edge = edge.Lnext()
			edgeList.append(edge)
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
			edgeList.append(edge)
		#####################
		for _ in range(5):
			edge = edge.Rprev()
			edgeList.append(edge)
			####################
			edge = edge.Lnext()
			edgeList.append(edge)
			####################
			edge = edge.Lnext()
			edgeList.append(edge)
			####################
	return edgeList
##################################################################
# Give list of all face arround primiordia
###################################################################
def getPrimordiaBoundaryFaceList(cell, targetid, large=False):
	def addFaceList(faceList,faceidlist, face):
		if not face.getID() in faceidlist:
			faceList.append(face)
			faceidlist.append(face.getID())
		return faceList,faceidlist
	########################################################
	"""
	Get a list of all Faces arround a primordia
	targetid : center of primordia
	large : if true, large primordia is calculated for (2 layers)
			if false, small primordia (1 layer)
	"""
	face = getFace(cell, targetid)
	edge = face.getEdge()
	faceidlist = []
	faceList = []
	########################################################
	# tracing the primordial boundary
	########################################################
	if not large:
		#########################################################
		# Smaller Promordia
		#########################################################
		###################
		# First Layer
		###################
		for _ in range(1):
			edge = edge.Rprev()
		###############################################################
		# listing Primordia vertex starts here
		###############################################################
		for _ in range(3):
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
		#####################
		for _ in range(5):
			edge = edge.Rprev()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
	else:
		#########################################################
		# Smaller Promordia
		#########################################################
		###################
		# First Layer
		###################
		for _ in range(1):
			edge = edge.Rprev()
		###############################################################
		# listing Primordia vertex starts here
		###############################################################
		for _ in range(3):
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
		#####################
		for _ in range(5):
			edge = edge.Rprev()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
		###############################################################
		###################
		# Second Layer
		###################
		face = getFace(cell, targetid)
		edge = face.getEdge()
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
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			#####################
			for _ in range(2):
				edge = edge.Lnext()
				faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			#####################
			edge = edge.Rprev()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			#####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
	return faceList
########################################################################################
# Getting list of all cells
########################################################################################
def getFaceList(cell):
	facelist = []
	faces1 = qd.CellFaceIterator(cell)
	face1 = faces1.next()
	while face1 != None:
		if face1.getID() == 1:
			face1 = faces1.next()
			continue
		facelist.append(face1)
		face1 = faces1.next()
	return facelist
########################################################################################
# Getting list of all primordia cells
########################################################################################
def getPrimordiaFaces(cell,targetid, large = False):
	########################################################
	"""
	Get a list of all faces in primordia
	targetid : center of primordia
	large : if true, large primordia is calculated for (2 layers)
			if false, small primordia (1 layer)
	"""
	face = getFace(cell, targetid)
	edge = face.getEdge()
	faceList = [face]
	########################################################
	# tracing the primordial boundary
	########################################################
	# Smaller Primordia
	#########################################################
	for _ in range(1):
		edge = edge.Rprev()
	###############################################################
	# listing Primordia vertex starts here
	###############################################################
	for _ in range(3):
		edge = edge.Lnext()
	faceList.append(edge.Left())
	#####################
	for _ in range(5):
		edge = edge.Rprev()
		faceList.append(edge.Left())
		####################
		edge = edge.Lnext()
		####################
		edge = edge.Lnext()
		####################
	if large:
		edge = face.getEdge()
		#########################################################
		# Larger Proimordia
		#########################################################
		for _ in range(2):
			edge = edge.Rprev()
		for _ in range(3):
			edge = edge.Lnext()
		####################################################
		# listing Primordia vertex starts here
		####################################################
		for _ in range(6):
			edge = edge.Rprev()
			#####################
			faceList.append(edge.Left())
			#####################
			for _ in range(2):
				edge = edge.Lnext()
			#####################
			edge = edge.Rprev()
			#####################
			faceList.append(edge.Left())
			#####################
			edge = edge.Lnext()
	return faceList
########################################################################################
#   Calculating Growth in the cells
# cell1 = cell_i
# cell2 = cell_{i+1}
########################################################################################
def getGrowthOfCells(cell1, cell2):
	faces1 = qd.CellFaceIterator(cell1)
	face1 = faces1.next()
	growthdict = {}
	##########################
	"""
	while face1 != None:
		if face1.getID() == 1:
			face1 = faces1.next()
			continue
		##########################
		faceid1 = face1.getID()
		face2 = getFace(cell2, faceid1)
		##########################
		cfm1 = [[qd.getCurrentFormMatrix(face1, 0, 0 ),qd.getCurrentFormMatrix(face1, 0, 1 )],
				[qd.getCurrentFormMatrix(face1, 1, 0 ),qd.getCurrentFormMatrix(face1, 1, 1 )]]
		##########################
		cfm2 = [[qd.getCurrentFormMatrix(face2, 0, 0 ),qd.getCurrentFormMatrix(face2, 0, 1 )],
				[qd.getCurrentFormMatrix(face2, 1, 0 ),qd.getCurrentFormMatrix(face2, 1, 1 )]]
		##########################
		#incomplete !!! 
		##########################
	"""
	return

########################################################################################
# Getting separate lists of faces in primorida
########################################################################################
def getSeparatePrimordiaBoundaryFaceList(cell, targetid, large=False):
	def addFaceList(faceList,faceidlist, face):
		if not face.getID() in faceidlist:
			faceList.append(face)
			faceidlist.append(face.getID())
		return faceList,faceidlist
	########################################################
	"""
	Get a list of all Faces arround a primordia
	targetid : center of primordia
	large : if true, large primordia is calculated for (2 layers)
			if false, small primordia (1 layer)
	"""
	face = getFace(cell, targetid)
	edge = face.getEdge()
	faceidlist = []
	faceList = []
	faceList2 = []
	########################################################
	# tracing the primordial boundary
	########################################################
	if not large:
		#########################################################
		# Smaller Promordia
		#########################################################
		###################
		# First Layer
		###################
		for _ in range(1):
			edge = edge.Rprev()
		###############################################################
		# listing Primordia vertex starts here
		###############################################################
		for _ in range(3):
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
		#####################
		for _ in range(5):
			edge = edge.Rprev()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
	else:
		#########################################################
		#Larger Promordia
		#########################################################
		###################
		# First Layer
		###################
		for _ in range(1):
			edge = edge.Rprev()
		###############################################################
		# listing Primordia vertex starts here
		###############################################################
		for _ in range(3):
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
		#####################
		for _ in range(5):
			edge = edge.Rprev()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
			edge = edge.Lnext()
			faceList,faceidlist = addFaceList(faceList,faceidlist, edge.Right())
			####################
		###############################################################
		###################
		# Second Layer
		###################
		face = getFace(cell, targetid)
		edge = face.getEdge()
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
			faceList2,faceidlist = addFaceList(faceList2,faceidlist, edge.Right())
			#####################
			for _ in range(2):
				edge = edge.Lnext()
				faceList2,faceidlist = addFaceList(faceList2,faceidlist, edge.Right())
			#####################
			edge = edge.Rprev()
			faceList2,faceidlist = addFaceList(faceList2,faceidlist, edge.Right())
			#####################
			edge = edge.Lnext()
			faceList2,faceidlist = addFaceList(faceList2,faceidlist, edge.Right())
	return [faceList,faceList2]
##############################################################################
# Calculate dilation in cells
##############################################################################
def calculateDilation(cell1,cell2):
	faces = qd.CellFaceIterator(cell1)
	face1 = faces.next()
	#######################################
	while face1 != None:
		if face1.getID() == 1:
			face1 = faces.next()
			continue
		#######################################
		face2 = getFace(cell2,face1.getID())
		#######################################
		#calculating the dilation in face
		face1.getRotatedGrowthMatrix(face2)
		#######################################
		face1 = faces.next()
	#######################################
	return
####################################################################################################################
# Calculating the max time step for target surface area
####################################################################################################################
def getTimeStep(targetArea, endStep, startStep=1, stepsize = 10,resetids = True):
	####################################################
	import gc
	for step in range(startStep, endStep+1,stepsize):
		if not os.path.isfile("qdObject_step=%03d.obj"%step):
			return step-stepsize, tissueSurfaceArea
		################################################
		cell = loadCellFromFile(step,resetids = resetids)
		################################################
		tissueSurfaceArea = getSurfaceArea(cell)
		if (tissueSurfaceArea > targetArea):
			gc.collect()
			for calstep in range(step-1,step-stepsize-1,-1):
					cell = loadCellFromFile(calstep,resetids = resetids)
					tissueSurfaceArea = getSurfaceArea(cell)
					if (tissueSurfaceArea <= targetArea):
						gc.collect()
						cell = loadCellFromFile(calstep+1,resetids = resetids)
						tissueSurfaceArea = getSurfaceArea(cell)
						return calstep+1,tissueSurfaceArea
		################################################
		gc.collect()
	return endStep,tissueSurfaceArea
####################################################################################################################
# Calculating the surface area for a given timestep
####################################################################################################################
def getSurfaceAreaTimeStep(step,resetids = True):
	####################################################
	if not os.path.isfile("qdObject_step=%03d.obj"%step):
		return 0.
	################################################
	cell = loadCellFromFile(step,resetids = resetids)
	################################################
	tissueSurfaceArea = getSurfaceArea(cell)
	################################################
	return tissueSurfaceArea
####################################################################################################################
# get Zone Number of a face, by default 3 zones
####################################################################################################################
def getZoneNum(cell, face,numOfLayer,zoneOffSet = 8,zones = 3.):
	#the radius of circle to be projected on
	Length=1.
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length
	# total cells
	totalCellNum = float(cell.countFaces())
	zoneCellNum = totalCellNum/zones
	zoneNum = int((face.getID()-zoneOffSet)/zoneCellNum)+1
	if zoneNum == 0:
		zoneNum = 1
	elif zoneNum>zones:
		zoneNum = zones
	return zoneNum
####################################################################################################################
# set Zone Number of the faces
####################################################################################################################
def setZoneNum(cell, numOfLayer, zoneOffSet  = 6,zones = 3.):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		face.setZone(getZoneNum(cell, face, numOfLayer,zoneOffSet=zoneOffSet, zones=zones))
		face = faces.next()
	return
####################################################################################################################
# set growth rates for one zone
####################################################################################################################
def setZonalKappa(cell,zone, kappa):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getZone() == zone:
			face.setKappa(kappa)
		face = faces.next()
	return
####################################################################################################################
# set division threshold for one zone
####################################################################################################################
def setZonalDivisionThreshold(cell,zone, divisionthreshold):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		if face.getZone() == zone:
			face.setDivisionThreshold(divisionthreshold)
		face = faces.next()
	return
##########################################################################################
#       Function to Plot division Threshold of the cells
##########################################################################################
def plotDivisionThresholdSurface(cell, numOfLayer, name=None, alpha = 0.5, Length=1.0,vmax=0,azim = 0, elev = 90):
	#import the libraries
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib as mpl
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection
	from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.colors as colors
	import matplotlib.cm as cmx
	#limits of the plot
	radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
	#plotting part
	fig = plt.figure(frameon=False,figsize=(10,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-.7*radius,.7*radius))
	ax.set_ylim((-.7*radius,.7*radius))
	ax.set_zlim((0*radius,1.4*radius))
	ax.axis('off')
	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	##########################################################
	#Get the Max/Min Of Division Threshold
	##########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	divisionthreshold = []
	while (face != None):
		faceid = face.getID()#grabbing face id
		if faceid ==1 :
			face = faces.next()
			continue
		divisionthreshold.append(face.getDivisionThreshold())
		face = faces.next()
	##########################################################
	vmax = max(divisionthreshold)
	vmin = min(divisionthreshold)
	if vmax == vmin:
		vmin = 0.
	##########################################################
	#### Making the COLOR BAR #########################
	##########################################################
	jet = cm = plt.get_cmap('spring') 
	cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	##########################################################
	##########################################################
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while (face != None):
		faceid = face.getID()#grabbing face id
		if faceid == 1:
			face  = faces.next()
			continue
		xlist = []
		ylist = []
		zlist = []
		xproj = []
		yproj = []
		zproj = []
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
		verts = [zip(xlist, ylist,zlist)]
		#adding to 3d plot
		dthreshold = face.getDivisionThreshold()
		color = scalarMap.to_rgba(dthreshold)
		pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		face = faces.next()
				#if face.getID() == 1: break
	#ax.axis("off")
	ax.view_init(elev=elev, azim=azim)
	########################################################################
	scalarMap._A = []
	cbar_ax = fig.add_axes([0.9, 0.25, 0.03, 0.5])
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax,
						  ticks = np.linspace(vmin,vmax,3),format = '%.1f')#,orientation='horizontal',cax = cbar_ax)
	########################################################################
	clrbar.set_label(r"Division Threshold")
	if name == None:#plot the figure
		plt.show()
	else:
		plt.savefig(name, transparent = True)
	plt.close()
	return
####################################################################################################################
# calculating growth ratio
####################################################################################################################
def getGrowthRatio(numOfLayer, targetid ,endStep ,startStep = 1,stepsize = 5, resetids = True, growthrates= False):
	if not os.path.isfile("qdObject_step=001.obj"):
			return [0.,0.,0.,0.,0.,0.,0.,0.,0.]
	cell = loadCellFromFile(1,resetids = resetids)
	####################################################################################################################
	# fit function
	####################################################################################################################
	def fitLinFunc(t,m,c):
		return m*t + c
	import scipy.optimize as sop
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
		cell = loadCellFromFile(step,resetids=resetids)
		################################################
		primordiafacelist  = getPrimordiaFaces(cell, targetid, large = False)
		primordiaarea = 0.
		for face in primordiafacelist:
			primordiaarea += face.getAreaOfFace()
		################################################
		tissueSurfaceArea = getSurfaceArea(cell)
		################################################
		primordialface = getFace(cell, targetid)
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
	if grrowthrates:
		return  [fastareafit[0]/slowareafit[0], fastareafit[0], slowareafit[0]]
	else:
		return fastareafit[0]/slowareafit[0]
####################################################################################################################
# aspectratio
####################################################################################################################
def getAspectRatio(face):
	targetformmatrix = np.array([[qd.getTargetFormMatrix(face, 0,0),qd.getTargetFormMatrix(face, 0,1)],
								 [qd.getTargetFormMatrix(face, 1,0),qd.getTargetFormMatrix(face, 1,1)]
								 ])
	eig_vec,eig_val,u = np.linalg.svd(targetformmatrix)
	return eig_val[0]/eig_val[1]
###################################################################
# To plot the ellipse surface of the cell
###################################################################
def plotMatrixSurface(cell, numOfLayer, step = None, 
    alpha = 0.5, Length=1.0,azim = -60.,elev=60,
    format='eps',name = None, ids = False):
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
    fig = plt.figure(frameon=False,figsize=(12,12))
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
            ax.text(xmean,ymean,zmean, face.getID(),fontsize = 20)
        face = faces.next()
        #if face.getID() == 1: break
    #plt.clf()
    ##PLOTTING ELLIPSE
    numofface = cell.countFaces()
    for i in range(2,numofface+1):
        # Target Form
        TFellipsepoints = ep.getTargetFormEllipsePoints(cell,i)
        verts = [zip(np.array(TFellipsepoints[0])[0],
                     np.array(TFellipsepoints[1])[0], 
                     np.array(TFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = 'r',zorder = 100,label="TFM" if i == 0 else "")
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # Current Form
        CFellipsepoints = ep.getCurrentFormEllipsePoints(cell,i)
        verts = [zip(np.array(CFellipsepoints[0])[0],np.array(CFellipsepoints[1])[0], np.array(CFellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = alpha,linewidths=1,zorder = 1, facecolor = 'b')
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # Strain 
        strainellipsepoints = ep.getStrainEllipsePoints(cell,i)
        verts = [zip(np.array(strainellipsepoints[0])[0],np.array(strainellipsepoints[1])[0], np.array(strainellipsepoints[2])[0])]
        pc = Poly3DCollection(verts,alpha = alpha,linewidths=1,zorder = 120, facecolor = 'k')
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
    ax.view_init(azim = azim,elev=elev)
    scatter1_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='r', marker = 's',ms=20)
    scatter2_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='b', marker = 's',ms=20)
    scatter3_proxy = mpl.lines.Line2D([0],[0], linestyle="none", c='k', marker = 's',ms=20)
    ax.legend([scatter1_proxy, scatter2_proxy,scatter3_proxy], ['Target Shape', 'Current Shape',"Strain"], numpoints = 1,loc = 0)
    if name == None:#plot the figure
        plt.show()
    else:
        if format == None:
            plt.savefig(name, transparent = True)
        else:
            plt.savefig(name+"."+format, transparent = True, format=format,dpi=500)
    #plt.suptitle("Step =%03d"%step,fontsize = 30)
    #plt.savefig('initial_TFM_layer8.png', transparent=True)
    return