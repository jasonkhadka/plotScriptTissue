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
import string
#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['axes.titlesize'] = 22


####################################################################################################################
# Add subplot annotation
####################################################################################################################
def addAnnotation(subplot,n = 1):
	subplot.text(-0.1,1.05,string.ascii_lowercase[n],transform = subplot.transAxes,size = 26, weight = 'bold')
	return
####################################################################################################################
# get mean radius of the vertices
####################################################################################################################
def getMeanRadius(cell):
	vertices = qd.CellVertexIterator(cell)
	vertex = vertices.next()
	countVertex = cell.countVertices()
	radiusArray = np.zeros(countVertex)
	counter = 0
	while vertex != None:
		radiusArray[counter] = np.sqrt(vertex.getXcoordinate()**2+vertex.getYcoordinate()**2)
		counter += 1
		vertex = vertices.next()
	return radiusArray
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
			return step-stepsize, tissueSurfaceArea
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
###############################################################################################################
# Calculate tissue height
###############################################################################################################
def getTissueHeight(cell):
    vertices = qd.CellVertexIterator(cell)
    vertex = vertices.next()
    vertnum = cell.countVertices()
    zArray = np.zeros(vertnum)
    counter = 0
    while vertex!= None:
        zArray[counter] = vertex.getZcoordinate()
        counter += 1
        vertex = vertices.next()
    ######################################
    return np.max(zArray)
###############################################################################################################
# Calculate distance to the tip of the tissue
###############################################################################################################
def getFaceDistance(face1, face2):
	x1 = face1.getXCentralised()
	y1 = face1.getYCentralised()
	z1 = face1.getZCentralised()

	x2 = face2.getXCentralised()
	y2 = face2.getYCentralised()
	z2 = face2.getZCentralised()

	return np.sqrt((x2-x1)**2+
					(y2-y1)**2+
					(z2-z1)**2
				   )
###############################################################################################################
# top initialise the dict
###############################################################################################################
def initialiseDict(mydict, maxNum):
    mydict = {x:[] for x in range(maxNum+1)}
    return mydict
###############################################################################################################
def addToDict(face,dicts, values):
	dict1 = dicts[0]
	if dict1.has_key(face.getID()):#face exists already
		for d,v in zip(dicts,values):
			d[face.getID()].append(v)
	else:
		for d,v in zip(dicts,values):
			d[face.getID()]=[v]
	return
####################################################################################################################
# plot the surface with growth rate and stress
####################################################################################################################
def plotSurfaceStressGrowthRate(plotData,numOfLayer, eta = etacurrent,name=None, alpha = 0.5, Length=1.0,
    ids=False, azim = 0, elev = 7,color = 'c',format = 'eps',
    zaxisoffset=0.2):
    #import the libraries
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib as mpl
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np
    import matplotlib.pyplot as plt
    #limits of the plot
    radius = (numOfLayer>1)*(np.sqrt(3.)*(numOfLayer-1)-Length)+Length#the radius of circle to be projected on
    ###############################################################################
    # getting data information 
    # format of plotData : [tissueSurfaceAreaArray,tissueSurfaceAreaArray2,
    #						 heightArray, timeArray,timeArray2, dictArray]
    # #  dictArray = [absStressDict, stressRadialDict,
    #				 stressOrthoradialDict,growthRateDict, tipDistanceDict]
    ###############################################################################
    laststep = plotData[4][-1]
    ############################################################
    #getting the dicts
    dictArray = plotData[5]
    absstressDict = dictArray[0]# stresses are computed as total abs stress on a cell
    growthrateDict = dictArray[3]
    ############################################################
    maxstress = max(absStressDict[max(absStressDict.keys(), key = lambda k : max(absStressDict[k]))])
    minstress = min(absStressDict[min(absStressDict.keys(), key = lambda k : min(absStressDict[k]))])
    ############################################################
    maxgrowthrate = max(growthrateDict[max(growthrateDict.keys(), key = lambda k : max(growthrateDict[k]))])
    mingrowthrate = min(growthrateDict[min(growthrateDict.keys(), key = lambda k : min(growthrateDict[k]))])
    ###############################################################################
    # colorbar
    ###############################################################################
    # for stress
    jet1 = cm = plt.get_cmap('inferno') 
    cNorm1  = colors.Normalize(vmin=minstress, vmax=maxstress)
    scalarMapStress = cmx.ScalarMappable(norm=cNorm1, cmap=jet1)
    # for growth
    jet2 = cm = plt.get_cmap('cool') 
    cNorm2  = colors.Normalize(vmin=mingrowthrate, vmax=maxgrowthrate)
    scalarMapGrowth = cmx.ScalarMappable(norm=cNorm2, cmap=jet2)
    ###############################################################################
    #plotting part
    ###############################################################################
    fig = plt.figure(1,frameon=False,figsize=(12,(1.-zaxisoffset)*10))
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


    fig2 = plt.figure(2, frameon=False,figsize=(12,(1.-zaxisoffset)*10))
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
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ##iterating through the cell##
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
        pc = Poly3DCollection(verts,alpha = alpha,linewidths=1, facecolor = color)
        pc.set_edgecolor('k')
        ax.add_collection3d(pc)
        # adding ids for face
        if ids: 
            ax.text(xmean,ymean,zmean,faceid)
        face = faces.next()
                #if face.getID() == 1: break
    #ax.axis("off")
    ax.view_init(azim = azim,elev=elev)
    if name == None:#plot the figure
        plt.show()
    else:
        if format == None:
            plt.savefig(name, transparent = True)
        else:
            plt.savefig(name+"."+format, transparent = True, format=format)
    #plt.clf()
    return
####################################################################################################################
# Calculating and plotting mean stress and growth
####################################################################################################################
def getGrowthRateStress(numOfLayer, endStep,eta,startStep=0,stepsize= 1,maxarea = None, areastep = 20,
	startarea = None, endarea = 850,resetids = False):
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
	topFaceID = 3*numOfLayer*(numOfLayer-1)+2#the id of cell at the top of the dome 
	#######################################################################
	# Starting the Calculation
	######################################################################
	#######################################################################
	laststep = 1
	plotargs = {"markersize": 10, "capsize": 10,"elinewidth":3,"markeredgewidth":2}
	#######################################################################
	heightArray = []
	tissueSurfaceAreaArray = []
	tissueSurfaceAreaArray2 = []
	timeArray = []
	timeArray2 = []
	dhdAArray = []
	volumeArray = []
	radiusMeanArray = []
	radiusVarArray = []
	topdistance = []
	################################################################
	absStressDict = {}
	stressRadialDict = {}
	stressOrthoradialDict = {}
	growthRateDict = {}
	tipDistanceDict = {}
	dictArray = [absStressDict, stressRadialDict, stressOrthoradialDict,growthRateDict, tipDistanceDict]
	###################################################
	if not startarea:#no startarea given
		startarea = int(initialTissueSurfaceArea)
	###################################################
	listsurfacearea = np.linspace(startarea,endarea,15)
	for steparea in listsurfacearea:
		step,tissueSurfaceArea = getTimeStep(steparea, endStep, laststep, stepsize = 10)
		step2,tissueSurfaceArea2 = getTimeStep(steparea+areastep, endStep, step, stepsize = 10)
		if step == step2:break
		########################################################################
		if not (os.path.isfile("qdObject_step=%03d.obj"%step) or os.path.isfile("qdObject_step=%03d.obj"%step2) ):#check if file exists
			break
		cell = sf.loadCellFromFile(step,resetids = resetids)
		cell2 = sf.loadCellFromFile(step2,resetids = resetids)
		################################################
		topFace = sf.getFace(cell2, topFaceID)
		################################################
		cell2.calculateStressStrain()
		cell2.setRadialOrthoradialVector(topFace)
		cell2.setRadialOrthoradialStress()
		################################################
		dTissueSurfaceArea = tissueSurfaceArea2-tissueSurfaceArea
		################################################
		#computing the cell growth rate
		################################################
		faces = qd.CellFaceIterator(cell2)
		face2 = faces.next()
		while face2 != None:
			faceid = face2.getID()
			face1 = sf.getFace(cell, faceid)
			############################################
			if face1 == None:
				absstress = np.abs(face2.getStressEigenValue1())+np.abs(face2.getStressEigenValue2())
				stressRadial = face2.getRadialStress()
				stressOrthoradial = face2.getOrthoradialStress()
				tipdistance = getFaceDistance(face2, topFace)
				growthrate  = np.nan
				############################################
				addToDict(face2, dictArray, [absstress, stressRadial, stressOrthoradial, growthrate, tipdistance])
				############################################
				face2 = faces.next()
				continue
			############################################
			facearea1 = face1.getAreaOfFace()
			facearea2 = face2.getAreaOfFace()
			############################################
			# checking if face2 has recently divided
			# if face divided recently, facearea is 
			# divided in near half
			############################################
			if (facearea2-facearea1) < -0.3*facearea1:
				# if this is true, then face1 has
				# seen drastic area decrease 
				# which would mean cell division
				growthrate = np.nan
			else:
				############################################
				#if no recent division
				############################################
				dfacearea = facearea2-facearea1
				growthrate = dfacearea/(facearea1*dTissueSurfaceArea)
			############################################
			absstress = np.abs(face2.getStressEigenValue1())+np.abs(face2.getStressEigenValue2())
			stressRadial = face2.getRadialStress()
			stressOrthoradial = face2.getOrthoradialStress()
			tipdistance = getFaceDistance(face2, topFace)
			############################################
			addToDict(face2, dictArray, [absstress, stressRadial, stressOrthoradial, growthrate, tipdistance])
			############################################
			face2 = faces.next()
			############################################
		################################################
		height = getTissueHeight(cell)
		heightArray.append(height)
		timeArray.append(step)
		timeArray2.append(step2)
		tissueSurfaceAreaArray.append(tissueSurfaceArea)
		tissueSurfaceAreaArray2.append(tissueSurfaceArea2)
		########################################################################
		print step2, tissueSurfaceArea, height
	########################################################################
	return [tissueSurfaceAreaArray,tissueSurfaceAreaArray2, heightArray, timeArray,timeArray2, dictArray]
####################################################################################################################################################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--startarea", help="Start of area",default =None, type = int)
parser.add_argument('-e',"--endarea", help="area end",default = 1000, type = int)
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
parser.add_argument('-d',"--areastep", help="area step for calculating the growth in cell area", type = int,default = 10)
parser.add_argument('-j',"--jobid", help="jobid", type = int,default = None)
parser.add_argument("-r","--resetids", help = "if option is used, the figures are not normalised", action= "store_true")

## Getting the arguments 
args = parser.parse_args()
#location = args.location
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
startarea = args.startarea
endarea =args.endarea
jobid = args.jobid
resetids = args.resetids

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
################################################################################
cwd = os.getcwd()#getting the current working directory
zone = cwd.split('/')[-2]
################################################################################
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
areaHeight = fig.add_subplot(221)
absstressplot = fig.add_subplot(222)
stressRadOrthplot = fig.add_subplot(223)
growthrateplot = fig.add_subplot(224)
#################################################################################
areaHeight.set_xlabel(r"$A_T$")
areaHeight.set_ylabel(r"$H_T$")

absstressplot.set_xlabel(r"$d_c$")
absstressplot.set_ylabel(r"$\langle \| \sigma_c \| \rangle_t$")

stressRadOrthplot.set_xlabel(r"$d_c$")
stressRadOrthplot.set_ylabel(r"$\langle\sigma_c\rangle_t$")

growthrateplot.set_xlabel(r"$d_c$")
growthrateplot.set_ylabel(r"$\langle\frac{dA_c}{dt}\rangle_t$")

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
	###########################################################
	# Growth Ratio calculation
	############################################################
	if False:#not fastkappaOption:
		print "\n ############################################################"
		print " "*15,"Growth Ratio = ",  sf.getGrowthRatio(numOfLayer = numOfLayer, targetid = targetid,
				endStep = endStep,startStep = startStep)
		print "############################################################"
	#print float(folderdict['n'])
	#print "\n",os.getcwd()
	plotData[etacurrent] = getGrowthRateStress(numOfLayer, endStep=endStep,eta=etacurrent,startStep = startStep,  stepsize = stepsize,
		maxarea = maxarea, areastep = areastep,startarea = startarea,resetids = resetids,
				endarea = endarea)
	############################################################
	plotSurfaceStressGrowthRate(plotData[etacurrent],numOfLayer, eta = etacurrent)
	#print sys.getsizeof(plotData)
	os.chdir("..")
	gc.collect()
	counter+= 1
###############################################################################
#plotting
###############################################################################
plotargs = {"linewidth":5}

for key,data in plotData.iteritems():
	color = scalarMap.to_rgba(key)
	################################## 
	#return 
	#[tissueSurfaceAreaArray,tissueSurfaceAreaArray2, heightArray, timeArray, dictArray]
	#  dictArray = [absStressDict, stressRadialDict, stressOrthoradialDict,growthRateDict, tipDistanceDict]
	##################################
	#area height
	##################################
	areaHeight.plot(data[0], data[2], c=color,**plotargs)
	##################################
	dictArray = data[5]
	absStressDict = dictArray[0]
	stressRadialDict = dictArray[1]
	stressOrthoradialDict = dictArray[2]
	growthRateDict = dictArray[3]
	tipDistanceDict = dictArray[4]
	##################################
	for faceid in absStressDict.keys():
		absstress = np.mean(absStressDict[faceid])
		stressradial = np.mean(stressRadialDict[faceid])
		stressorthoradial = np.mean(stressOrthoradialDict[faceid])
		growthrate = np.nanmean(growthRateDict[faceid])
		tipdistance = np.mean(tipDistanceDict[faceid])
		##################################
		absstressplot.plot(tipdistance, absstress, c= color, marker = 'o',alpha = 0.9,**plotargs)
		stressRadOrthplot.plot(tipdistance, stressradial, c= color, marker = 'o',alpha = 0.9,**plotargs)
		stressRadOrthplot.plot(tipdistance, stressorthoradial, c= color, marker = '<',alpha = 0.9,**plotargs)
		growthrateplot.plot(tipdistance, growthrate,  c= color, marker = 'o',alpha = 0.9,**plotargs)
	##################################
############################################################
# Legend of the plot
############################################################

from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], linestyle = "-", color='k',marker = 'o', ms = 20, label=r"$\sigma_{r}$",**plotargs),
				   Line2D([0], [0], linestyle = "-", color='k',marker = '<', ms = 20, label=r"$\sigma_{o}$",**plotargs),
				   ]
stressRadOrthplot.legend(handles = legend_elements)

"""
# make legend handles for all new plots
growth_legend_elements = [
			Line2D([0], [0], linestyle = "-.", color='k', label=r"$\langle g_{r} \rangle_c $",**plotargs),
			Line2D([0], [0],  color='k', label=r"$\langle g_{o} \rangle_c $",**plotargs),
			Line2D([0], [0], linestyle = "-.", color='k', label=r"$\langle g_{1} \rangle_c $",**plotargs),
			Line2D([0], [0],  color='k',  label=r"$\langle g_{2} \rangle_c $",**plotargs)]
meangrowth12plot.legend(handles = growth_legend_elements[2:])
meangrowthROplot.legend(handles = growth_legend_elements[:2])

stress_legend_elements = [
			Line2D([0], [0], linestyle = "-.", color='k', label=r"$\langle \sigma_{r} \rangle_c $",**plotargs),
			Line2D([0], [0],  color='k', label=r"$\langle \sigma_{o} \rangle_c $",**plotargs),
			Line2D([0], [0], linestyle = "-.", color='k', label=r"$\langle \sigma_{1} \rangle_c $",**plotargs),
			Line2D([0], [0],  color='k',  label=r"$\langle \sigma_{2} \rangle_c $",**plotargs)]
stress12plot.legend(handles = stress_legend_elements[2:])
stressROplot.legend(handles = stress_legend_elements[:2])
"""
###############################################################################
#color bar fig
###############################################################################
fig.tight_layout()
scalarMap._A = []
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.15, 0.07, 0.7, 0.03])
clrbar = plt.colorbar(scalarMap,orientation='horizontal',cax = cbar_ax)

clrbar.set_label(r"Mechanical Feedback, $\eta$")
################################################################################
minarea = min(plotData.values()[0][0])
"""
if False:#startarea:#start area != none
	meangrowth12plot.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
	meangrowthROplot.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
	stressROplot.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
	stress12plot.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
	boundaryareaplot.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
	boundaryareaplot1.set_xticks(np.linspace(startarea,endarea,3).astype('int'))
else:
	meangrowth12plot.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
	meangrowthROplot.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
	stress12plot.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
	stressROplot.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
	boundaryareaplot.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
	boundaryareaplot1.set_xticks(np.linspace(minarea,endarea,3).astype('int'))
"""
################################################################################

if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	fig.savefig(saveDirectory+r"/plot_growthratio_meanstress_meangrowth_targetface=%d.png"%(endStep,targetid),transparent = True, bbox_inches="tight")
	#fig2.savefig(saveDirectory+r"/plot_eta%d_romeangrowthstress_areastep=%d_targetface=%d.eps"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig2.savefig(saveDirectory+r"/plot_growthratio_romeangrowthstress_areastep=%d_targetface=%d.png"%(areastep,targetid),transparent = True, bbox_inches="tight")
	#fig3.savefig(saveDirectory+r"/plot_eta%d_12meangrowthstress_areastep=%d_targetface=%d.eps"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig3.savefig(saveDirectory+r"/plot_growthratio_12meangrowthstress_areastep=%d_targetface=%d.png"%(areastep,targetid),transparent = True, bbox_inches="tight")
	#fig4.savefig(saveDirectory+r"/plot_eta%d_boundaryarea_time=%d_targetface=%d.eps"%(maxeta,endStep,targetid),transparent = True, bbox_inches="tight")
else:
	fig.savefig(saveDirectory+r"/plot_cellgrowthratestress_%s_height_time=%d_area=%d.png"%(zone, endStep,endarea),transparent = True, bbox_inches="tight")
	#fig2.savefig(saveDirectory+r"/plot_scatterheightgrowthrate_%s_height_time=%d_area=%d.png"%(zone, endStep,endarea),transparent = True, bbox_inches="tight")
	#fig3.savefig(saveDirectory+r"/plot%s_mean_radius_height_time=%d_area=%d.png"%(zone, endStep,endarea),transparent = True, bbox_inches="tight")

	#fig2.savefig(saveDirectory+r"/plot_eta%d_romeangrowthstress_areastep=%d_targetface=%d.eps"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig2.savefig(saveDirectory+r"/plot_eta%d_romeangrowthstress_areastep=%d_targetface=%d.png"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig3.savefig(saveDirectory+r"/plot_eta%d_12meangrowthstress_areastep=%d_targetface=%d.eps"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig3.savefig(saveDirectory+r"/plot_eta%d_12meangrowthstress_areastep=%d_targetface=%d.png"%(maxeta,areastep,targetid),transparent = True, bbox_inches="tight")
	#fig4.savefig(saveDirectory+r"/plot_eta%d_boundaryarea_time=%d_targetface=%d.eps"%(maxeta,endStep,targetid),transparent = True, bbox_inches="tight")



#fig1.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_areaPrimodia_%d.png"%endStep,transparent = True)
#fig2.savefig(saveDirectory+r"/plot_eta_vs_curvature_height_surfaceratio_%d.png"%endStep,transparent = True)
#fig3.savefig(saveDirectory+r"/plot_eta_vs_sphericity_%d.png"%endStep,transparent = True)
plt.close('all')
### Saving Data Dictionary ###

### Saving Data Dictionary ###
if fastkappaOption:# if true calculate with respect to changing fastkappa, else Eta
	print " nothing "
else:
	if jobid:
		np.save('job=%d_%s_cellgrowthratestress_time=%d_endarea=%d.npy.npy'%(jobid, zone,endStep,endarea),plotData)
	else:
		np.save('%s_hcellgrowthratestress_time=%d_endarea=%d.npy'%(zone,endStep,endarea),plotData)
################################################################################
print '\n',15*" "+"################################################################################"
print 45*" "+"DONE "
print 15*" "+"################################################################################"

