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
		if sf.checkExternalFace(face):
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
		ratio = abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2))/max(abs(eigenvalue1),abs(eigenvalue2))
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
	maxvalue = 1#maxTrace
	minvalue = -1#minTrace
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
		if sf.checkExternalFace(face):
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
	clrbar.ax.tick_params(labelsize=20)
	clrbar.set_label("Rescaled Magnitude of Strain", fontsize = 30)
	#ax.set_title("Time %d : Magnitude of Strain : Max %.4f   ; Min %.4f"%(step,maxTrace,minTrace), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	if save:
		saveDirectory = "strainFigures_u%02d_v%02d"%(azim,elev)+r"/strainMagnitude"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		plt.savefig(saveDirectory+r"/strainMagnitude_Time=%03d.png"%(step),transparent=True)
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
		if sf.checkExternalFace(face):
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
		ratio = abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2))/max(abs(eigenvalue1),abs(eigenvalue2))
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
		ratio = abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2))/(normMax*max(abs(eigenvalue1),abs(eigenvalue2)))
		if sf.checkExternalFace(face):
			ratio = 0.
		#print "face ID : ", face.getID(), " ratio : ", ratio
		color = scalarMap.to_rgba(ratio)
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax.add_collection3d(pc)
		ax.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
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
		saveDirectory = "strainFigures_u%02d_v%02d"%(azim,elev)+r"/stressSurface"
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
	maxMagnitude = 1#(max(areaArray))
	minMagnitude = 0#(min(areaArray))
	normMax = 12.
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
		color = scalarMap.to_rgba(face.getAreaOfFace()/normMax)
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
	cbar_ax = fig.add_axes([0.873, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar.set_label("Area of Cell", fontsize = 18)
	clrbar.ax.tick_params(labelsize=14)
	#ax.set_title("Time %d : Area of Cell : Max %.4f Min %.4f"%(step,maxMagnitude,minMagnitude), fontsize = 20)
	ax.view_init(azim = azim, elev= elev)
	if save:
		saveDirectory = "strainFigures_u%02d_v%02d"%(azim,elev)+r"/StressfaceArea"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		plt.savefig(saveDirectory+r"/faceArea_Time=%03d.png"%(step),transparent=True)
		plt.close()
	#plt.close("all")
	#return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	return
##########################################################################################
def printFormMatrix(cell):
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face!=None:
		faceid = face.getID()
		print " ========================================================================================"
		print faceid
		print "MU  ", [face.getMu1(),face.getMu2(),face.getMu3(),face.getMu4()]
		print "CF  ",[qd.getCurrentFormMatrix(face,0,0),(qd.getCurrentFormMatrix(face,0,1)),qd.getCurrentFormMatrix(face,1,0),(qd.getCurrentFormMatrix(face,1,1))]
		print "TF  ",[qd.getTargetFormMatrix(face,0,0),(qd.getTargetFormMatrix(face,0,1)),qd.getTargetFormMatrix(face,1,0),(qd.getTargetFormMatrix(face,1,1))]
		face = faces.next()
	return
###############################################################################################################
###         Loading the TargetFormMatrix from file and setting it to the faces for given growthstep
###     For older simulation when xx and yy needs to be switched
###############################################################################################################
def setOldTargetFormMatrix(cell, growthcounter):
	### Loading the TargetFormMatrix
	loadedDictionary = np.load("TargetFormMatrix_step=%d.npy"%growthcounter).item()
	faces = qd.CellFaceIterator(cell)
	face = faces.next()
	while face != None:
		faceid = face.getID()
		matrixArray = loadedDictionary[faceid]
		qd.setTargetFormMatrix(face, 0,0, matrixArray[1][1])
		qd.setTargetFormMatrix(face, 1,1, matrixArray[0][0])
		qd.setTargetFormMatrix(face, 0,1, matrixArray[0][1])
		qd.setTargetFormMatrix(face, 1,0, matrixArray[1][0])
		face = faces.next()
	cell.setParameters()
	######
	return cell
################################################################################################
#			Plotting Part 
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
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
								   default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
								   default = 1., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
								   default = 0.0, type = float)
parser.add_argument("-n","--angle",help = "value to set for convex angle threshold, default = 360",
								   default = 360., type = float)
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
beta = args.beta
zeta  = args.zeta
pressure = args.pressure
numOfLayer = args.layer
gamma = args.gamma
anglethreshold = args.angle
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
	plotStrainMagnitude(cell,numOfLayer,step = step, targetface = targetface,save = True,azim = azim, elev = elev)
	plotStrainSurface(cell,numOfLayer,step = step, save = True,azim = azim, elev = elev)
	plotFaceArea(cell,numOfLayer,step = step, save = True,azim = azim, elev = elev)
	# Directory to save the figures
	saveDirectory = "strainFigures_u%02d_v%02d"%(azim,elev)+r"/surfacePlot"
	if not os.path.exists(saveDirectory):
		os.makedirs(saveDirectory)
	latdev.plotSurface(cell,numOfLayer,name=saveDirectory+r"/surface_%03d.png"%step,alpha = 0.8, azim = azim, elev = elev)
	plt.close("all")
################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

