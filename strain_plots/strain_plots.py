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



#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 24.
plt.rcParams['ytick.labelsize'] = 24.
plt.rcParams['axes.labelsize'] = 24.
plt.rcParams['legend.fontsize'] = 24.
plt.rcParams['axes.titlesize'] = 30




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
	 fig = plt.figure(frameon=False,figsize=(10,10))
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
		  """if sf.checkExternalFace(face):
				face  = faces.next()
				continue
		  """
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
		  """if sf.checkExternalFace(face):
				ratio = 0.
		  """
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
# Stress difference
##########################################################################################
def getDifference(lambda1, lambda2):
	lambdamax = max(lambda1, lambda2)
	lambdamin = min(lambda1, lambda2)
	return (lambdamax-lambdamin)

##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotStrainDifference(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False,azim = -70, elev=50,
					  colormap = 'inferno'):
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
	fig = plt.figure(frameon=False,figsize=(10,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.5*radius,0.9*radius))
	ax.set_ylim((-0.5*radius,0.9*radius))
	ax.set_zlim((-0.5*radius,1.*radius))
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
		eigenvec1 = face.getStressEigenVector1()
		eigenvec2 = face.getStressEigenVector2()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = getDifference(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
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
	jet = cm = plt.get_cmap(colormap) 
	maxvalue = maxEigenValueRatio
	minvalue = minEigenValueRatio
	normMax = 1 # value to normalize by
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
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = getDifference(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		eigenvalueratioarray.append(ratio)
		if sf.checkExternalFace(face):
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
	cbar_ax = fig.add_axes([0.78, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax,ticks=np.linspace(minvalue,maxvalue,5))#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label(r"Stress difference $(s_d)$")
	#clrbar.ax.tick_params(labelsize=24)
	#ax.set_title("Time %d : Anisotrophy of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	ax.view_init(azim = azim, elev= elev)
	if save:
		saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressSurface"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		#plt.tight_layout()
		#plt.savefig(saveDirectory+r"/strainSurface_Time=%03d.eps"%(step),dpi = 500,format='eps', transparent=True)
		plt.savefig(saveDirectory+r"/stressDifference_Time=%03d.png"%(step),dpi = 500,format='png', transparent=True)
		plt.close()
	return# eigenvalueratioarray
##########################################################################################
# Anistropic measure
##########################################################################################
def getAnistropy(lambda1, lambda2):
	lambda1 = abs(lambda1)
	lambda2 = abs(lambda2)
	lambdamax = max(lambda1, lambda2)
	lambdamin = min(lambda1, lambda2)
	return (lambdamax-lambdamin)/(lambdamax+lambdamin)

##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################

def plotStrainAnisotropy(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False,azim = -70, elev=50,
					  colormap = 'plasma'):
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
	fig = plt.figure(frameon=False,figsize=(10,10))
	fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax = Axes3D(fig)
	ax.set_xlim((-0.5*radius,0.9*radius))
	ax.set_ylim((-0.5*radius,0.9*radius))
	ax.set_zlim((-0.5*radius,1.*radius))
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
		eigenvec1 = face.getStressEigenVector1()
		eigenvec2 = face.getStressEigenVector2()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = getAnistropy(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
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
	jet = cm = plt.get_cmap(colormap) 
	maxvalue = 1.#maxEigenValueRatio
	minvalue = 0. #minEigenValueRatio
	normMax = 1 # value to normalize by
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
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(max(eigen1,eigen2))
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = getAnistropy(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		eigenvalueratioarray.append(ratio)
		if sf.checkExternalFace(face):
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
	cbar_ax = fig.add_axes([0.78, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	clrbar = plt.colorbar(scalarMap,cax = cbar_ax,ticks=np.linspace(minvalue,maxvalue,5))#,orientation='horizontal',cax = cbar_ax)
	clrbar.set_label(r"Stress Anisotrophy $(s_a)$")
	#ax.set_title("Time %d : Anisotrophy of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	ax.view_init(azim = azim, elev= elev)
	if save:
		saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressSurface"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		#plt.tight_layout()
		#plt.savefig(saveDirectory+r"/strainSurface_Time=%03d.eps"%(step),dpi = 500,format='eps', transparent=True)
		plt.savefig(saveDirectory+r"/stressAnisotropy_Time=%03d.png"%(step),dpi = 500,format='png', transparent=True)
		plt.close()
	return# eigenvalueratioarray


def getAbsAnistropy(lambda1, lambda2):
	lambda1 = abs(lambda1)
	lambda2 = abs(lambda2)
	lambdamax = max(lambda1, lambda2)
	lambdamin = min(lambda1, lambda2)
	return (lambdamax-lambdamin)/(lambdamax+lambdamin)
################################################################
def getDifference(lambda1, lambda2):
	lambdamax = max(lambda1, lambda2)
	lambdamin = min(lambda1, lambda2)
	return (lambdamax-lambdamin)#/(lambdamax+lambdamin)
################################################################
def plotAbsAnisotropyStress(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False,azim = -60, elev=50,
                      colormap = 'magma'):
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
	fig = plt.figure(frameon=True,figsize=(22,10))
	#fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	ax1 = fig.add_subplot(121,projection='3d')
	ax2 = fig.add_subplot(122,projection='3d')
	#ax = Axes3D(fig)
	xlim = 0.5
	ax1.set_xlim((-xlim*radius,xlim*radius))
	ax1.set_ylim((-xlim*radius,xlim*radius))
	ax1.set_zlim((-0.3*radius,0.7*radius))
	ax1.axis('off')
	ax1.xaxis.pane.set_edgecolor('black')
	ax1.yaxis.pane.set_edgecolor('black')
	#################################################
	ax2.set_xlim((-xlim*radius,xlim*radius))
	ax2.set_ylim((-xlim*radius,xlim*radius))
	ax2.set_zlim((-0.3*radius,.7*radius))
	ax2.axis('off')
	ax2.xaxis.pane.set_edgecolor('black')
	ax2.yaxis.pane.set_edgecolor('black')
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
		eigenvec1 = face.getStressEigenVector1()
		eigenvec2 = face.getStressEigenVector2()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		ratio = getAbsAnistropy(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
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
	############################################
	#                 Plotting the Cell        #
	############################################
	######### Color Map A ######################
	jet1 = cm1 = plt.get_cmap(colormap) 
	maxvalue = 1.
	minvalue = 0.
	normMax = 1 # value to normalize by
	cNorm1  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap1 = cmx.ScalarMappable(norm=cNorm1, cmap=jet1)
	######### Color Map B ######################
	jet2 = cm2 = plt.get_cmap(colormap) 
	maxvalue = maxEigenValueRatio
	minvalue = minEigenValueRatio
	normMax = 1 # value to normalize by
	cNorm2  = colors.Normalize(vmin=minvalue, vmax=maxvalue)
	scalarMap2 = cmx.ScalarMappable(norm=cNorm2, cmap=jet2)
	###################
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
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(max(eigen1,eigen2))
		ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		eigenvec1 = face.getStressEigenVector1()
		eigenvec2 = face.getStressEigenVector2()
		eigenvalue1 = face.getStressEigenValue1()
		eigenvalue2 = face.getStressEigenValue2()
		eigenvalue.append(eigenvalue1)
		eigenvalue.append(eigenvalue2) 
		eigenvalue1array.append(eigenvalue1)
		eigenvalue2array.append(eigenvalue2)
		eigen1 = eigenvalue1
		eigen2 = eigenvalue2
		#ratio = (max(eigen1,eigen2)- min(eigen1,eigen2))/max(eigen1,eigen2)
		absratio = getAbsAnistropy(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		diffratio = getDifference(eigen1,eigen2)#(max(eigen1,eigen2)- min(eigen1,eigen2))/(eigen1+eigen2)
		#####################################################################
		eigenvalueratioarray.append(ratio)
		if sf.checkExternalFace(face):
			absratio = 0.
			diffratio = 0.
		#print "face ID : ", face.getID(), " ratio : ", ratio
		color1 = scalarMap1.to_rgba(absratio)
		color2 = scalarMap2.to_rgba(diffratio)
		#print face.getID(), ratio
		#print face.getZCentralised(), alpha_fac
		#ax.add_collection3d(arrow(xcen-0.5,ycen-0.5,zcen-0.5,xcen+0.5,ycen+0.5,zcen+0.5))
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color1,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax1.add_collection3d(pc)
		ax1.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
		###########################################################################################
		pc = Poly3DCollection(verts,alpha = alpha,facecolor = color2,linewidths=1,zorder=0)
		pc.set_edgecolor('k')
		ax2.add_collection3d(pc)
		ax2.scatter(xcenarray[-1], ycenarray[-1],zcenarray[-1],c='r')
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
		plot1 = ax1.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=colorvec,length = veclength,pivot='tail',zorder=-1)
		plot2 = ax2.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color=colorvec,length = veclength,pivot='tail',zorder=-1)
		#ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
	scalarMap1._A = []
	scalarMap2._A = []
	#cbar_ax = fig.add_axes([0.78, 0.2, 0.04, 0.55])
	#clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	#clrbar = plt.colorbar(scalarMap,cax = cbar_ax,ticks=np.linspace(minvalue,maxvalue,5))#,orientation='horizontal',cax = cbar_ax)
	clrbar1 = plt.colorbar(scalarMap1, ax=ax1,shrink = 0.5,aspect = 10,ticks=np.linspace(0,1,5))#,orientation='horizontal',cax = cbar_ax)
	clrbar1.set_label(r"Stress Anisotrophy $(s_a)$")
	clrbar2 = plt.colorbar(scalarMap2, ax=ax2,shrink = 0.5,aspect = 10, ticks=np.linspace(minvalue,maxvalue,5))#,orientation='horizontal',cax = cbar_ax)
	clrbar2.set_label(r"Stress Difference $(s_d)$")
	#ax.set_title("Time %d : Anisotrophy of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
	#print xcenarray-0.5
	#plt.close("all")
	ax1.view_init(azim = azim, elev= elev)
	ax2.view_init(azim = azim, elev= elev)
	if save:
		saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/stressSurface"
		import os
		if not os.path.exists(saveDirectory):
			os.makedirs(saveDirectory)
		#plt.tight_layout()
		#plt.savefig(saveDirectory+r"/strainSurface_Time=%03d.eps"%(step),dpi = 500,format='eps', transparent=True)
		plt.savefig(saveDirectory+r"/stressAnisotropy_Time=%03d.png"%(step),dpi = 500,format='png', transparent=True)
		plt.close()
	return# eigenvalueratioarray




##########################################################################################
#       Function to Plot Magnitude of Normal Forces on the Faces of the Cell
##########################################################################################
def plotFaceArea(cell, numOfLayer, step = None, alpha = 0.8, Length=1.0,save = False,azim = -70, elev=50,
				colormap = 'viridis'):
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
	 fig = plt.figure(frameon=False,figsize=(10,10))
	 fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
	 ax = Axes3D(fig)
	 ax.set_xlim((-0.5*radius,0.9*radius))
	 ax.set_ylim((-0.5*radius,0.9*radius))
	 ax.set_zlim((-0.5*radius,1.*radius))
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
	 jet = cm = plt.get_cmap(colormap) 
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
	 cbar_ax = fig.add_axes([0.78, 0.2, 0.04, 0.55])
	 #clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
	 clrbar = plt.colorbar(scalarMap,cax = cbar_ax)#,orientation='horizontal',cax = cbar_ax)
	 clrbar.set_label(r"$A_c$")
	 #clrbar = plt.colorbar(scalarMap, format=ticker.FuncFormatter(fmt))
	 #clrbar.set_label("Area of Cell")
	 #ax.set_title("Time %d : Area of Cell : Max %.4f Min %.4f"%(step,maxMagnitude,minMagnitude), fontsize = 20)
	 ax.view_init(azim = azim, elev= elev)
	 if save:
		  saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)+r"/faceArea"
		  import os
		  if not os.path.exists(saveDirectory):
				os.makedirs(saveDirectory)
		  plt.savefig(saveDirectory+r"/faceArea_Time=%03d.png"%(step),dpi = 500,format='png', transparent=True)
		  plt.close()
	 #plt.close("all")
	 #return eigenvalueratioarray, eigenvalue1array, eigenvalue2array
	 return
################################################################################################
#        Plotting Part 
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
											  default = .8, type = float)
parser.add_argument("-n","--nonNormalize", help = "if option is used, the figures are not normalised", action= "store_false")
parser.add_argument("-b","--beta",help = "value to set for Beta (weigth of second term of Energy), default = 0",
											  default = 0., type = float)
parser.add_argument("-z","--zeta",help = "value to set for Zeta (weigth of fourth term of Energy), default = 0.",
											  default = 0., type = float)
parser.add_argument("-p","--pressure",help = "value to set for Pressure (wegith of third term of Energy), default = 0.001",
											  default = 0.0, type = float)
parser.add_argument("-o","--angle",help = "value to set for convex angle threshold, default = 360",
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
norm = args.nonNormalize

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
	cell.calculateStressStrain()
	#print " Energy: ", cell.getEnergyCartesianVolume()
	###########################################################
	#plotting and Savnig
	############################################################
	#plotStrainMagnitude(cell,numOfLayer,step = step, alpha = alpha, targetface = targetface,save = True,azim = azim, elev = elev)
	plotAbsAnisotropyStress(cell,numOfLayer,step = step, alpha = alpha,  save = True,colormap='inferno',azim = azim, elev = elev)
	#plotStrainAnisotropy(cell,numOfLayer,step = step, alpha = alpha,  save = True,azim = azim, elev = elev)
	#plotStrainDifference(cell,numOfLayer,step = step, alpha = alpha,  save = True,azim = azim, elev = elev)
	plotFaceArea(cell,numOfLayer,step = step, save = True,azim = azim, elev = elev)
	### Plotting surface
	saveDirectory = "strainFigures_u=%03d_v=%03d"%(azim,elev)
	surfaceSaveDirectory = saveDirectory+r"/surfaceFigures"
	if not os.path.exists(surfaceSaveDirectory):
			os.makedirs(surfaceSaveDirectory)
	latdev.plotSurface(cell,numOfLayer, name = surfaceSaveDirectory+r"/surface_u=%02d_v=%02d_plot_=%03d.png"%(azim,elev,step), alpha = alpha,  azim = azim, elev = elev)
	plt.close('all')


################################################################################
print "\n ################################################################################"
print " DONE "
print "################################################################################"

