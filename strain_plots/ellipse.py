import sys
sys.path.append("/Users/jasonkhadka/Documents/git/plantdev")
import quadedge as qd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

###############################################################
# Function to get the FORM Matrix of given cells
###############################################################
def getFormMatrix(cell,targetfaceid):
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face!=None:
        faceid = face.getID()
        if faceid != targetfaceid:
            face= faces.next()
            continue
        targetForm = np.array([[qd.getTargetFormMatrix(face,0,0),qd.getTargetFormMatrix(face,0,1)],
                               [qd.getTargetFormMatrix(face,1,0), qd.getTargetFormMatrix(face,1,1)]])
        currentForm = np.array([[qd.getCurrentFormMatrix(face,0,0),qd.getCurrentFormMatrix(face,0,1)],
                               [qd.getCurrentFormMatrix(face,1,0), qd.getCurrentFormMatrix(face,1,1)]])
        break
    return targetForm, currentForm
###############################################################
def getIntrinsicCoordinates(cell,targetfaceid):
    faces = qd.CellFaceIterator(cell)
    face = faces.next()
    while face!=None:
        faceid = face.getID()
        if faceid != targetfaceid:
            face= faces.next()
            continue
        edges = qd.FaceEdgeIterator(face)
        edge = edges.next()
        x = []
        y = []
        while edge != None:
            vertex = edge.Dest()
            ## Getting coordinates
            x.append(vertex.getProjectedXcoordinate(targetfaceid))
            y.append(vertex.getProjectedYcoordinate(targetfaceid))
            edge = edges.next()
        break
    x.append(x[0])
    y.append(y[0])
    return x,y
#import pylab as mp
def plot_eigenvalues(cov, ax,plot_kwargs=None):
    eig_vec,eig_val,u = np.linalg.svd(cov)
    eig_vec = [eig_val[0]*eig_vec[0],eig_val[1]*eig_vec[1]]
    #print [0,eig_vec[0][0]],[0,eig_vec[0][1]], [0,eig_vec[1][0]],[0,eig_vec[1][1]]
    if plot_kwargs is None:
        ax.plot([0,eig_vec[0][0]],[0,eig_vec[0][1]])
        ax.plot([0,eig_vec[1][0]],[0,eig_vec[1][1]])
    else:
        ax.plot([0,eig_vec[0][0]],[0,eig_vec[0][1]],**plot_kwargs)
        ax.plot([0,eig_vec[1][0]],[0,eig_vec[1][1]],**plot_kwargs)
    return ax
def plot_ellipse(semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e3,ax=None,plot_kwargs=None,\
                    fill=False,fill_kwargs=None,data_out=False,cov=None,mass_level=0.68,norm = False):
    # Get Ellipse Properties from cov matrix
    #print " plott ellipse"
    if cov is not None:
        eig_vec,eig_val,u = np.linalg.svd(cov)
        # Make sure 0th eigenvector has positive x-coordinate
        #if eig_vec[0][0] < 0:
        #    eig_vec[0] *= -1
        #semimaj = np.sqrt(eig_val[0])
        #semimin = np.sqrt(eig_val[1])
        #if mass_level is None:
        #    multiplier = np.sqrt(2.279)
        #else:
        #    distances = np.linspace(0,20,20001)
        #    chi2_cdf = chi2.cdf(distances,df=2)
        #    multiplier = np.sqrt(distances[np.where(np.abs(chi2_cdf-mass_level)==np.abs(chi2_cdf-mass_level).min())[0][0]])
        #semimaj *= multiplier
        #semimin *= multiplier
        unitx = np.array([1.,0.])
        dot = np.dot(unitx,eig_vec[0])
        det = unitx[0]*eig_vec[0][1]-unitx[1]*eig_vec[0][0]
        phi = np.fmod(((np.arctan2(det,dot))+(np.pi*2.)),np.pi*2)
        #print "check this : : ", phi
        #phi = np.arccos(np.dot(eig_vec[0],np.array([1,0])))
        #if eig_vec[0][1] < 0 and phi > 0:
        #    phi *= -1
        semimaj = eig_val[0]
        semimin = eig_val[1]
        #print "maj :", semimaj,"min :", semimaj
        #print "maj :", semimaj,"min :", semimin
    # Generate data for ellipse structure
    # Normalising
    if norm:
        normfac = max(semimaj,semimin)
        semimaj /= normfac
        semimin /= normfac
    #print "maj :", semimaj,"min :", semimin
    #print "---------------------------------"
    theta = np.linspace(0,2*np.pi,theta_num)
    r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    data = np.array([x,y])
    S = np.array([[semimaj,0],[0,semimin]])
    R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    T = np.dot(R,S)
    data = np.dot(T,data)
    data[0] += x_cent
    data[1] += y_cent

    # Output data?
    if data_out == True:
        return data

    # Plot!
    return_fig = False
    if ax is None:
        return_fig = True
        fig,ax = plt.subplots()
        #plotting the eigen vectors 
    if plot_kwargs is None:
        ax.plot(data[0],data[1],color='b',linestyle='-')
    else:
        ax.plot(data[0],data[1],**plot_kwargs)

    if fill == True:
        ax.fill(data[0],data[1],**fill_kwargs)

    if return_fig == True:
        return fig
###################################################################################
#### get Target Form Matrix Ellipse Points ####
###################################################################################
def getTargetFormMatrixEllipsePoints(targetface = 10):
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
    ##### getting data from ellipse
    data = plot_ellipse(cov=targetformmatrix, data_out=True)
    ################################
    return data
##########################################################################################
#       Function to Plot STRAIN on the Surface of the Tissue
##########################################################################################
def plotStrainSurface(cell, numOfLayer,step = None, alpha = 0.8, Length=1.0,save=False):
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
        if False :#checkExternalFace(face):
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
    maxEigenValue = max(eigenvalue)
    #print "Max Eigen Value Ratio :", maxEigenValueRatio
    #print "Min Eigen Value Ratio :", minEigenValueRatio
    ########    ########    ########    ########    ########
    #                 Plotting the Cell                    #
    ########    ########    ########    ########    ########
    ######### Color Map
    jet = cm = plt.get_cmap('plasma') 
    maxvalue = maxEigenValueRatio
    minvalue = minEigenValueRatio
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
        ratio = abs(max(eigenvalue1, eigenvalue2)- min(eigenvalue1, eigenvalue2))/max(abs(eigenvalue1),abs(eigenvalue2))
        if False:#checkExternalFace(face):
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
        ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',length = veclength,pivot='tail',zorder=-1)
        #ax.quiver(X[i], Y[i], Z[i], U[i], V[i], W[i],color='k',pivot='tail',zorder=4)# for older matplotlib
    scalarMap._A = []
    clrbar = plt.colorbar(scalarMap,shrink=0.5, aspect=10)
    clrbar.set_label("Magnitude of Strain Anisotrophy")
    ax.set_title("Time %d : Magnitude of Strain : Max %.4f Min %.4f"%(step,maxEigenValueRatio,minEigenValueRatio), fontsize = 20)
    #print xcenarray-0.5
    #plt.close("all")
    #ax.view_init(azim = 90, elev=90)
    ax.view_init(azim = -70, elev=50)
    if save:
        saveDirectory = "strainFigures"
        import os
        if not os.path.exists(saveDirectory):
            os.makedirs(saveDirectory)
        plt.savefig(saveDirectory+r"/strainSurface_Time=%03d.png"%(step))
        plt.close()
    return
################################################################################################
# Function to plot ellipse
# Fixing the orientation of CFM and TFM
################################################################################################
def plot2DEllipse(cell,targetid,othermatrix = None):
    ############
    # plotting #
    ############
    fig = plt.figure(2,figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_title("Corrected Matrix Calculation")
    #ax.axis('equal')
    ax.set_xlim(-2.5,2.5)
    ax.set_ylim(-2.5,2.5)
    ## Initial shape
    cell.setParameters()
    ## Plotting Coordinates of Face ##
    #ax.plot(x,y,label='coordinate %d'%0)
    ## Geting the matrices
    targetForm, currentForm = getFormMatrix(cell,targetid)
    ## Getting coordinates
    x,y = getIntrinsicCoordinates(cell,targetid)
    ## Plotting Coordinates of Face ##
    ax.plot(x,y,label='coordinate')
    ## Getting Ellipse of CFM
    print "CFM : "
    data = plot_ellipse(cov=currentForm,data_out=True)
    ## plot Ellipse
    ax.plot(data[0],data[1],'-.',label='CFM',zorder=10,c='orange')
    ax = plot_eigenvalues(cov=currentForm,ax = ax,plot_kwargs = {'ls':'-.','zorder':1,'color':'orange'})
    ## Getting Ellipse of TFM
    print "TFM : "
    data = plot_ellipse(cov=targetForm,data_out=True)
    ax.plot(data[0],data[1],'k',label='TFM ')
    ax = plot_eigenvalues(cov=targetForm,ax = ax,plot_kwargs = {'ls':'-.','zorder':10,'color':'k'})
    ## Strain matrix
    strain = np.subtract(currentForm, targetForm)/(targetForm[0,0]+targetForm[1,1])
    print "Strain : "
    data = plot_ellipse(cov=strain,data_out=True)
    ax.plot(data[0],data[1],label='strain',c='r')
    ax = plot_eigenvalues(cov=strain,ax = ax,plot_kwargs = {'zorder':9,'color':'r'})
    if othermatrix != None:
        othermatcounter = 1
        for matrix in othermatrix:
            c=np.random.rand(3,)
            data = plot_ellipse(cov=matrix,data_out=True)
            ## plot Ellipse
            ax.plot(data[0],data[1],'-.',label='OtherMatrix %d'%othermatcounter,zorder=10,c=c)
            ax = plot_eigenvalues(cov=matrix,ax = ax,plot_kwargs = {'ls':'-.','zorder':1,'color':c})
    ###################################
    ax.legend(loc='best')
    plt.show()
    return
################################################################################################
# Growing Targeted Face with Feedback Strain
################################################################################################
def targetedFeedbackStrainGrowFaces(cell,targetid):
    #cell.calculateVertexForce()# Vertex Force calculation 
    cell.calculateStressStrain()# Stress-Strain calculation
    faces = qd.CellFaceIterator(cell)
    #performing growth in the face
    face = faces.next()
    while face != None:
        if face.getID() != targetid:
            face = faces.next()
            continue
        face.feedbackStrainGrow()
        break
    #now setting new parameters on the cell
    #cell = settingParameters(cell)
    cell.setParameters()
    return
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