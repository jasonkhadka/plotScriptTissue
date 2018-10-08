import matplotlib.pyplot as plt
import sys
import numpy as np
# Cluster Paths
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
sys.path.append('/home/jkhadka/transferdata/scripts/strain_plots/')
sys.path.append('/home/jkhadka/plantdev')
sys.path.append('/home/jkhadka/plantdev/python_quadedge')
# Local pc paths
sys.path.append('/Users/jasonkhadka/Documents/git/plantdev/python_quadedge')
sys.path.append("/Users/jasonkhadka/Documents/git/plantdev")
sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scriptsTransferData/scripts/strain_plots/')
sys.path.append('/Users/jasonkhadka/Documents/git/simulations/cluster_simulation/scriptsTransferData/scripts/simulation_functions/')
import quadedge as qd
import Quadedge_lattice_development as latdev
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from simulation_functions import *
import simulation_functions as sf
import ellipse as ep
#%matplotlib inline
import os

import scipy.optimize as sop

#plt.rcParams['figure.figsize'] = (20.0, 10.0)
plt.rcParams['xtick.labelsize'] = 24.
plt.rcParams['ytick.labelsize'] = 24.
plt.rcParams['axes.labelsize'] = 30.
plt.rcParams['legend.fontsize'] = 30.
plt.rcParams['axes.titlesize'] = 30.
plt.rcParams['figure.autolayout']=False




######################################################################
# Function to plot data on the different feedback data
######################################################################
def plotFeedbackPrimordiaHeight(meancurvefile,xlim = None, maxeta=None,
    fitlength = None,plotlen =None,maxarea=None):
    ####################################################################################
    # Plotting 
    ####################################################################################
    curvefiledata = np.load(meancurvefile).item()
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    jet = cm = plt.get_cmap('viridis') 
    if not maxeta:
        maxeta = max(curvefiledata.keys())
    minvalue = min(curvefiledata.keys())
    ##############################################################################
    #print "Max value", maxeta, " minvalue", minvalue
    cNorm  = colors.Normalize(vmin=minvalue, vmax=maxeta)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    #################################################################################
    #     Making the Plots
    #################################################################################
    fig = plt.figure(1,figsize=(20,20))
    #fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    #ax5 = fig.add_subplot(325)
    #ax6 = fig.add_subplot(326)
    #fig.set_aspect(aspect='equal', adjustable='box')
    #ax.axis('equal')
    #################################################################################
    # Min Gaussian curvature
    ##################################
    """
    ax1.set_title("Mean Gaussian Curvature")
    ax1.set_xlabel("Total Surface Area")
    ax1.set_ylabel("Mean Gaussian Curvature")
    """
    ax1.set_title("Primordia height to Area")
    ax1.set_xlabel("Primordia Area")
    ax1.set_ylabel("Primordia Height")
    #ax1.set_ylim(-0.2,0.)
    ###################################
    # Height of Primodia
    ###################################
    ax2.set_title("Height of Primodia vs surface area")
    ax2.set_xlabel("Total Surface Area")
    ax2.set_ylabel("Height of Primodia")
    ###################################
    # Primordial area
    ###################################
    ax3.set_title("Primordia area  vs surface area ")
    ax3.set_xlabel(r"Surface Area of Tissue")
    ax3.set_ylabel(r"Primordia Area")
    ###################################
    # Surface area vs time
    ###################################
    ax4.set_title("Surface area vs time ")
    ax4.set_xlabel(r"$t$")
    ax4.set_ylabel("Total Surface Area")
    """
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax5.set_title("Area vs time ")
    ax5.set_xlabel(r"$t$")
    ax5.set_ylabel(r"$Area$")
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax6.set_title("Area vs scaled time ")
    ax6.set_xlabel(r"$\lambda t$")
    ax6.set_ylabel(r"$Area$")
    """
    ###################################
    ###################################
    #print facefile, curvefile
    curvefiledata = np.load(meancurvefile).item()
    newfilename = meancurvefile.split('/')[-1]
    ########################################################
    facefile = newfilename.split("_")
    facefile = [facefile[1],facefile[2]]
    #print facefile
    filedict = dict(item.split("=") for item in facefile)
    fk = float(filedict['fk'])
    sk = float(filedict['sk'])
    ########################################################
    #print maxarea
    ########################################################
    for eta in sorted(curvefiledata.keys()):
        if eta > maxeta:continue
        curvedata = curvefiledata[eta]
        ###########################################
        heightdata = curvedata[0]
        gaussiancurvedata = curvedata[1]
        primordiaareadata = curvedata[2]
        surfaceareadata = np.array(curvedata[3])
        timedata = curvedata[-1]
        ###########################################
        c = scalarMap.to_rgba(eta)
        ###########################################
        # Plotting
        ###########################################
        if maxarea: 
            #if maxarea is given:
            # plotting till this max area for surfacearea
            plotlen = getPlotlenMaxArea(surfaceareadata,maxarea)
            #print plotlen
        ###########################################
        if plotlen == None:#if no fitlen defined
            # plot of meanGuassianCurvature vs surfacearea
            ax1.plot(primordiaareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # plot of primordial height vs surfacearea
            ax2.plot(surfaceareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # primordia area vs surfacearea
            ax3.plot(surfaceareadata,primordiaareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # surfacearea vs time
            ax4.plot(timedata,surfaceareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
        else:
            # plot of meanGuassianCurvature vs surfacearea
            ax1.plot(primordiaareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # plot of primordial height vs surfacearea
            ax2.plot(surfaceareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # primordia area vs surfacearea
            ax3.plot(surfaceareadata[:plotlen],primordiaareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # surfacearea vs time
            ax4.plot(timedata[:plotlen],surfaceareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk)) 
    ##########################################################################################
    # Making color map
    ##########################################################################################
    if xlim != None:
        ax4.set_xlim(xlim)
        ax6.set_xlim(xlim)
    scalarMap._A = []
    fig.subplots_adjust(bottom=0.12)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.4])
    clrbar = plt.colorbar(scalarMap,orientation='vertical',cax = cbar_ax,ticks=np.linspace(0,maxeta,3))
    #clrbar.ax.set_yticklabels(['Low','Medium','High'],rotation='vertical')
    clrbar.ax.tick_params(labelsize=24) 
    #clrbar.set_label(r"$\eta$",fontsize=24)
    clrbar.set_label(r"Feedback",fontsize=36)
    plt.tight_layout(rect=[0.1,0.1,.8,0.9])
    #plt.savefig(growthtype+"_fk=%.3f_sk=%.3f_feedback.png"%(fk, sk),transparent=True)
    plt.savefig("heightplot_"+"%s.eps"%newfilename[:-4], format='eps',dpi=1200)
    return
######################################################################
# Function to plot data on the different feedback data
######################################################################
def plotFeedbackPrimordiaHeightAgainstSurface(meancurvefile,xlim = None, maxeta=None,
    fitlength = None,plotlen =None,maxarea=None,growthtype = None):
    ####################################################################################
    # Plotting 
    ####################################################################################
    curvefiledata = np.load(meancurvefile).item()
    import matplotlib.colors as colors
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.cm as cmx
    from matplotlib.ticker import FormatStrFormatter
    ############################################################
    jet = cm = plt.get_cmap('plasma') 
    if not maxeta:
        maxeta = max(curvefiledata.keys())
    minvalue = min(curvefiledata.keys())
    ##############################################################################
    #print "Max value", maxeta, " minvalue", minvalue
    cNorm  = colors.Normalize(vmin=minvalue, vmax=maxeta)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    #################################################################################
    #     Making the Plots
    #################################################################################
    fig = plt.figure(1,figsize=(6.5,6.5))
    #fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
    ax1 = fig.add_subplot(111)
    #ax5 = fig.add_subplot(325)
    #ax6 = fig.add_subplot(326)
    #fig.set_aspect(aspect='equal', adjustable='box')
    #ax.axis('equal')
    #################################################################################
    ###################################
    # Height of Primodia
    ###################################
    #ax1.set_title("Height of Primodia vs surface area")
    #if growthtype:
    #    ax1.set_title(growthtype)
    ax1.set_xlabel(r"$A_{T}$")
    ax1.set_ylabel(r"$h$")
    #ax1.set_ylim((0.1,1.0))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ###################################
    """
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax5.set_title("Area vs time ")
    ax5.set_xlabel(r"$t$")
    ax5.set_ylabel(r"$Area$")
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax6.set_title("Area vs scaled time ")
    ax6.set_xlabel(r"$\lambda t$")
    ax6.set_ylabel(r"$Area$")
    """
    ###################################
    ###################################
    #print facefile, curvefile
    curvefiledata = np.load(meancurvefile).item()
    newfilename = meancurvefile.split('/')[-1]
    ########################################################
    facefile = newfilename.split("_")
    facefile = [facefile[1],facefile[2]]
    #print facefile
    filedict = dict(item.split("=") for item in facefile)
    fk = float(filedict['fk'])
    sk = float(filedict['sk'])
    ########################################################
    #print maxarea
    ########################################################
    for eta in sorted(curvefiledata.keys()):
        if eta > maxeta:continue
        curvedata = curvefiledata[eta]
        ###########################################
        heightdata = curvedata[0]
        gaussiancurvedata = curvedata[1]
        primordiaareadata = curvedata[2]
        surfaceareadata = np.array(curvedata[3])
        timedata = curvedata[-1]
        ###########################################
        c = scalarMap.to_rgba(eta)
        ###########################################
        # Plotting
        ###########################################
        if maxarea: 
            #if maxarea is given:
            # plotting till this max area for surfacearea
            plotlen = getPlotlenMaxArea(surfaceareadata,maxarea)
            #print plotlen
        ###########################################
        if plotlen == None:#if no fitlen defined
            # plot of meanGuassianCurvature vs surfacearea
            #ax1.plot(primordiaareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # plot of primordial height vs surfacearea
            ax1.plot(surfaceareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # primordia area vs surfacearea
            #ax3.plot(surfaceareadata,primordiaareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # surfacearea vs time
            #ax4.plot(timedata,surfaceareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
        else:
            # plot of meanGuassianCurvature vs surfacearea
            #ax1.plot(primordiaareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # plot of primordial height vs surfacearea
            ax1.plot(surfaceareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # primordia area vs surfacearea
            #ax3.plot(surfaceareadata[:plotlen],primordiaareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            # surfacearea vs time
            #ax4.plot(timedata[:plotlen],surfaceareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk)) 
    ##########################################################################################
    # Making color map
    ##########################################################################################
    #if xlim != None:
    #    ax4.set_xlim(xlim)
    #    ax6.set_xlim(xlim)
    scalarMap._A = []
    #fig.subplots_adjust(right=0.9,bottom=0.15,left=0.15)
    #print fig.get_axes()
    #cax = fig.add_axes([0.92, 0.2, 0.05, 0.75])
    #divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #clrbar = plt.colorbar(scalarMap,orientation='vertical',cax = cax,ticks=np.linspace(0,maxeta,3),format = '%3s')
    #clrbar.ax.set_yticklabels(['Low','Medium','High'],rotation='vertical')
    #clrbar.ax.tick_params(labelsize=plt.rcParams['xtick.labelsize']) 
    axpos = ax1.get_position()
    #print axpos
    plt.tight_layout(rect=[0.,0.,.9,.9])
    clrbarpos = [axpos.x0+axpos.width,axpos.y0,0.04,axpos.height]
    #cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.62])
    cbar_ax = fig.add_axes(clrbarpos)
    clrbar = plt.colorbar(scalarMap,cax = cbar_ax,ticks=np.linspace(minvalue,maxeta,3))
    clrbar.set_label(r"$\eta$")
    #clrbar.set_label(r"Feedback",fontsize=plt.rcParams['legend.fontsize'] ,labelpad = 0)
    
    #plt.savefig(growthtype+"_fk=%.3f_sk=%.3f_feedback.png"%(fk, sk),transparent=True)
    if growthtype:
        plt.savefig("%s_height_vs_tissueArea_plot.eps"%growthtype, format='eps',dpi=300, bbox_inches="tight")
    else:
        plt.savefig("height_vs_tissueArea_plot_"+"%s.eps"%newfilename[:-4], format='eps',dpi=300)
    return

######################################################################
# Function to get plotlen with max area
######################################################################

def getPlotlenMaxArea(surfaceareadata,maxarea):
    surfaceareadata = np.array(surfaceareadata)
    try:
        plotlen = np.min(np.nonzero(surfaceareadata>maxarea))
    except ValueError:
        plotlen = getPlotlenMaxArea(surfaceareadata,maxarea-10.)
    return plotlen

######################################################################
# Function plot comparative plots
######################################################################
def plotComparativeFeedbackPrimordiaHeight(meancurvefileA,meancurvefileB,xlim = None,
 maxeta=None,fitlength = None,plotlen =None,maxarea=None,save=False):
    ####################################################################################
    # Plotting 
    ####################################################################################
    files = [meancurvefileA,meancurvefileB]
    ####################################################################################
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    ##############################################################################
    #print "Max value", maxeta, " minvalue", minvalue
    curvefiledataA = np.load(meancurvefileA).item()
    curvefiledataB = np.load(meancurvefileB).item()
    if not maxeta:
        maxetaA = max(curvefiledataA.keys())
        maxetaB = max(curvefiledataB.keys())
        maxeta = max(maxetaA,maxetaB)
    else:
        maxetaA = maxeta
        maxetaB = maxeta
    minvalueA = min(curvefiledataA.keys())
    jetA = cm = plt.get_cmap('YlOrRd') 
    cNormA  = colors.Normalize(vmin=minvalueA, vmax=maxetaA)
    scalarMapA = cmx.ScalarMappable(norm=cNormA, cmap=jetA)
    ##############################################################################
    minvalueB = min(curvefiledataB.keys())
    jetB = cm = plt.get_cmap('PuBu') 
    cNormB  = colors.Normalize(vmin=minvalueB, vmax=maxetaB)
    scalarMapB = cmx.ScalarMappable(norm=cNormB, cmap=jetB)
    colormaps = [scalarMapA, scalarMapB]
    #################################################################################
    #     Making the Plots
    #################################################################################
    fig = plt.figure(1,figsize=(20,20))
    #fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    #ax5 = fig.add_subplot(325)
    #ax6 = fig.add_subplot(326)
    #fig.set_aspect(aspect='equal', adjustable='box')
    #ax.axis('equal')
    #################################################################################
    # Min Gaussian curvature
    ##################################
    ax1.set_title("Primordia height to Area")
    ax1.set_xlabel("Primordia Area")
    ax1.set_ylabel("Primordia Height")
    #ax1.set_ylim(-0.2,0.)
    ###################################
    # Height of Primodia
    ###################################
    ax2.set_title("Height of Primodia vs surface area")
    ax2.set_xlabel("Total Surface Area")
    ax2.set_ylabel("Height of Primodia")
    ###################################
    # Primordial area
    ###################################
    ax3.set_title("Primordia area  vs surface area ")
    ax3.set_xlabel(r"Surface Area of Tissue")
    ax3.set_ylabel(r"Primordia Area")
    ###################################
    # Surface area vs time
    ###################################
    ax4.set_title("Surface area vs time ")
    ax4.set_xlabel(r"$t$")
    ax4.set_ylabel("Total Surface Area")
    """
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax5.set_title("Area vs time ")
    ax5.set_xlabel(r"$t$")
    ax5.set_ylabel(r"$Area$")
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax6.set_title("Area vs scaled time ")
    ax6.set_xlabel(r"$\lambda t$")
    ax6.set_ylabel(r"$Area$")
    """
    ###################################
    # Starting the calculation
    ###################################
    for meancurvefile,scalarMap in zip(files, colormaps):
        curvefiledata = np.load(meancurvefile).item()
        newfilename = meancurvefile.split('/')[-1]
        ########################################################
        facefile = newfilename.split("_")
        facefile = [facefile[1],facefile[2]]
        #print facefile
        filedict = dict(item.split("=") for item in facefile)
        fk = float(filedict['fk'])
        sk = float(filedict['sk'])
        ########################################################
        #print maxarea
        ########################################################
        for eta in sorted(curvefiledata.keys()):
            if eta > maxeta:continue
            curvedata = curvefiledata[eta]
            ###########################################
            heightdata = curvedata[0]
            gaussiancurvedata = curvedata[1]
            primordiaareadata = curvedata[2]
            surfaceareadata = np.array(curvedata[3])
            timedata = curvedata[-1]
            ###########################################
            c = scalarMap.to_rgba(eta)
            ###########################################
            # Plotting
            ###########################################
            if maxarea: 
                #if maxarea is given:
                # plotting till this max area for surfacearea
                plotlen = getPlotlenMaxArea(surfaceareadata,maxarea)
                #print plotlen

            ###########################################
            if plotlen == None:#if no fitlen defined
                # plot of meanGuassianCurvature vs surfacearea
                ax1.plot(primordiaareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # plot of primordial height vs surfacearea
                ax2.plot(surfaceareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # primordia area vs surfacearea
                ax3.plot(surfaceareadata,primordiaareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # surfacearea vs time
                ax4.plot(timedata,surfaceareadata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            else:
                # plot of meanGuassianCurvature vs surfacearea
                ax1.plot(primordiaareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # plot of primordial height vs surfacearea
                ax2.plot(surfaceareadata[:plotlen],heightdata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # primordia area vs surfacearea
                ax3.plot(surfaceareadata[:plotlen],primordiaareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # surfacearea vs time
                ax4.plot(timedata[:plotlen],surfaceareadata[:plotlen],c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk)) 
    ##########################################################################################
    # Making color map
    ##########################################################################################
    simNameA  = meancurvefileA.split('/')[-1][:-4]#.split('_')[0]
    simNameB = meancurvefileB.split('/')[-1][:-4]#.split('_')[0]
    if xlim != None:
        ax4.set_xlim(xlim)
        ax6.set_xlim(xlim)
    scalarMapA._A = []
    scalarMapB._A = []
    #fig.subplots_adjust(bottom=0.12)
    ##########################################################################################
    cbar_ax = fig.add_axes([0.85, 0.01, 0.02, 0.4])
    clrbar = plt.colorbar(scalarMapA,orientation='vertical',cax = cbar_ax,ticks=np.linspace(0,maxetaA,4))
    #clrbar.ax.set_yticklabels(['Low','Medium','High'],rotation='vertical')
    clrbar.ax.tick_params(labelsize=20) 
    #clrbar.set_label(r"$\eta$",fontsize=24)
    clrbar.set_label(r"%s"%simNameA,fontsize=24)
    ##########################################################################################
    cbar_ax = fig.add_axes([0.85, 0.51, 0.02, 0.4])
    clrbar = plt.colorbar(scalarMapB,orientation='vertical',cax = cbar_ax,ticks=np.linspace(0,maxetaB,4))
    #clrbar.ax.set_yticklabels(['Low','Medium','High'],rotation='vertical')
    clrbar.ax.tick_params(labelsize=20) 
    #clrbar.set_label(r"$\eta$",fontsize=24)
    clrbar.set_label(r"%s"%simNameB,fontsize=24)
    ##########################################################################################
    plt.tight_layout(rect=[0.,0.,.8,0.9])
    #plt.savefig(growthtype+"_fk=%.3f_sk=%.3f_feedback.png"%(fk, sk),transparent=True)
    if save:
        plt.savefig("heightplot_"+"%s.eps"%newfilename[:-4], format='eps',dpi=400)
    return
######################################################################
# Function plot comparative plots
######################################################################
def plotComparativeNoFeedbackPrimordiaHeight(meancurvefileA,meancurvefileB,xlim = None, maxeta=None,
    fitlength = None,plotlen =None,maxarea=None,filename = None):
    ####################################################################################
    # Plotting 
    ####################################################################################
    files = [meancurvefileA,meancurvefileB]
    linestyles = [ '--', ':']
    colorlist = ['g','k']
    lw = 3
    newfilename = meancurvefileA.split('/')[-1][:-4].split('_')[0]
    simName = newfilename[0].upper()+newfilename[1:]+" growth"
    ####################################################################################
    curvefiledata = np.load(meancurvefileA).item()
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    if not maxeta:
        maxeta = max(curvefiledata.keys())
    minvalue = min(curvefiledata.keys())
    #################################################################################
    #     Making the Plots
    #################################################################################
    fig = plt.figure(1,figsize=(15,15))
    #fig.suptitle("Time Step = %d"%endStep,fontsize = 40)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    #ax5 = fig.add_subplot(325)
    #ax6 = fig.add_subplot(326)
    #fig.set_aspect(aspect='equal', adjustable='box')
    #ax.axis('equal')
    #################################################################################
    # Min Gaussian curvature
    ##################################
    """
    ax1.set_title(simName)
    ax1.set_xlabel("Total Surface Area")
    ax1.set_ylabel("Mean Gaussian Curvature")
    """
    ax1.set_title("Primordia height to Area")
    ax1.set_xlabel("Primordia Area")
    ax1.set_ylabel("Primordia Height")
    #ax1.set_ylim(-0.2,0.)
    ###################################
    # Height of Primodia
    ###################################
    ax2.set_title(simName)
    ax2.set_xlabel("Total Surface Area")
    ax2.set_ylabel("Height of Primodia")
    ###################################
    # Primordial area
    ###################################
    ax3.set_title(simName)
    ax3.set_xlabel(r"Surface Area of Tissue")
    ax3.set_ylabel(r"Primordia Area")
    ###################################
    # Surface area vs time
    ###################################
    ax4.set_title("Surface area vs time ")
    ax4.set_xlabel(r"$t$")
    ax4.set_ylabel("Total Surface Area")
    """
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax5.set_title("Area vs time ")
    ax5.set_xlabel(r"$t$")
    ax5.set_ylabel(r"$Area$")
    ###################################
    # Mean Area of Cells non scaled
    ###################################
    ax6.set_title("Area vs scaled time ")
    ax6.set_xlabel(r"$\lambda t$")
    ax6.set_ylabel(r"$Area$")
    """
    ###################################
    # Starting the calculation
    ###################################
    for meancurvefile,c,linestyle in zip(files, colorlist,linestyles):
        curvefiledata = np.load(meancurvefile).item()
        newfilename = meancurvefile.split('/')[-1]
        ########################################################
        facefile = newfilename.split("_")
        facefile = [facefile[1],facefile[2]]
        #print facefile
        filedict = dict(item.split("=") for item in facefile)
        fk = float(filedict['fk'])
        sk = float(filedict['sk'])
        ########################################################
        #print maxarea
        ########################################################
        for eta in sorted(curvefiledata.keys()):
            if eta != 0.:continue
            curvedata = curvefiledata[eta]
            ###########################################
            heightdata = curvedata[0]
            gaussiancurvedata = curvedata[1]
            primordiaareadata = curvedata[2]
            surfaceareadata = np.array(curvedata[3])
            timedata = curvedata[-1]
            ###########################################
            # Plotting
            ###########################################
            if maxarea: 
                #if maxarea is given:
                # plotting till this max area for surfacearea
                plotlen = getPlotlenMaxArea(surfaceareadata,maxarea)
                #print plotlen
            ###########################################
            if plotlen == None:#if no fitlen defined
                # plot of meanGuassianCurvature vs surfacearea
                ax1.plot(primordiaareadata,heightdata,c=c,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # plot of primordial height vs surfacearea
                ax2.plot(surfaceareadata,heightdata,c=c,linestyle = linestyle,lw = lw, label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # primordia area vs surfacearea
                ax3.plot(surfaceareadata,primordiaareadata,c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # surfacearea vs time
                ax4.plot(timedata,surfaceareadata,c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
            else:
                # plot of meanGuassianCurvature vs surfacearea
                ax1.plot(primordiaareadata[:plotlen],heightdata[:plotlen],c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # plot of primordial height vs surfacearea
                ax2.plot(surfaceareadata[:plotlen],heightdata[:plotlen],c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # primordia area vs surfacearea
                ax3.plot(surfaceareadata[:plotlen],primordiaareadata[:plotlen],c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk))
                # surfacearea vs time
                ax4.plot(timedata[:plotlen],surfaceareadata[:plotlen],c=c,linestyle = linestyle, lw = lw,label = r'$\kappa_{fast}=%.3f; \kappa_{slow}=%.3f$'%(fk,sk)) 
            ###################################################################################################
            break
    ##########################################################################################
    # Making color map
    ##########################################################################################
    if xlim != None:
        ax4.set_xlim(xlim)
        ax6.set_xlim(xlim)
    ###############################################################
    # Legend
    ###############################################################
    ax2.legend()
    ###############################################################
    #plt.tight_layout(rect=[0.1,0.1,.8,0.9])
    plt.tight_layout()
    ###############################################################
    # Save fig
    ###############################################################
    # Save just the portion _inside_ the second axis's boundaries
    extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # Pad the saved area by 10% in the x-direction and 20% in the y-direction
    fig.savefig("height_vs_area_"+"%s_growthratecomparision.eps"%newfilename[:-4], 
        format='eps',dpi=1200, bbox_inches=extent.expanded(1.38, 1.26))
    #plt.savefig(growthtype+"_fk=%.3f_sk=%.3f_feedback.png"%(fk, sk),transparent=True)
    if filename == None:
        plt.savefig("heightplot_"+"%s.eps"%newfilename[:-4], format='eps',dpi=1200)
    else:
        plt.savefig(filename, format='eps',dpi=1200)
    return