import numpy as np
import pandas as pd
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
import matplotlib.tri as tri
import scipy.stats as st
from scipy.interpolate import griddata



# Binary diagrams: -------------------------------------------------------------------------

def plot_histo(dataset, mask=None, legend=None, xlabel= None, ylabel = None, linestyle=None, linewidth=None):
    import matplotlib.pyplot as plt
    import numpy as np
    histo=[]
    plt.figure()
    fig, ax = plt.subplots()
    for i in range (len(dataset)):
        #dataset[i].data[mask[i].data] = np.nan
        histo.append(dataset[i].get_histogram())
        ax.plot(histo[i].axes_manager[0].axis, (histo[i].data*100/np.sum(histo[i].data)), label=legend[i], linestyle=linestyle[i], linewidth=linewidth[i])

    leg = ax.legend();
    plt.xlim([None, None]) 
    plt.grid()
    plt.xlabel(xlabel, size=14)
    plt.ylabel(ylabel, size=14)
    plt.legend()


    
# Ternary diagrams: ------------------------------------------------------------------------

def triplot(dataset, mask=None, legend= None, type=None, markersize = 1.2, color=None):
    """Plots a ternary diagram with dots : Each dot represent a pixel of the given data"""
    tri_y=[]
    tri_x=[]
    for i in range (len(dataset)):
        tri_y.append(dataset[i][0]/(dataset[i][0]+dataset[i][1]+dataset[i][2]))
        tri_x.append(np.array(-0.5)*(np.array(1.)-tri_y[i])+(dataset[i][2]/(dataset[i][0]+dataset[i][1]+dataset[i][2])))
    
    #### Plot initialisation: --------------------------------------------------------------
    plt.figure()
    plt.plot([-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0],  [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], color = 'black', marker="_")
    plt.plot([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],  [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0], color='black', marker="_")
    plt.plot([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], [0,0,0,0,0,0,0,0,0,0,0], color='black', marker="|")
    
    if type=='silicate':
        plt.plot([0.3, -0.3], [0.4, 0.4], linestyle='--', color='black')
        plt.annotate('Serpentine', xy=(0.3, 0.4), xytext = (0.33, 0.4), size=10) 
        plt.plot([0.21, -0.21], [0.57, 0.57], linestyle='--', color='black')
        plt.annotate('Saponite', xy=(0.21, 0.57), xytext = (0.24, 0.57), size=10) 
        #plt.plot([0.33, -0.33], [0.33, 0.33], linestyle=':', color='black', linewidth = 0.2) #olivine
        #plt.plot([0.25, -0.25], [0.5, 0.5], linestyle=':', color='black', linewidth = 0.2) #pyroxene
        plt.annotate('Si+Al', xy=(0., 1.), xytext = (-0.1, 1.03), size=13) 
        plt.annotate('Mg', xy=(-0.55,0.), xytext = (-0.64,-0.03), size=13) 
        plt.annotate('Fe', xy=(0.55,0.), xytext = (0.52, -0.03), size=13)
    
    if type=='sulfide':
        plt.plot(-0.25, 0.5, 'ro', markersize=5)
        plt.plot(0.0, 0.47, 'bo', markersize = 5)
        plt.plot([-0.25, -0.22], [0.5, 0.56], color='red', linewidth = 7)
        #plt.plot(-0.17,0.67, 'go',  markersize = 5)
        plt.annotate('S', xy=(0., 1.), xytext = (-0.02, 1.03), size=15) 
        plt.annotate('Fe', xy=(-0.5,0.), xytext = (-0.6,-0.03), size=15) 
        plt.annotate('Ni', xy=(0.55,0.), xytext = (0.52, -0.03), size=15)
        plt.annotate('Troilite', xy=(-0.3, 0.5), xytext = (-0.46, 0.48), size=11) 
        plt.annotate('Pentlandite', xy=(-0.1, 0.47), xytext = (-0.06,0.5), size=11) 
        plt.annotate('Pyrrhotite', xy=(-0.25, 0.55), xytext = (-0.5,0.57), size=11) 
        #plt.annotate('Pyrite', xy=(-0.22,0.67), xytext = (-0.36,0.67), size=11) 
        
    #### Deal with the case when there is no masks: ----------------------------------------
    if mask==None: 
        mask = []
        for i in range(len(dataset)):
            mask.append(dataset[i][0])
            mask[i].data[:,:] = True
            
    #### Plot the data in the ternary: ----------------------------------------
    for i in range(len(dataset)):
        plt.plot(tri_x[i].data[mask[i]].flatten(), tri_y[i].data[mask[i]].flatten(), '.', markersize = markersize, label = legend[i], color = color[i])

    plt.legend()
    plt.axis('off')
    
#-----------------------------------------------------------------------------------------
    

def Ternary_Contour(dataset, type=None, cmap=None, legend=None, mask=None, levels=None, alpha=None, ContourValues=None, ContLines=None, ContColourFill=None, DataPointDisp=None):
    """ Plots a contour map of the given data : Gives an idea of the density of the points in the ternary"""
        
    def Tern_Base():  
        """ Uses Corentin Le Guillou scripts to plot a ternay diagram"""
        #### Plot initialisation: ---------------------------------------------------------
        ax.plot([-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0],  [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], color = 'black', marker="_")
        ax.plot([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],  [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0], color='black', marker="_")
        ax.plot([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], [0,0,0,0,0,0,0,0,0,0,0], color='black', marker="|")
        if type == 'silicate':
            ax.plot([0.3, -0.3], [0.4, 0.4], linestyle='--', color='black')
            ax.annotate('Serpentine', xy=(0.3, 0.4), xytext = (0.33, 0.4), size=10)
            ax.plot([0.21, -0.21], [0.57, 0.57], linestyle='--', color='black')
            ax.annotate('Saponite', xy=(0.21, 0.57), xytext = (0.24, 0.57), size=10)
            ax.annotate('Si+Al', xy=(0., 1.), xytext = (-0.10, 1.03), size=14)
            ax.annotate('Mg', xy=(-0.55,0.), xytext = (-0.64,-0.03), size=14)
            ax.annotate('Fe', xy=(0.55,0.), xytext = (0.53, -0.03), size=14)
        if type == 'sulfide':
            ax.plot(-0.25, 0.5, 'ro', markersize=5)
            ax.plot(0.0, 0.47, 'bo', markersize = 5)
            ax.plot([-0.25, -0.22], [0.5, 0.56], color='purple', linewidth = 7)
            ax.plot(-0.17,0.67, 'go',  markersize = 5)
            ax.annotate('S', xy=(0., 1.), xytext = (-0.05, 1.03), size=14)
            ax.annotate('Fe', xy=(-0.55,0.), xytext = (-0.56,-0.03), size=14)
            ax.annotate('Ni', xy=(0.55,0.), xytext = (0.52, -0.03), size=14)
            ax.annotate('Troilite', xy=(-0.25, 0.5), xytext = (-0.4, 0.5), size=11)
            ax.annotate('Pentlandite', xy=(0.0, 0.47), xytext = (0.03,0.47), size=11)
            ax.annotate('Pyrrhotite', xy=(-0.235, 0.55), xytext = (-0.42,0.55), size=11)
            ax.annotate('Pyrite', xy=(-0.17,0.67), xytext = (-0.32,0.67), size=11)
    
    #### Initialisation: ----------------------------------------------------------------- 
    tri_x=[]
    tri_y=[]
    
    x = [[] for i in range(len(dataset))]
    y = [[] for i in range(len(dataset))]
    
    xmin, xmax = -0.6, 0.6
    ymin, ymax = -0.1, 1.1
    
    fig, ax = plt.subplots()
    Tern_Base()
    
    ContourLabelTextSize = 0.5
    ContourLineThickness = 0.4
    ContourLineStyle = '-'
    DataPointSize = 0.05
  
    #### Run the for loop to plot the different data in a ternary: ----------------------
    for i in range (len(dataset)):
        tri_y.append(dataset[i][0].data[mask[i]].flatten()/(dataset[i][0].data[mask[i]].flatten()+dataset[i][1].data[mask[i]].flatten()+dataset[i][2].data[mask[i]].flatten()))
        tri_x.append(np.array(-0.5)*(np.array(1.)-tri_y[i])+(dataset[i][2].data[mask[i]].flatten()/(dataset[i][0].data[mask[i]].flatten()+dataset[i][1].data[mask[i]].flatten()+dataset[i][2].data[mask[i]].flatten())))

        for j in range(len(tri_x[i])):
            if math.isnan(tri_x[i][j])==False:
                x[i].append(tri_x[i][j])
                y[i].append(tri_y[i][j])

        # Peform the kernel density estimate for each data
        X, Y = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([x[i],y[i]])
        kernel = st.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

        if ContColourFill == 'y':
            ContColour_Ternary = plt.contourf(X,Y, Z, cmap=cmap[i], alpha=alpha, locator = ticker.MaxNLocator(prune = 'lower', nbins=levels), zorder=5)
            plt.colorbar(label = legend[i])
        if ContLines == 'n' and ContourValues == 'y':
            ContourLineThickness = 0
            cset_Ternary = plt.contour(X, Y, Z, colors='k', alpha=1, linewidths = ContourLineThickness, linestyles = ContourLineStyle, locator = ticker.MaxNLocator(prune = 'lower', nbins=levels), zorder=10) # Drawing contour lines.
            if ContourValues == 'y':
                ax.clabel(cset_Ternary, inline=1, fontsize=ContourLabelTextSize, zorder=15) # Labelling contour levels within the contour lines.
        if ContLines == 'y':
            cset_Ternary = plt.contour(X, Y, Z, colors='k', alpha=1, linewidths = ContourLineThickness, linestyles = ContourLineStyle, locator = ticker.MaxNLocator(prune = 'lower', nbins=levels), zorder=10) # Drawing contour lines.
            if ContourValues == 'y':
                ax.clabel(cset_Ternary, inline=1, fontsize=ContourLabelTextSize, zorder=15) # Labelling contour levels within the contour lines.
        if DataPointDisp=='y':
            ax.scatter(x[i], y[i], color='grey', alpha=0.45, s=DataPointSize, zorder=5)
            
    ax.axis('off')