def plot_color_map(N, G, legend, colorscale, save = False, filename='filename'):
    """
    ** N is the number of phases that you want to plot
    ** legend is their names (strings)
    ** coloscale must be given in the format : colorscale = [[1,0,0], [0,1,0], [0,0,1]]
    """
    
    
    from skimage import color
    import matplotlib.pyplot as plt
    import numpy as np
    RGB = []
    cmap=plt.cm.hsv
    cmaplist=np.zeros([N,4])
    RGB_final = 0

    for i in range (0,N):
        cmaplist[i]=(colorscale[i,:][0],colorscale[i,:][1],colorscale[i,:][2],1)

    #create the next map
    cmap = cmap.from_list('custum cmap', cmaplist, N)


    for i in range (N):
         RGB.append(0)
         RGB[i] = color.gray2rgb(G.data[i,:,:])
         for j in range (3):
             RGB[i][:,:,j] = RGB[i][:,:,j]*colorscale[i][j]
         RGB_final += RGB[i]

    plt.figure()
    plt.imshow(RGB_final,cmap=cmap)

    cb=plt.colorbar()

    step = 1/N
    ticks = [0.5/(N)]
    for i in range(N-1):
        ticks.append(ticks[i]+step)
    cb.set_ticks(ticks)

    #cb.ax.autoscale(True)
    #cb.set_ticks([0.125, 0.375, 0.75, 1.])
    cb.set_ticklabels(legend)
    
    if save==True:
        plt.gcf()
        plt.savefig(filename+'.png', bbox_inches='tight')
        plt.savefig(filename+'.pdf', bbox_inches='tight')


def triplot(dataset, mask=None, legend= None, type='Silicate'):
    """ Uses matplotlib to plot a ternay diagram
    ** dataset should be an array containing the values to plot, in the form [A, B, C] where A is the top (Si+Al), B is the left corner ('Fe') and C the right corner ('Mg') of the triangle.
    Several dataset can be plotted, then dataset in in the Form [[A1, B1, C1], [A2, B2, C2], ...]
    ** mask : can be applied to the dataset 
    ** legend : A list for legending different datasets
    ** type : Either 'silicate' [Si+Al, Mg, Fe] or 'sulfide' [S, Ni, Fe]. 
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    tri_y=[]
    tri_x=[]
    for i in range (len(dataset)):
        tri_y.append(dataset[i][0]/(dataset[i][0]+dataset[i][1]+dataset[i][2]))
        tri_x.append(np.array(-0.5)*(np.array(1.)-tri_y[i])+(dataset[i][2]/(dataset[i][0]+dataset[i][1]+dataset[i][2])))
        
    plt.figure()
    plt.plot([-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0],  [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], color = 'black', marker="_")
    plt.plot([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],  [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0], color='black', marker="_")
    plt.plot([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], [0,0,0,0,0,0,0,0,0,0,0], color='black', marker="|")
    
    if type=='silicate':
        plt.plot([0.3, -0.3], [0.4, 0.4], linestyle='--', color='black')
        plt.annotate('Serpentine', xy=(0.3, 0.4), xytext = (0.33, 0.4), size=10) 
        plt.plot([0.21, -0.21], [0.57, 0.57], linestyle='--', color='black')
        plt.annotate('Saponite', xy=(0.21, 0.57), xytext = (0.24, 0.57), size=10) 
        plt.plot([0.33, -0.33], [0.33, 0.33], linestyle=':', color='black', linewidth = 0.2) #olivine
        plt.plot([0.25, -0.25], [0.5, 0.5], linestyle=':', color='black', linewidth = 0.2) #pyroxene
        plt.annotate('Si+Al', xy=(0., 1.), xytext = (-0.07, 1.03), size=17)
        plt.annotate('Fe', xy=(-0.55,0.), xytext = (-0.58,-0.03), size=17) 
        plt.annotate('Mg', xy=(0.55,0.), xytext = (0.52, -0.03), size=17)
    
    if type=='Classif_phyllo':            
        plt.annotate('M+', xy=(0., 1.), xytext = (-0.05, 1.03), size=15)
        plt.annotate('4Si', xy=(-0.55,0.), xytext = (-0.6,-0.03), size=15) 
        plt.annotate('3R2+', xy=(0.55,0.), xytext = (0.52, -0.03), size=15)
    
    if type=='sulfide':
        plt.plot(-0.25, 0.5, 'ro', markersize=5)
        plt.plot(0.0, 0.47, 'bo', markersize = 5)
        plt.plot([-0.25, -0.22], [0.5, 0.56], color='purple', linewidth = 7)
        plt.plot(-0.17,0.67, 'go',  markersize = 5)
        plt.annotate('S', xy=(0., 1.), xytext = (-0.05, 1.03), size=13) 
        plt.annotate('Fe', xy=(-0.55,0.), xytext = (-0.56,-0.03), size=13) 
        plt.annotate('Ni', xy=(0.55,0.), xytext = (0.52, -0.03), size=13)
        plt.annotate('Troilite', xy=(-0.25, 0.5), xytext = (-0.4, 0.5), size=11) 
        plt.annotate('Pentlandite', xy=(0.0, 0.47), xytext = (0.03,0.47), size=11) 
        plt.annotate('Pyrrhotite', xy=(-0.235, 0.55), xytext = (-0.42,0.55), size=11) 
        plt.annotate('Pyrite', xy=(-0.17,0.67), xytext = (-0.32,0.67), size=11) 
        
    if type=='Organics':            
    	plt.annotate('Aromatics+Olefinics', xy=(0., 1.), xytext = (-0.25, 1.03), size=15)
    	plt.annotate('Ketones+Phenols+Nitriles', xy=(-0.40,0.), xytext = (-0.7,-0.075), size=15) 
    	plt.annotate('Aliphatics', xy=(0.55,0.), xytext = (0.4, -0.075), size=15)
        
    if mask==None: 
        mask = []
        for i in range(len(dataset)):
            mask.append(dataset[i][0])
            mask[i].data[:,:] = True
            
    for i in range(len(dataset)):
        plt.plot(tri_x[i].data[mask[i]].flatten(), tri_y[i].data[mask[i]].flatten(), '.', markersize = 0.6, label = legend[i])

    lgnd=plt.legend(fontsize = 13)
    for i in range (len(dataset)):lgnd.legendHandles[i]._legmarker.set_markersize(6)

    plt.legend()
    plt.axis('off')


def triplot(dataset, mask=None, legend= None, type='Silicate'):
    """ Uses matplotlib to plot a ternay diagram
    ** dataset should be an array containing the values to plot, in the form [A, B, C] where A is the top (Si+Al), B is the left corner ('Fe') and C the right corner ('Mg') of the triangle.
    Several dataset can be plotted, then dataset in in the Form [[A1, B1, C1], [A2, B2, C2], ...]
    ** mask : can be applied to the dataset 
    ** legend : A list for legending different datasets
    ** type : Either 'silicate' [Si+Al, Mg, Fe] or 'sulfide' [S, Ni, Fe]. 
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    tri_y=[]
    tri_x=[]
    for i in range (len(dataset)):
        tri_y.append(dataset[i][0]/(dataset[i][0]+dataset[i][1]+dataset[i][2]))
        tri_x.append(np.array(-0.5)*(np.array(1.)-tri_y[i])+(dataset[i][2]/(dataset[i][0]+dataset[i][1]+dataset[i][2])))
        
    plt.figure()
    plt.plot([-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0],  [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], color = 'black', marker="_")
    plt.plot([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],  [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0], color='black', marker="_")
    plt.plot([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5], [0,0,0,0,0,0,0,0,0,0,0], color='black', marker="|")
    
    if type=='silicate':
        plt.plot([0.3, -0.3], [0.4, 0.4], linestyle='--', color='black')
        plt.annotate('Serpentine', xy=(0.3, 0.4), xytext = (0.33, 0.4), size=10) 
        plt.plot([0.21, -0.21], [0.57, 0.57], linestyle='--', color='black')
        plt.annotate('Saponite', xy=(0.21, 0.57), xytext = (0.24, 0.57), size=10) 
        plt.plot([0.33, -0.33], [0.33, 0.33], linestyle=':', color='black', linewidth = 0.2) #olivine
        plt.plot([0.25, -0.25], [0.5, 0.5], linestyle=':', color='black', linewidth = 0.2) #pyroxene
        plt.annotate('Si+Al', xy=(0., 1.), xytext = (-0.07, 1.03), size=17)
        plt.annotate('Fe', xy=(-0.55,0.), xytext = (-0.58,-0.03), size=17) 
        plt.annotate('Mg', xy=(0.55,0.), xytext = (0.52, -0.03), size=17)
    
    if type=='Classif_phyllo':            
        plt.annotate('M+', xy=(0., 1.), xytext = (-0.05, 1.03), size=15)
        plt.annotate('4Si', xy=(-0.55,0.), xytext = (-0.6,-0.03), size=15) 
        plt.annotate('3R2+', xy=(0.55,0.), xytext = (0.52, -0.03), size=15)
    
    if type=='sulfide':
        plt.plot(-0.25, 0.5, 'ro', markersize=5)
        plt.plot(0.0, 0.47, 'bo', markersize = 5)
        plt.plot([-0.25, -0.22], [0.5, 0.56], color='purple', linewidth = 7)
        plt.plot(-0.17,0.67, 'go',  markersize = 5)
        plt.annotate('S', xy=(0., 1.), xytext = (-0.05, 1.03), size=13) 
        plt.annotate('Fe', xy=(-0.55,0.), xytext = (-0.56,-0.03), size=13) 
        plt.annotate('Ni', xy=(0.55,0.), xytext = (0.52, -0.03), size=13)
        plt.annotate('Troilite', xy=(-0.25, 0.5), xytext = (-0.4, 0.5), size=11) 
        plt.annotate('Pentlandite', xy=(0.0, 0.47), xytext = (0.03,0.47), size=11) 
        plt.annotate('Pyrrhotite', xy=(-0.235, 0.55), xytext = (-0.42,0.55), size=11) 
        plt.annotate('Pyrite', xy=(-0.17,0.67), xytext = (-0.32,0.67), size=11) 
        
    if type=='Organics':            
            plt.annotate('Aromatics+Olefinics', xy=(0., 1.), xytext = (-0.25, 1.03), size=15)
            plt.annotate('Ketones+Phenols+Nitriles', xy=(-0.40,0.), xytext = (-0.7,-0.075), size=15) 
            plt.annotate('Aliphatics', xy=(0.55,0.), xytext = (0.4, -0.075), size=15)
            
    if mask==None: 
        mask = []
        for i in range(len(dataset)):
            mask.append(dataset[i][0])
            mask[i].data[:,:] = True
            
    for i in range(len(dataset)):
        plt.plot(tri_x[i].data[mask[i]].flatten(), tri_y[i].data[mask[i]].flatten(), '.', markersize = 0.6, label = legend[i])

    lgnd=plt.legend(fontsize = 13)
    for i in range (len(dataset)):lgnd.legendHandles[i]._legmarker.set_markersize(6)

    plt.legend()
    plt.axis('off')


def Ternary_Contour(dataset, type=None, cmap=None, legend=None, mask=None, levels=None, alpha=None, ContourValues=None, ContLines=None, ContColourFill=None, DataPointDisp=None):
    """ Plots a contour map of the given data : Gives an idea of the density of the points in the ternary"""
    
    import numpy as np
    import math
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    import matplotlib.ticker as ticker
    import matplotlib.tri as tri
    import scipy.stats as st
    from scipy.interpolate import griddata        
    
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
        
        if type=='Classif_phyllo':            
            plt.annotate('M+', xy=(0., 1.), xytext = (-0.05, 1.03), size=15)
            plt.annotate('4Si', xy=(-0.55,0.), xytext = (-0.6,-0.03), size=15) 
            plt.annotate('3R2+', xy=(0.55,0.), xytext = (0.52, -0.03), size=15)

            ax.plot(np.array(-0.5+(1.5/(0.5+1.5))), np.array(0.), 'ro', markersize=8)
            ax.plot(np.array(-0.5*(1-(0.5/(3.5/4+1+0.5)))+(1/(3.7/4+1+0.3))), np.array(0.5/(3.5/4+1+0.5)), 'ro', markersize=8)
            ax.plot(np.array(-0.5*(1-(0.3/(1+0.3)))), np.array(0.3/(1+0.3)), 'ro', markersize=8)
            
            ax.plot(np.array(0.), np.array(0.), 'ro', markersize=8)
            
            ax.annotate('serpentine', xy=(0.14, -0.08), xytext = (0.14, -0.08), size=11)
            ax.annotate('saponite', xy=(-0.03, 0.24), xytext = (-0.03, 0.24), size=11)
            ax.annotate('nontronite (FeIII)', xy=(-0.56, 0.28), xytext = (-0.56, 0.28), size=11)
            ax.annotate('Talc', xy=(0.0, -0.08), xytext = (0.0, -0.08), size=11)
            
        if type=='Organics':            
            plt.annotate('Aromatics+Olefinics', xy=(0., 1.), xytext = (-0.25, 1.03), size=15)
            plt.annotate('Ketones+Phenols+Nitriles', xy=(-0.40,0.), xytext = (-0.7,-0.075), size=15) 
            plt.annotate('Aliphatics', xy=(0.55,0.), xytext = (0.4, -0.075), size=15)

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
    if mask==None: 
        mask = []
        for i in range(len(dataset)):
            mask.append(dataset[i][0])
            mask[i].data[:,:] = True
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


def plot_histo(dataset, mask=None, legend=None, xlabel= None, ylabel = None, linestyle=None):
    import matplotlib.pyplot as plt
    import numpy as np
    histo=[]
    D=[]
    plt.figure()
    for i in range (len(dataset)):
        D.append(dataset[i].deepcopy())
        if mask != None: 
            if mask[i] != None: D[i].data[mask[i].data] = np.nan
        histo.append(D[i].get_histogram())
        plt.plot(histo[i].axes_manager[0].axis, (histo[i].data*100/np.sum(histo[i].data)), 
                label=legend[i],
                linestyle = linestyle[i]
               )
    
    #if linestyle != None: plt.set_linestyle(linestyle)
    #leg = ax.legend();
    plt.xlim([None, None]) 
    plt.grid()
    plt.xlabel(xlabel, size=14)
    plt.ylabel(ylabel, size=14)
    plt.legend()
