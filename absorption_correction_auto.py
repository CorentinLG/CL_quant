# -*- coding: utf-8 -*-

from math import sin
from scipy import pi
import numpy as np
import copy
import hyperspy.api as hs
from water_content import water_content
from Correction import correction


def compare (Quant1, Quant2, line1, line2, Elt_rat, mt, d, Full_mask, D, Crit_mt, Min_thickness):
    # recalculate the difference between the two lines of interest for next iteration             
    SubQuant = []
    for i in range (len(Quant1)):
        if line1 in Quant1[i].metadata.Sample.xray_lines : SubQuant.append(Quant1[i])
    for i in range (len(Quant2)):
        if line2 in Quant2[i].metadata.Sample.xray_lines : SubQuant.append(Quant2[i])
        
    Dev = (SubQuant[0].data - Elt_rat * SubQuant[1].data)/SubQuant[0].data
    rat = SubQuant[0].data/(Elt_rat*SubQuant[1].data)
    rat = np.nan_to_num(rat)
    print('line1/line2: mean deviation (%) =', np.nanmean(Dev)*100,'max=',np.nanmax(Dev)*100, 'min=', np.nanmin(Dev)*100, '...Modifying thickness for next iteration')
    
    converged_mask = hs.signals.BaseSignal(np.full((D[:-1]), False))
    thin_mask =  hs.signals.BaseSignal(np.full((D[:-1]), False))

    Thin_pix = []
    if len(D)==1:
        if (mt/(d*10**-7)).data<Min_thickness:Thin_pix = [1]
    if len(D)==2:
        for i in range(D[0]):
            if Full_mask.data[i]== False:
                if (mt.inav[i]/(d.inav[i]*10**-7)).data<Min_thickness:
                    Thin_pix.append((i))
                    thin_mask.data[i] = True
                    Dev[i] = np.nan
                if abs(Dev[i])<(Crit_mt/100) : 
                    converged_mask.data[i] = True
    if len(D)==3:
        for i in range(D[0]):
            for j in range (D[1]):
                if Full_mask.data[i][j] == False:
                    if (mt.inav[j,i]/(d.inav[j,i]*10**-7)).data<Min_thickness:
                        Thin_pix.append((i, j))
                        thin_mask.data[i][j] = True
                        Dev[i][j]= np.nan
                    if abs(Dev[i][j])<(Crit_mt/100) : 
                        converged_mask.data[i][j] = True
    if Thin_pix!= []:
        print('pixel(s) '+str(Thin_pix)+' seem(s) to have a thickness lower than ' + str(Min_thickness)+' nm. Correction is at its minimum, Issue with SNR? k-factors?, fitting model? element ratio?')
    return Dev, rat, thin_mask, converged_mask    

def absorption_correction_auto(result, s, kfactors, 
                               line1 = 'Fe_Ka', line2 = 'Fe_La', Elt_rat = 1, 
                               d = 3, t = 150, tilt_stage = 0, 
                               Crit_AC = 0.1, Crit_mt = 0.5, 
                               navigation_mask = None, composition_units = 'atomic',
                               water = False, valence = None, Min_thickness = 10):
    
    """
    This function quantifies TEM-EDS data based on the Cliff-Lorimer method and takes absorption correction into account (a function of sample thickness, 
    density, and take-off angle). Most importantly, it offers the possibility to determine the mass thickness automatically, if the ratio of two elements is known (for instance: the K and L lines of iron has a ratio of 1; the Si/O ratio in olivine (Mg2SiO4) is (4*16/1*28)= 2.29)

    This is a double iterative process : quantification is ran, then absorption is calculated iteratively, and then mass absorption is modified untill similar quantification is obtained for the two selected lines.
    Convergence criteria can be given for absorption correction itaration (Ac_Crit) and mass thickness convergence (Crit_mt). This sequence is performed until two successive quantification are similar within the converegnce criteria. 
	For calculation principles, See William, D. B., & Carter, C. B. (1996). Transmission electron microscopy: a textbook for materials science. Edit. Plenum Press, New York & London.]    
    
    Parameters
    ----------
    result : a list of signals, the integrated gaussian intensitiy of each X-ray lines.
	
    s: the signal being quantified. Its metadata (elements and xray_lines) must be documented.
    
    kfactors : a list of k-factors (must have the same length as xray_lines from s.metadata and as intensities)
    
    line1 : the less absorbed x-ray line (usually the one with the highest energy)

    line2 : the most absorbed x-ray line (usually the one with the lowest energy)
    
    Elt_rat : the known mass ratio of two elements (see above)
    
    d : density in g.cm**3 (can be a single value or have the same navigation dimension as the dataset. In this case a density map is required as well.
    
    t : thickness is in nm (can be a single value or have the same navigation dimension as the dataset, but must have the same type and dimension as d).
    
    tilt_stage : the sample holder tilt stage, used to calculate the absorption pathways (in degrees).
    
    Crit_AC : the convergence criterium for the absorption correction.
    
    Crit_mt : the convergence criterium for the mass thickness determination
    
    navigation_mask : None By default, a BaseSignal of booleans otherwise.
    
	composition_units : the output units of the quantification
    
    water : if == True ; water content is recalculated based on electroneutrality assumption. Requires to input the 'valence' argument
    
    valence : a list of valence for each element
    
    Min_thickness : the minimum thickness value below which mass thickness determination is stopped.

    Returns
    ------
    A list of quantified elemental composition of the sample (in atomic percent) taking into account the absorption correction.
    The water content (in wt.%); only if water = True
    The mass thickness

    Examples
    --------
    ### if the Fe_K and Fe_L lines are used to constrain the mass thickness : 
    result = m.get_lines_intensity()
    factors = kfactors(result_cor)
    val = [0, 2, 0, 2.6, 2.6, 2, -2, 4] #Avec Fe_L
    Quant, H2O, mt = absorption_correction_auto (result, s, factors, line1 = 'Fe_Ka', line2 = 'Fe_La', Elt_rat = 1, d = d, t = t_map, tilt_stage = 31, Crit_AC = 0.1, Crit_mt = 0.5, water = True, valence = val)
    
    """
    
    
    if len(s.data.shape) ==4: raise ValueError('The absorption correction routine is not implemented for 2D signals')  # dimensions are gonna that of a 3D signals, 3 navigation dimension.
    #if t== 0 or d==0 :raise ("thickness and density can't be equal to 0")
    
    D = []
    for i in range (len(s.data.shape)-1):
        D.append(s.data.shape[i])
    D.append(len(s.metadata.Sample.xray_lines))
    
    mt = hs.signals.Signal1D(np.zeros(D[:-1]+[1]))
    if navigation_mask is None:
        navigation_mask = hs.signals.BaseSignal(np.full((D[:-1]), False))
    if type(d)==float or type(d)==int and type(t)==float or type(t)==int:
        mt.isig[0] = d*(t*10**-7)
        t = mt/(d*10**-7)
        d = mt/(t*10**-7)
    elif type(t)!=(int and float) and type(d)!= (int and float):
        mt.isig[0] = d.isig[0]*t.isig[0]*10**-7
    mt.data[navigation_mask.data] = np.nan
    
    mt_f = copy.deepcopy(mt)
    Dev = np.full((D[:-1]), 1.)
    Dev_f = np.full((D[:-1]), 1.)
    rat = np.ones((D[:-1]))
    dif = np.ones((D[-1:]+D[:-1]))

    Full_mask = copy.deepcopy(navigation_mask)
    
    s.metadata.Acquisition_instrument.TEM.tilt_stage = tilt_stage
    alpha = (s.metadata.Acquisition_instrument.TEM.Detector.EDS.elevation_angle + 
         s.metadata.Acquisition_instrument.TEM.tilt_stage)*pi/180
                
    elts = []
    for i in range(0, len(s.metadata.Sample.xray_lines)):
        elts.append(result[i].metadata.Sample.elements)   
    print('elts', elts)
    
    if len(s.data.shape)!=1:
        if len(s.data.shape)==2: 
            a= D[0]     
        if len(s.data.shape)==3: 
            a= D[0]*D[1]
        l = a-np.count_nonzero(navigation_mask.data)
        print('number of pixels to deal with: ', l, 'total number of pixels: ', a)        

    ### If the two elements used to calculate the mt value are the same, i.e. if ones uses Fe_Ka and Fe_La for instance, then two 
    ## quantifications are made and the intensity of the other line is brought to zero.
    if line1[:2] == line2[:2]:
        result_int_K = copy.deepcopy(result) # since result is manipulated many times, better copy it deeply before.
        result_int_L = copy.deepcopy(result)
        result_mod_K = copy.deepcopy(result) 
        result_mod_L = copy.deepcopy(result)
        Ac_K = hs.signals.Signal1D(np.ones((D)))
        wt_K = hs.signals.Signal1D(np.ones((D)))
        Ac_L = hs.signals.Signal1D(np.ones((D)))
        wt_L = hs.signals.Signal1D(np.ones((D)))
        for i in range (0,len(result_int_K)):
            if line1 in result_int_L[i].metadata.Sample.xray_lines : result_int_L[i].data = np.zeros((D[:-1]))
            elif line2 in result_int_K[i].metadata.Sample.xray_lines : result_int_K[i].data = np.zeros((D[:-1]))
        Quant_K = s.quantification(method="CL", intensities=result_int_K, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)     
        Quant_L = s.quantification(method="CL", intensities=result_int_L, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)
        Quant_K2 = Quant_K
        Quant_L2 = Quant_L
        Quant_f = copy.deepcopy(Quant_K)

        while (abs(Dev)>Crit_mt/100).any():                    
            mt.isig[0] = rat*mt.isig[0]
            dif = np.ones((D[-1:]+D[:-1]))
			
            while (abs(dif) > Crit_AC/100).any():
                Quant_K = Quant_K2
                #Calculation of intensities corrected for absorption
                result_mod_K, Ac_K = correction (s, elts, Quant_K, result_int_K, result_mod_K, alpha, mt, Full_mask, D, Dev, dif, Crit_mt, Ac_K, wt_K)
                #New quantification using corrected intensities
                Quant_K2 = s.quantification(method="CL", intensities=result_mod_K, factors=kfactors, composition_units='weight', navigation_mask = Full_mask, plot_result=False)
                #Compares the relative difference between previous and last quantification (convergence test value)
                for i in range (len(Quant_K)):
                    dif[i]= abs(Quant_K[i].data-Quant_K2[i].data)/Quant_K[i].data
                print('AbsCor_high energy line')        
            dif = np.ones((D[-1:]+D[:-1]))
			
            while (abs(dif) > Crit_AC/100).any():
                Quant_L = Quant_L2
                result_mod_L, Ac_L  = correction (s, elts, Quant_L, result_int_L, result_mod_L, alpha, mt, Full_mask, D, Dev, dif, Crit_mt, Ac_L, wt_L)
                Quant_L2 = s.quantification(method="CL", intensities=result_mod_L, factors=kfactors, composition_units='weight', navigation_mask = Full_mask, plot_result=False)
                for i in range (len(Quant_L)):
                    dif[i]= abs(Quant_L[i].data-Quant_L2[i].data)/Quant_L[i].data     
                print('AbsCor_low energy line')
            # recalculate the difference between the two lines of interest and correct the thickness accordingly for next iteration             
            Dev, rat, thin_mask, converged_mask = compare(Quant_K, Quant_L, line1, line2, Elt_rat, mt, d, Full_mask, D, Crit_mt, Min_thickness)
            
            if len(s.data.shape)==1:
                for k in range (D[-1]): 
                    Quant_f[k] = Quant_K[k]
                mt_f.data = mt.data
            else:
                Dev_f[thin_mask+converged_mask] = Dev[thin_mask+converged_mask]
                mt_f.data[thin_mask+converged_mask] = mt.data[thin_mask+converged_mask]  
                for k in range (D[-1]): Quant_f[k].data[thin_mask+converged_mask] = Quant_K[k].data[thin_mask+converged_mask]
                Full_mask = Full_mask+thin_mask+converged_mask
                print('% of pixels left to analyze', (a-np.count_nonzero(Full_mask.data))*100/l)
            
        ###########################################################################
        ######  Further quantification using WATER corrected Abs Correction factors
        if water == True:
            if valence == None : raise ValueError('a list of valence for each xray_lines or element should be provided as argument')
            Quant_f = hs.material.weight_to_atomic(Quant_f, elements='auto')
            Quant_f, H2O_f = water_content(s, Quant_f, valence)
            Quant_f = hs.material.atomic_to_weight(Quant_f, elements='auto')
            print('Water computed')

        if composition_units=='atomic': Quant_f = hs.material.weight_to_atomic(Quant_f, elements='auto')
        
###########################################################################################################        
### If two different elements are used, then only one quantification is performed and the convergence criteria are internal to this calculation
    if line1[:2] != line2[:2]:      
        result_int = copy.deepcopy(result)
        result_mod = copy.deepcopy(result)
        Ac = hs.signals.Signal1D(np.ones((D)))
        wt = hs.signals.Signal1D(np.ones((D)))
        
        Quant = s.quantification(method="CL", intensities=result_int, factors=kfactors, composition_units='weight', navigation_mask=navigation_mask, plot_result=False)   
        Quant2 = Quant
        Quant_f = copy.deepcopy(Quant)
        
        while (abs(Dev)>Crit_mt/100).any():                    
            mt.isig[0] = rat*mt.isig[0]
            dif = np.ones((D[-1:]+D[:-1]))
            while (abs(dif) > Crit_AC/100).any():
                Quant = Quant2
                #Calculation of intensities corrected for absorption
                result_mod, Ac = correction(s, elts, Quant, result_int, result_mod, alpha, mt, Full_mask, D, Dev, dif, Crit_mt, Ac, wt)
                #New quantification using corrected intensities
                Quant2 = s.quantification(method="CL", intensities=result_mod, factors=kfactors, composition_units='weight', navigation_mask=Full_mask, plot_result=False)          
                #Compares the relative difference between two successive quantification (convergence test)
                for i in range (len(Quant)):
                    dif[i]= abs(Quant[i].data-Quant2[i].data)/Quant[i].data
                print('AbsCor')
            # recalculate the difference between the two lines of interest anc orrect the thickness accordingly for next iteration             
            Dev, rat, thin_mask, converged_mask = compare(Quant2, Quant2, line1, line2, Elt_rat, mt, d, Full_mask, D, Crit_mt, Min_thickness)
            
            if len(s.data.shape)==1:
                for k in range (D[-1]): Quant_f[k] = Quant_K[k]
                mt_f.data = mt.data
            else:
                Dev_f[thin_mask+converged_mask] = Dev[thin_mask+converged_mask]
                mt_f.data[thin_mask+converged_mask] = mt.data[thin_mask+converged_mask]  
                for k in range (D[-1]): Quant_f[k].data[thin_mask+converged_mask] = Quant2[k].data[thin_mask+converged_mask]           
                Full_mask = Full_mask+thin_mask+converged_mask
                print('% of pixels left to analyze', (a-np.count_nonzero(Full_mask.data))*100/l)
            
        #######################################
        #New quantification using WATER corrected Abs Correction factors
        ############################################
        if water == True:
            if valence == None : raise ValueError('a list of valence for each xray_lines or element should be provided as argument')
            Quant_f = hs.material.weight_to_atomic(Quant_f, elements='auto')
            Quant_f, H2O_f = water_content (s, Quant_f, valence)
            print('Water computed')
            Quant_f = hs.material.atomic_to_weight(Quant_f, elements='auto')

        if composition_units=='atomic': Quant_f = hs.material.weight_to_atomic(Quant_f, elements='auto')

    if water==False : H2O_f = np.zeros((D[:-1]))
    if len(s.data.shape)==2:                             
        for i in range(D[0]):
            if navigation_mask.data[i]== True:
                mt_f.inav[i].data = np.nan
                H2O_f[i] = np.nan
                for k in range (len(Quant_f)):
                    Quant_f[k].data[i]=np.nan
    if len(s.data.shape)==3:                             
        for i in range(D[0]):
            for j in range (D[1]):
                if navigation_mask.data[i][j]== True:
                    mt_f.inav[j,i].data = np.nan
                    H2O_f[i][j] = np.nan
                    for k in range (len(Quant_f)):
                        Quant_f[k].data[i][j]=np.nan
    if water == True: return Quant_f, H2O_f, mt_f
    if water == False: return Quant_f, mt_f