"""Moving to 2D"""

# -*- coding: utf-8 -*-

from math import exp, sin
from scipy import pi
import numpy as np
import copy
import hyperspy.api as hs
from water_content import water_content


def correction(s, elts, Quant, result_int, result_mod, alpha, mt, navigation_mask, D, Dev, Crit, Ac, wt):
    
    F_Si = 0
    F_S = 0.01
    F_C = 0  #0.001
    F_Ca = 0.01 # 0.01
    F_O = 0.003
    F_Fe = 0.03
    
    wt_Si, wt_S, wt_C, wt_Ca, wt_O = wt.isig[0], wt.isig[0], wt.isig[0], wt.isig[0], wt.isig[0]
    for k in range(D[-1]):                  
        if 'Si_Ka' in Quant[k].metadata.Sample.xray_lines: wt_Si.data = Quant[k].data
        if 'S_Ka' in Quant[k].metadata.Sample.xray_lines:wt_S.data = Quant[k].data
        if 'C_Ka' in Quant[k].metadata.Sample.xray_lines:wt_C.data = Quant[k].data
        if 'Ca_Ka' in Quant[k].metadata.Sample.xray_lines:wt_Ca.data = Quant[k].data
        if 'O_Ka' in Quant[k].metadata.Sample.xray_lines:wt_O.data = Quant[k].data 

    if len(s.data.shape)==1: 
        for k in range(D[-1]):
            wt.isig[k] = Quant[k]
            Ac.data[k] = hs.material.mass_absorption_mixture(wt.data, elts, energies = s.metadata.Sample.xray_lines[k])    
            if 'C_Ka' in Quant[k].metadata.Sample.xray_lines: 
                Ac.data[k] = (1+F_Si*wt_Si.data+F_S*wt_S.data+F_C*wt.data[k])*Ac.data[k]
            if 'Ca_La' in Quant[k].metadata.Sample.xray_lines: 
                Ac.data[k] = (1+F_Si*wt_Si.data+F_S*wt_S.data+F_C*wt_C.data+F_Ca*wt.data[k])*Ac.data[k]
            if 'O_Ka' in Quant[k].metadata.Sample.xray_lines: 
                Ac.data[k] = (1+F_Si*wt_Si.data+F_S*wt_S.data+F_C*wt_C.data+F_Ca*wt_Ca.data+F_O*wt.data[k])*Ac.data[k]
            if 'Fe_La' in Quant[k].metadata.Sample.xray_lines: 
                Ac.data[k] = (1+F_Si*wt_Si.data+F_S*wt_S.data+F_C*wt_C.data+F_Ca*wt_Ca.data+F_O*wt_O.data+F_Fe*wt.data[k])*Ac.data[k]
    if len(s.data.shape)==2: 
        for i in range(D[0]):
            if navigation_mask.data[i]== False:
                if Dev[i]>Crit/100:
                    for k in range(D[-1]):
                        wt.isig[k] = Quant[k]
                        Ac.inav[i].data[k] = hs.material.mass_absorption_mixture(wt.inav[i].data, elts, energies = s.metadata.Sample.xray_lines[k])    
                        if 'C_Ka' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[i].data[k] = (1+F_Si*wt_Si.inav[i].data+F_S*wt_S.inav[i].data+F_C*wt.inav[i].data[k])*Ac.inav[i].data[k]
                        if 'Ca_La' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[i].data[k] = (1+F_Si*wt_Si.inav[i].data+F_S*wt_S.inav[i].data+F_C*wt_C.inav[i].data+F_Ca*wt.inav[i].data[k])*Ac.inav[i].data[k]
                        if 'O_Ka' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[i].data[k] = (1+F_Si*wt_Si.inav[i].data+F_S*wt_S.inav[i].data+F_C*wt_C.inav[i].data+F_Ca*wt_Ca.inav[i].data+F_O*wt.inav[i].data[k])*Ac.inav[i].data[k]
                        if 'Fe_La' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[i].data[k] = (1+F_Si*wt_Si.inav[i].data+F_S*wt_S.inav[i].data+F_C*wt_C.inav[i].data+F_Ca*wt_Ca.inav[i].data+F_O*wt_O.inav[i].data+F_Fe*wt.inav[i].data[k])*Ac.inav[i].data[k]
    if len(s.data.shape)==3:                             
        for i in range(D[0]):
            for j in range (D[1]):
                if navigation_mask.data[i][j]== False:
                    if Dev[i][j]>Crit/100:
                        for k in range (D[-1]):
                            wt.isig[k] = Quant[k]
                            Ac.inav[j,i].data[k] = hs.material.mass_absorption_mixture(wt.inav[j,i].data, elts, energies = s.metadata.Sample.xray_lines[k])    
                            if 'C_Ka' in Quant[k].metadata.Sample.xray_lines: 
                                Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt.inav[j,i].data[k])#*Ac.inav[j,i].data[k]
                            if 'Ca_La' in Quant[k].metadata.Sample.xray_lines: 
                                Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
                            if 'O_Ka' in Quant[k].metadata.Sample.xray_lines: 
                                Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt_Ca.inav[j,i].data+F_O*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
                            if 'Fe_La' in Quant[k].metadata.Sample.xray_lines: 
                                Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt_Ca.inav[j,i].data+F_O*wt_O.inav[j,i].data+F_Fe*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
    #Calculate the corrected intensities thanks to the abs. correction factors           
    for k in range (len(s.metadata.Sample.xray_lines)):
        result_mod[k] = result_int[k]*Ac.isig[k]/(np.ones((D[:-1]))-np.exp(-Ac.isig[k]*mt.isig[0]/sin(alpha)))
    return result_mod, Ac

def compare (Quant1, Quant2, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D):
    # recalculate the difference between the two lines of interest for next iteration             
    SubQuant = []
    for i in range (len(Quant1)): 
        if line1 in Quant1[i].metadata.Sample.xray_lines : SubQuant.append(Quant1[i])
    for i in range (len(Quant2)):
        if line2 in Quant2[i].metadata.Sample.xray_lines : SubQuant.append(Quant2[i])
        
    rat = SubQuant[0].data/(Elt_rat*SubQuant[1].data)
    rat = np.nan_to_num(rat)
    Dev = (SubQuant[0].data - Elt_rat * SubQuant[1].data)/SubQuant[0].data
    print('line1/line2 deviation (%) =', np.nanmedian(abs(Dev))*100,'  processing next iteration')
    
    if len(D)==1:
        if (mt/(d*10**-7)).data<5: 
            print('pixel ([i]) '+str(i)+' has a thickness lower than 10 nm. Correction is at its minimum, Issue with SNR? k-factors?, fitting model? element ratio?')
    if len(D)==2:
        for i in range(D[0]):
            if navigation_mask.data[i]== False:
                if (mt.inav[i]/(d.inav[i]*10**-7)).data<5: 
                    print('pixel ([i]) '+str(i)+' has a thickness lower than 10 nm. Correction is at its minimum, Issue with SNR? k-factors?, fitting model? element ratio?')
                    #navigation_mask.data[i] = True
                    Dev[i] = np.nan
    if len(D)==3:
        for i in range(D[0]):
            for j in range (D[1]):
                if navigation_mask.data[i][j]== False:
                    if (mt.inav[j,i]/(d.inav[j,i]*10**-7)).data<5: 
                        print('pixel ([i][j]) '+str(i)+'_'+str(j)+' has a thickness lower than 10 nm. Correction is at its minimum, Issue with SNR? k-factors?, fitting model? element ratio?')
                        #navigation_mask.data[i][j] = True
                        Dev[i][j] = np.nan
    return Dev, rat    

def absorption_correction_auto(result, s, kfactors, line1 = 'Fe_Ka', line2 = 'Fe_La', Elt_rat = 1, d = 3, t = 150, tilt_stage = 0, navigation_mask = None, Crit = 1, water = False, valence = None):
    if len(s.data.shape) ==4: raise ValueError('The absorption correction routine is not implemented for 2D signals')  # dimensions are gonna that of a 3D signals, 3 navigation dimension.
    #if t== 0 or d==0 :raise ("thickness and density can't be equal to 0")
    
    D = []
    for i in range (len(s.data.shape)-1):
        D.append(s.data.shape[i])
    D.append(len(s.metadata.Sample.xray_lines))
    
    mt = hs.signals.Signal1D(np.zeros(D[:-1]+[1]))
    if navigation_mask is None:
        navigation_mask = hs.signals.BaseSignal(np.ndarray((D[:-1]), bool)) 
        navigation_mask.data[:] = False
    if type(d)==float or type(d)==int and type(t)==float or type(t)==int:
        mt.isig[0] = d*(t*10**-7)
        t = (mt/(d*10**-7))
        d = mt/(t*10**-7)
    elif type(t)!=(int and float) and type(d)!= (int and float):
        mt.isig[0] = d.isig[0]*t.isig[0]*10**-7
    mt.data[navigation_mask.data] = np.nan

    Dev = np.ones((D[:-1]))
    rat = np.ones((D[:-1]))
        
    s.metadata.Acquisition_instrument.TEM.tilt_stage = tilt_stage
    alpha = (s.metadata.Acquisition_instrument.TEM.Detector.EDS.elevation_angle + 
         s.metadata.Acquisition_instrument.TEM.tilt_stage)*pi/180
                
    elts = []
    for i in range(0, len(s.metadata.Sample.xray_lines)):
        elts.append(result[i].metadata.Sample.elements)   
    print('elts', elts)
    
    a= copy.deepcopy(navigation_mask)
    a.unfold_signal_space()
    l=0
    for i in range (len(a.data)):
        if a.data[i] == False:l=l+1
    print('number of pixels to deal with: ', l, 'total number of pixels: ', len(a.data))        
    
    ### If the two elements used to calculate the mt value are the same, i.e. if ones uses Fe_Ka and Fe_La for instance, then two 
    ## quantifications are made and the intesity of the other line is brought to zero.
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
        
        while (np.nanmedian(abs(Dev)) > Crit/100):
            #Calculation of intensities corrected for absorption
            result_mod_K, Ac_K = correction (s, elts, Quant_K, result_int_K, result_mod_K, alpha, mt, navigation_mask, D, Dev, Crit, Ac_K, wt_K)
            result_mod_L, Ac_L  = correction (s, elts, Quant_L, result_int_L, result_mod_L, alpha, mt, navigation_mask, D, Dev, Crit, Ac_L, wt_L)
            #New quantification using corrected intensities
            Quant_K = s.quantification(method="CL", intensities=result_mod_K, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)
            Quant_L = s.quantification(method="CL", intensities=result_mod_L, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)
            # recalculate the difference between the two lines of interest for next iteration             
            Dev, rat = compare (Quant_K, Quant_L, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D)
            # Calculates a new mass-thickness signal for use in the next Abs.Correction             
            mt.isig[0] = rat*mt.isig[0]
    
        #####################################################################
        ######  Further quantification using WATER corrected Abs Correction factors
        if water == True:
            if valence == None : raise ValueError('a list of valence for each xray_lines or element should be provided as argument')
            Quant_K = hs.material.weight_to_atomic(Quant_K, elements='auto')
            Quant_K, H2O = water_content(s, Quant_K, valence)
            Quant_K = hs.material.atomic_to_weight(Quant_K, elements='auto')

            Quant_L = hs.material.weight_to_atomic(Quant_L, elements='auto')
            Quant_L, H2O = water_content(s, Quant_L, valence)
            Quant_L = hs.material.atomic_to_weight(Quant_L, elements='auto')

            """Dev, rat = compare (Quant_K, Quant_L, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D)
            mt.isig[0] = rat*mt.isig[0]"""
            while (np.nanmedian(abs(Dev)) > Crit/100):
                result_mod_K, Ac_K = correction (s, elts, Quant_K, result_int_K, result_mod_K, alpha, mt, navigation_mask, D, Dev, Crit, Ac_K, wt_K)
                result_mod_L, Ac_L = correction (s, elts, Quant_L, result_int_L, result_mod_L, alpha, mt, navigation_mask, D, Dev, Crit, Ac_L, wt_L)

                #New quantification using corrected Abs Correction factors
                Quant_K = s.quantification(method="CL", intensities=result_mod_K, factors=kfactors, composition_units='atomic', navigation_mask = navigation_mask, plot_result=False)
                Quant_K, H2O = water_content (s, Quant_K, valence)
                Quant_K = hs.material.atomic_to_weight(s, Quant_K, elements='auto')

                Quant_L = s.quantification(method="CL", intensities=result_mod_L, factors=kfactors, composition_units='atomic', navigation_mask = navigation_mask, plot_result=False)
                Quant_L, H2O = water_content (Quant_L, valence)
                Quant_L = hs.material.atomic_to_weight(Quant_L, elements='auto')

                Dev, rat = compare (Quant_K, Quant_L, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D)
                mt.isig[0] = rat*mt.isig[0]
            print('Water computed')
        Quant3 = hs.material.weight_to_atomic(Quant_K, elements = 'auto')

###########################################################################################################        
### If two different elements are used, then only one quantification is performed and the convergence criteria are internal to this calculation
    if line1[:2] != line2[:2]:      
        result_int = copy.deepcopy(result)
        result_mod = copy.deepcopy(result)
        Ac = hs.signals.Signal1D(np.ones((D)))
        wt = hs.signals.Signal1D(np.ones((D)))
        Quant2 = s.quantification(method="CL", intensities=result_int, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)
              
        while (np.nanmedian(abs(Dev)) > Crit/100):                  
            result_mod, Ac = correction (s, elts, Quant2, result_int, result_mod, alpha, mt, navigation_mask, D, Dev, Crit, Ac, wt)
    
            #New quantification using corrected Abs Correction factors
            Quant2 = s.quantification(method="CL", intensities=result_mod, factors=kfactors, composition_units='weight', navigation_mask = navigation_mask, plot_result=False)

            # recalculate the difference between the two lines of interest for next iteration             
            Dev, rat = compare (Quant2, Quant2, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D)
            mt.isig[0] = rat*mt.isig[0]

        ###########################################
        #New quantification using WATER corrected Abs Correction factors
        ############################################
        if water == True:
            if valence == None : raise ValueError('a list of valence for each xray_lines or element should be provided as argument')
            Quant2 = hs.material.weight_to_atomic(Quant2, elements='auto')
            Quant2, H2O = water_content (s, Quant2, valence)
            Quant2 = hs.material.atomic_to_weight(Quant2, elements='auto')

            while (np.nanmedian(abs(Dev)) > Crit/100):
                result_mod, Ac = correction (s, elts, Quant2, result_int, result_mod, alpha, mt, navigation_mask, D, Dev, Crit, Ac, wt)

                #New quantification using corrected Abs Correction factors
                Quant2 = s.quantification(method="CL", intensities=result_mod, factors=kfactors, composition_units='atomic', navigation_mask = navigation_mask, plot_result=False)
                Quant2, H2O = water_content (s, Quant2, valence)
                Quant2 = hs.material.atomic_to_weight(Quant2, elements='auto')

                # recalculate the difference between the two lines of interest for next iteration             
                Dev, rat = compare (Quant2, Quant2, Dev, rat, line1, line2, Elt_rat, mt, d, navigation_mask, D)
                mt.isig[0] = rat*mt.isig[0]
            print('Water computed')
        Quant3 = hs.material.weight_to_atomic(Quant2, elements='auto')
    
    if water==False : H2O = np.zeros((D[:-1]))
    if len(s.data.shape)==2:                             
        for i in range(D[0]):
            if navigation_mask.data[i]== True:
                mt.inav[i].data = np.nan
                H2O[i] = np.nan
                for k in range (len(Quant3)):
                    Quant3[k].data[i]=np.nan
    if len(s.data.shape)==3:                             
        for i in range(D[0]):
            for j in range (D[1]):
                if navigation_mask.data[i][j]== True:
                    mt.inav[j,i].data = np.nan
                    H2O[i][j] = np.nan
                    for k in range (len(Quant3)):
                        Quant3[k].data[i][j]=np.nan
    if water == True: return Quant3, H2O, mt, Dev
    if water == False: return Quant3, mt, Dev