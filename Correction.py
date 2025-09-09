from math import sin
import numpy as np
import hyperspy.api as hs
import exspy 

def correction(s, elts, Quant, result_int, result_mod, alpha, mt, Full_mask, D, Dev, dif, Crit_mt, Ac, wt):
       
    wt_Si, wt_S, wt_C, wt_Ca, wt_O = wt.isig[0], wt.isig[0], wt.isig[0], wt.isig[0], wt.isig[0]


    if len(s.data.shape)==1: 
        for k in range(D[-1]):
            wt.isig[k] = Quant[k]
            Ac.data[k] = exspy.material.mass_absorption_mixture(wt.data, elts, energies = s.metadata.Sample.xray_lines[k])    
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k].transpose()/(np.ones((D[:-1]))-np.exp(-Ac.isig[k].transpose()*mt.isig[0].transpose()/sin(alpha)))

    if len(s.data.shape)==2: 
        for i in range(D[0]):
            if Full_mask.data[i] == False:
                for k in range(D[-1]):
                    wt.isig[k] = Quant[k]
                    Ac.inav[i].data[k] = exspy.material.mass_absorption_mixture(wt.inav[i].data, elts, energies = s.metadata.Sample.xray_lines[k])    
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k]/(np.ones((D[:-1]))-np.exp(-Ac.isig[k]*mt.isig[0]/sin(alpha)))
            
    if len(s.data.shape)==3:                             
        for i in range(D[0]):
            for j in range (D[1]):
                if Full_mask.data[i][j] == False:
                    for k in range (D[-1]):
                        wt.isig[k] = Quant[k]
                        Ac.inav[j,i].data[k] = exspy.material.mass_absorption_mixture(wt.inav[j,i].data, elts, energies = s.metadata.Sample.xray_lines[k])    
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k]/(np.ones((D[:-1]))-np.exp(-Ac.isig[k]*mt.isig[0]/sin(alpha)))
            
    return result_mod, Ac
