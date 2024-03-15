from math import sin
import numpy as np
import hyperspy.api as hs

def correction(s, elts, Quant, result_int, result_mod, alpha, mt, Full_mask, D, Dev, dif, Crit_mt, Ac, wt):
    F_Si = 0
    F_S = 0.01
    F_C = 0  #0.001
    F_Ca = 0.01 # 0.01
    F_O = 0.00
    F_Fe = 0.024
       
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
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k].transpose()/(np.ones((D[:-1]))-np.exp(-Ac.isig[k].transpose()*mt.isig[0].transpose()/sin(alpha)))

    if len(s.data.shape)==2: 
        for i in range(D[0]):
            if Full_mask.data[i] == False:
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
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k]/(np.ones((D[:-1]))-np.exp(-Ac.isig[k]*mt.isig[0]/sin(alpha)))
    if len(s.data.shape)==3:                             
        for i in range(D[0]):
            for j in range (D[1]):
                if Full_mask.data[i][j] == False:
                    for k in range (D[-1]):
                        wt.isig[k] = Quant[k]
                        Ac.inav[j,i].data[k] = hs.material.mass_absorption_mixture(wt.inav[j,i].data, elts, energies = s.metadata.Sample.xray_lines[k])    
                        if 'C_Ka' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
                        if 'Ca_Ka' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
                        if 'O_Ka' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt_Ca.inav[j,i].data+F_O*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
                        if 'Fe_La' in Quant[k].metadata.Sample.xray_lines: 
                            Ac.inav[j,i].data[k] = (1+F_Si*wt_Si.inav[j,i].data+F_S*wt_S.inav[j,i].data+F_C*wt_C.inav[j,i].data+F_Ca*wt_Ca.inav[j,i].data+F_O*wt_O.inav[j,i].data+F_Fe*wt.inav[j,i].data[k])*Ac.inav[j,i].data[k]
        #Calculate the corrected intensities thanks to the abs. correction factors           
        for k in range (len(s.metadata.Sample.xray_lines)):
            result_mod[k] = result_int[k]*Ac.isig[k]/(np.ones((D[:-1]))-np.exp(-Ac.isig[k]*mt.isig[0]/sin(alpha)))
    return result_mod, Ac