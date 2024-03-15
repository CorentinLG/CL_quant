# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 16:26:24 2017

@author: Corentin
"""

import hyperspy.api as hs

def water_content (s, quantification, valence):
    """This routine calculate the water content of a quantified spectrum based on the electroneutrality assumption,
    and add H to the list of elements.
    
    quantification : the output of the quantification routine
    valence : water calculation requires the list of valences of each element to determine the charge balance. 
    It should  be a list of values of the same length as the quantification output. Ex : for Fe3Si2O5(OH)4 if xray_lines
    are 'Fe_Ka', 'Si_Ka' and 'O_Ka', valence  = [2, 4, -2].
    """    
    
    import copy

    Quant_H2O = copy.deepcopy(quantification)
    H = copy.deepcopy(Quant_H2O[0])      # H (at.%) correspond to the charge excess, if any. 
    H.metadata = copy.deepcopy(Quant_H2O[0].metadata)
    H.metadata.General.title = 'atomic percent of H'
    H.metadata.Sample.elements = ['H']
    H.metadata.Sample.xray_lines = ['H_Ka']
    
    D = []
    for i in range (len(s.data.shape)-1):
        D.append(s.data.shape[i])
    D.append(len(s.metadata.Sample.xray_lines))
    
    if len(s.data.shape)==1: H.data = 0
    if len(s.data.shape)==2: 
        for i in range (D[0]):
            H.data[i] = 0
    if len(s.data.shape)==3: 
        for i in range (D[0]):
            for j in range (D[1]):
                H.data[i][j] = 0   
        
    for i in range (len(Quant_H2O)):
        H = H + Quant_H2O[i]*(-1)*valence[i]

    #Given that the totals are being modified, new at.% must be recalcutated for all elements.
    for i in range (len(Quant_H2O)):
        Quant_H2O[i] = (Quant_H2O[i]*100)/(H+100)
    H = (H*100)/(H+100)

    Quant_H2O.append(H)
    Quant_H2O_wt = hs.material.atomic_to_weight(Quant_H2O)
    
    H2O = []
    for i in range (len(Quant_H2O)):
        if 'H_Ka' in Quant_H2O[i].metadata.Sample.xray_lines: H2O = Quant_H2O_wt[i].data*9

    return Quant_H2O, H2O