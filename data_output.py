# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 11:47:28 2018

@author: Corentin
"""

import xlsxwriter
import numpy as np

def data_output (Quant, result, H2O_List = None, mt_List = None, density_or_thickness = None, name = 'data.xlsx'):    #H2O_List, mt_List, density_or_thickness, Dev
    
    workbook = xlsxwriter.Workbook(name, {'nan_inf_to_errors': True})
    worksheet_table = workbook.add_worksheet(name = 'Quantifed data')
    worksheet_error = workbook.add_worksheet(name = 'Errors')
    
    Q = np.ndarray((len(Quant), len((Quant[0].data).flatten())))
    for i in range (len (Quant)):
        Q[i] = (Quant[i].data).flatten()
        Q[i] = np.nan_to_num(Q[i])

    R = np.ndarray((len(result), len((result[0].data).flatten())))
    for i in range (len (result)):
        R[i] = (result[i].data).flatten()
        R[i] = np.nan_to_num(R[i])
    
    if mt_List is not None : 
        mt = mt_List.data.flatten()
        mt = np.nan_to_num(mt)
    if H2O_List is not None : 
        H2O = H2O_List.flatten()
        H2O = np.nan_to_num(H2O)
    if density_or_thickness is not None:
        d_or_t = np.ones(shape=mt.shape, dtype=float)
        for i in range (len(d_or_t)):
            d_or_t[i]= density_or_thickness

    row = []
    for i in range (len(Quant)):
        row.append(0)
        row [i] = str(Quant[i].metadata.Sample.xray_lines)
    worksheet_table.write_row ('B1', row)
    worksheet_error.write_row ('B1', row)
    worksheet_table.write (0, len(row)+1, 'H2O')
    worksheet_table.write(0, len(row)+2, 'Deduced t or d')
    worksheet_table.write(0, len(row)+3, 'Chosen t or d')
    worksheet_table.write(0, len(row)+4, 'Auto Abs . Cor. deviation')
    
    for i in range(len(Q)):
        k=0
        for j in range (len(Q[0])):
            if Q[:,j].any() != 0:
                worksheet_table.write (k+1,0, str(k))
                worksheet_table.write(k+1,i+1, Q[i][j])
                if H2O_List is not None : worksheet_table.write (k+1, len(row)+1, H2O[j])
                if mt_List is not None :worksheet_table.write (k+1, len(row)+2, mt[j]/(d_or_t[j]*10**-7))
                if density_or_thickness is not None : worksheet_table.write (k+1, len(row)+3, d_or_t[j])
                #worksheet_table.write (i+1, len(row)+4, Dev[i])
                k=k+1
    for i in range(len(R)):
        k=0
        for j in range (len(R[0])):
            if Q[:,j].any() != 0 and R[i][j]!= 0 and Q[i][j]!= np.isnan:
                worksheet_error.write (k+1,0, str(k))
                worksheet_error.write(k+1,i+1, (Q[i][j]*R[i][j]**0.5)/(R[i][j]))
                k=k+1
        
    workbook.close()