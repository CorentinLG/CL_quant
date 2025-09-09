# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:10:00 2020

@author: Corentin
"""

from math import sin
from math import pi
from Correction import correction
import numpy as np
import copy
import hyperspy.api as hs
import exspy

def absorption_correction (result, s, kfactors, d = 3, t = 100, tilt_stage = 0, navigation_mask = None, Crit = 0.5, composition_units='atomic'):
    """
    This function quantifies TEM-EDS data accounting for absorption (a function of sample thickness, 
    density, and take-off angle) based on the Cliff-Lorimer method. It is based on the hyperspy quantification function. 
    Mass attenuation coefficients are from Chantler et al. (2005).
    This is an iterative process : quantification is ran, then absorption is calculated based on this first quantification, and then an updated quantification is retrieved. 
    This sequence is performed until two successive quantifications are similar within 0.5% (by default). 
	For calculation principles, See William, D. B., & Carter, C. B. (1996). Transmission electron microscopy: a textbook for materials science. Edit. Plenum Press, New York & London.]    
    CAUTION : x-ray lines are taken from the signal metadata. There must be only one family line (K, L or M) per element.
    
    Parameters
    ----------
    result : a list of signals, the integrated gaussian intensitiy of each X-ray lines.
	s: the signal being quantified. Its metadata (elements and xray_lines) must be documented.
    kfactors : a list of k-factors (must have the same length as xray_lines from s.metadata and as intensities)
    d : density in g.cm**3 (can be a single value or have the same navigation dimension as the dataset. In this case a density map is required as well.
    t : thickness is in nm (can be a single value or have the same navigation dimension as the dataset, but must have the same type and dimension as d).
    tilt_stage : the sample holder tilt stage, used to calculate the absorption pathways (in degrees).
    navigation_mask : None By default, a BaseSignal of booleans otherwise.
    Crit  : The convergence criterium in percentage of the difference between 2 successive iteration. 0.5% by default.
	composition_units : the output units of the quantification

    Returns
    ------
    A list of quantified elemental composition of the sample (in atomic percent)
    taking into account the absorption correction.

    Examples
    --------
    
    
    """
    
    if len(s.data.shape) ==4: raise ValueError('The absorption correction routine is not implemented for 2D signals')  # dimensions are gonna that of a 3D signals, 3 navigation dimension.
    #if t== 0 or d==0 :raise ("thickness and density can't be equal to 0")
    
    ### First step is to define the dimension of the dataset (single spectrum, line profile or 2D map).
	# D is a variable containing the navigation dimension and the number of x-ray lines
	
    D = []
    for i in range (len(s.data.shape)-1):
        D.append(s.data.shape[i])
    D.append(len(s.metadata.Sample.xray_lines))
    
	### To mask part of a dataset, a navigation_mask can be provided. If None, then a default one (containing all pixels) is created.
    if navigation_mask is None:
        navigation_mask = hs.signals.BaseSignal(np.ndarray((D[:-1]), bool)) 
        navigation_mask.data[:] = False
    ### To be able to deal with a sample of variable thickness and density, one can provide a density or a thickness map. They are combined here under the 'mt' variable.
 	# If a single value is given as inputs, then it is converted to a signal that has the same dimension as the raw data. 
	
    if type(d)==float or type(d)==int and type(t)==float or type(t)==int:
        mt = hs.signals.Signal1D(np.zeros(D[:-1]+[1]))
        mt.isig[0] = d*(t*10**-7)
        mt.data[navigation_mask.data] = np.nan
    elif type(t)!=(int and float) and type(d)!= (int and float):
            mt = d*t*10**-7
    
	### dif is a variable constructed as the difference between the composition for each X-ray lines at each position from two successive iteration.
	# It is used to check if convergence as been reached.

    dif = np.ones((D[-1:]+D[:-1]))
    Dev = np.ones((D[:-1]))

	### Absorption depends on the take off angle, which is calculated here. (Actually should be better to use the take_off_angle function from hyperspy, which requiers the metadata to be carefully written)
    s.metadata.Acquisition_instrument.TEM.tilt_stage = tilt_stage
    alpha = (s.metadata.Acquisition_instrument.TEM.Detector.EDS.elevation_angle + 
         s.metadata.Acquisition_instrument.TEM.tilt_stage)*pi/180 
    
	### The function later uses the mass_absorption_mixture function from HS, which requires a list of elements.
    elts = []
    for i in range(len(s.metadata.Sample.xray_lines)):
        elts.append(result[i].metadata.Sample.elements)		
    
    if len(s.data.shape)!=1:
        if len(s.data.shape)==2: 
            a= D[0]     
        if len(s.data.shape)==3: 
            a= D[0]*D[1]
        l = a-np.count_nonzero(navigation_mask.data)
        print('number of pixels to deal with: ', l, 'total number of pixels: ', a)       

	### since the "result" variable is manipulated several times, better copy it deeply before.
    result_int = copy.deepcopy(result) 
    result_mod = copy.deepcopy(result)   
    
	### Ac contains the absorption correction factor for each X-ray lines at each position, and wt is the quantified composition at each pixels.
	
    Ac = hs.signals.Signal1D(np.zeros((D)))
    wt = hs.signals.Signal1D(np.zeros((D)))
	
	### A first quantification is ran without correction, to get things started. Quant2 is needed for later comparison of successive iterations.
	# Given that the mass_absorption_mixture function requires weight percent, this is the unit of choice all along.
    Quant = s.quantification(method="CL", intensities=result_int, factors=kfactors, composition_units='weight', navigation_mask=navigation_mask, plot_result=False)   
    Quant2 = Quant
	

    while (abs(dif) > Crit/100).any():
        Quant = Quant2
        
        #Calculation of intensities corrected for absorption
        result_mod, Ac = correction(s, elts, Quant, result_int, result_mod, alpha, mt, navigation_mask, D, Dev, dif, 0.0001, Ac, wt)
        
        #New quantification using corrected intensities
        Quant2 = s.quantification(method="CL", intensities=result_mod, factors=kfactors, composition_units='weight', navigation_mask=navigation_mask, plot_result=False)          
          
        #Compares the relative difference between previous and last quantification (convergence test value)
        for i in range (len(Quant)):
            dif[i]= abs(Quant[i].data-Quant2[i].data)/Quant[i].data

    if composition_units == 'atomic': Quant3 = exspy.material.weight_to_atomic(Quant2, elements='auto')
    if composition_units == 'weight': Quant3 = Quant2

    return Quant3
