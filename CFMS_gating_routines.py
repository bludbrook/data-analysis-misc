# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:52:07 2016

@author: ludbroba

A collection of functions to load, process, analyze and plot Hall effect measurements.
"""
import pandas as pd 
import os
import glob
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
plt.style.use('ggplot')

"""
Define some functions to import data
"""
def load_file_press(file_str,Vg,t):
    """ Load data measured using a press-puck
        arguments:  file_str    the file name
                    Vg          applied gate voltage
                    t           sample thickness (nm) 
    """
    
    file_name = file_str
    thick = t
    gate = Vg

    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21),skip_header=32,names=None)
    c_names = ['Temp1','Field','Pos','Resistance','AHE']
    
#    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21,34),skip_header=32,names=None)
#    c_names = ['Temp1','Field','Pos','Resistance','AHE','Temp2']
    df1 = pd.DataFrame(data[:].copy(), columns=c_names)
    df1['Thick'] = thick
    df1['Gate'] = gate
    h_max = df1['AHE'].max()
    h_min = df1['AHE'].min()
    shift = h_min + (h_max-h_min)/2
    df1['AHE'] = df1['AHE']*(-1) + shift
    return df1

def load_file_press_T2(file_str,Vg,t):
""" Load data measured using a press-puck with T2 activated
        arguments:  file_str    the file name
                    Vg          applied gate voltage
                    t           sample thickness (nm) 
    """
    
    file_name = file_str
    thick = t
    gate = Vg

    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21,34),skip_header=32,names=None)
    c_names = ['Temp1','Field','Pos','Resistance','AHE','Temp2']
    df1 = pd.DataFrame(data[:].copy(), columns=c_names)
    df1['Thick'] = thick
    df1['Gate'] = gate
    h_max = df1['AHE'].max()
    h_min = df1['AHE'].min()
    shift = h_min + (h_max-h_min)/2
    df1['AHE'] = df1['AHE']*(-1) + shift
    return df1

def load_file(file_str,Vg,t):
    """ Load data measured using a Hall-bar
        arguments:  file_str    the file name
                    Vg          applied gate voltage
                    t           sample thickness (nm) 
    """
    
    file_name = file_str
    thick = t
    gate = Vg

    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21),skip_header=32,names=None)
    c_names = ['Temp1','Field','Pos','AHE','Resistance']
    
#    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21,34),skip_header=32,names=None)
#    c_names = ['Temp1','Field','Pos','Resistance','AHE','Temp2']
    df1 = pd.DataFrame(data[:].copy(), columns=c_names)
    df1['Thick'] = thick
    df1['Gate'] = gate
    h_max = df1['AHE'].max()
    h_min = df1['AHE'].min()
    shift = h_min + (h_max-h_min)/2
    df1['AHE'] = df1['AHE']*(-1) + shift
    return df1

def load_file_T2(file_str,Vg,t):
    """ Load data measured using a Hall-bar with T2 activated
        arguments:  file_str    the file name
                    Vg          applied gate voltage
                    t           sample thickness (nm) 
    """
    
    file_name = file_str
    thick = t
    gate = Vg

    data =  np.genfromtxt(file_name,delimiter=',',usecols=(3,4,5,20,21,34),skip_header=32,names=None)
    c_names = ['Temp1','Field','Pos','Resistance','AHE','Temp2']
    df1 = pd.DataFrame(data[:].copy(), columns=c_names)
    df1['Thick'] = thick
    df1['Gate'] = gate
    h_max = df1['AHE'].max()
    h_min = df1['AHE'].min()
    shift = h_min + (h_max-h_min)/2
    df1['AHE'] = df1['AHE']*(-1) + shift
    return df1
    

"""
Some functions to process the data
"""

def procT_df(df_name):
    """
    Split data measured at different temperatures into separate dataframes
    Correct for arbitrary data offset
    """
    df1 = df_name
    df1['Temp2_diff']=df1['Temp2'].diff()
    df_filt = df1[df1['Temp2_diff'].abs() > 4]

    ends = df_filt.index.tolist()
    start_0 = 0
    ends.insert(0,start_0)
    ends.append(len(df1))
    dataframes = {} # Create a dict to store all the dataframes
    #
    for i in range(len(df_filt)+1):
        start = ends[i]
        end = ends[i+1]-1
        dataframes[i] = df1.ix[start:end,:]

    corr_data = {}

    for i in range(len(df_filt)+1):
        df = dataframes[i].copy()
        
        df['Pos_diff']=df['Pos'].diff()
        df_filt1 = df[df['Pos_diff'] > 10]
        flip = df_filt1.index.tolist()[0]-1
        dfa = df.ix[:flip,:]
        dfc = dfa.copy()
        h_max = dfc['AHE'].max()
        h_min = dfc['AHE'].min()
        shift = h_min + (h_max-h_min)/2
        dfc['AHE'] = dfc['AHE']*(-1) + shift

        corr_data[i]=dfc
    return corr_data

def procT_NF_df(df_name):
     """
    Split data measured at different temperatures into separate dataframes
    Correct for arbitrary data offset
    """
    df1 = df_name
    df1['Temp2_diff']=df1['Temp2'].diff()
    df_filt = df1[df1['Temp2_diff'].abs() > 2]

    ends = df_filt.index.tolist()
    start_0 = 0
    ends.insert(0,start_0)
    ends.append(len(df1))
    dataframes = {} # Create a dict to store all the dataframes
    #
    for i in range(len(df_filt)+1):
        start = ends[i]
        end = ends[i+1]-1
        dataframes[i] = df1.ix[start:end,:]

    corr_data = {}

    for i in range(len(df_filt)+1):
        df = dataframes[i].copy()
        
        dfc = df.copy()
        h_max = dfc['AHE'].max()
        h_min = dfc['AHE'].min()
        shift = h_min + (h_max-h_min)/2
        dfc['AHE'] = dfc['AHE']*(-1) + shift

        corr_data[i]=dfc
    return corr_data
    
def split_T(df_im):
     """
    Split data measured at different temperatures into separate dataframes
    Correct for arbitrary data offset
    """
    df = df_im.copy()
    df['Temp2_diff']=df['Temp2'].diff()
    df_filt = df[df['Temp2_diff'].abs() > 2]
    
    ends = df_filt.index.tolist()
    start_0 = 0
    ends.insert(0,start_0)
    ends.append(len(df))
    dataframes = {} # Create a dict to store all the dataframes
    
    for i in range(len(df_filt)+1):
        start = ends[i]
        end = ends[i+1]-1
        dataframes[i] = df.ix[start:end,:]
    
    corr_data = {} # Create a dict to store the corrected dataframes
    
    for i in range(len(df_filt)+1):
        df1 = dataframes[i].copy()
    
        dfc = df1.copy()
        h_max = dfc['AHE'].max()
        h_min = dfc['AHE'].min()
        shift = h_min + (h_max-h_min)/2
        dfc['AHE'] = dfc['AHE']*(-1) + shift
    
        corr_data[i]=dfc
    print('Dataframe has been split by temperature. Data is stored in a dict')
    return corr_data
    
def analyse_T(dict_sample):
     """
    Perform analysis on Hall data to extract key parameters
    Takes a dict of temperature-split data
    """
    corr_data = dict_sample
    Thickness = []
    Resistivity = []
    R_eff = []
    C_conc = []
    C_sheet = []
    AHE = []
    AHE_fit = []
    Temperature = []
    H_coercive = []
    H_c1 = []
    H_c2 = []
    Slope = []
    Rem = []
    
    
    for temp, df1 in corr_data.items():
    
        thick = df1['Thick'].iloc[0]
        temperature = df1['Temp2'].iloc[1]
        resistance = df1['Resistance'].mean()
        resistivity = df1['Resistivity'].mean()
        r_eff = df1['Resistivity_eff'].mean()
        #================================================
        # Get the normal Hall slope
        #================================================        
        df_aa = df1.head(50)
        df_a = df_aa[(df_aa['Field'] > 5000)]
        x1 = df_a['Field'].values
        y1 = df_a['AHE'].values
        fit1 = np.polyfit(x1,y1,1)
        polynomial1 = np.poly1d(fit1)
        y1s = polynomial1(x1)
        df_bb = df1.tail(50)
        df_b = df_bb[(df_bb['Field'] < -5000)]
        x2 = df_b['Field'].values
        y2 = df_b['AHE'].values
        fit2 = np.polyfit(x2,y2,1)
        polynomial2 = np.poly1d(fit2)
        y2s = polynomial2(x2)
        slope = (fit1[0]+fit2[0])/2
        #================================================
        # Calculate the carrier concentration
        #================================================      
        q = 1.6E-19
        const = 1E-8
        n_Pd = 1E21
        t_Pd = 2.5E-7
        RH_Pd = 1/(q*n_Pd)
        rho_Pd = 60E-6
        t_Tot = thick*(1E-7)+t_Pd
        rho_Tot = resistance*t_Tot
    
        c1 = RH_Pd*((1/rho_Pd)**2)*t_Pd
        t1 = slope*t_Tot*((1/rho_Tot)**2)*t_Tot/const
        t2 = ((1/resistivity)**2)*thick*(1E-7)
    
        RH = (t1-c1)/t2
        n = -1/(q*RH) 
        
        n_s = -1/(q*(slope/const))
        
        #================================================
        # Get the AHE        
        #================================================
        x=0
        df2 = df1.iloc[(df1.Field-x).abs().argsort()[:4]]
        ahe = df2['AHE'].abs().mean()
        #        
        sign = df1.iloc[-2]['AHE']
        sign = sign/abs(sign)
        ahe=ahe*sign
        #================================================
        # Get the AHE from the fits       
        #================================================
        ahe_fit = (abs(fit2[1])+abs(fit2[1]))/2
        #================================================
        # Get the Hc     
        #================================================
        df1['AHE_diff']=df1['AHE'].diff().abs()
        df_hc1 = df1.sort_values('AHE_diff',ascending=False).head(2)
        Hc1 = df_hc1.Field.abs().mean()
    
        x=0
        df_hc2 = df1.iloc[(df1.AHE-x).abs().argsort()[:4]]
        Hc2 = df_hc2['Field'].abs().mean()
        
        Hc = (Hc1+Hc2)/2
        #================================================
        # Get the Remnance     
        #================================================
        x=0
        df_rem = df1.iloc[(df1.Field-x).abs().argsort()[:2]]
        rem = df_rem['AHE'].abs().mean()
        remnance = rem/ahe_fit
        #================================================
        # Compile all the results       
        #================================================
        C_conc.append(n)        
        Thickness.append(thick)
        Resistivity.append(resistivity)
        Temperature.append(temperature)
        AHE.append(ahe)
        AHE_fit.append(ahe_fit)
        H_coercive.append(Hc)
        H_c1.append(Hc1)
        H_c2.append(Hc2)
        Slope.append(slope)
        R_eff.append(r_eff)
        C_sheet.append(n_s)
        Rem.append(remnance)
    data_all = pd.DataFrame({'Temperature': Temperature,'R_eff':R_eff, 
                             'Thickness':Thickness, 'Resistivity':Resistivity, 
                             'AHE': AHE, 'AHE_fit': AHE_fit, 'C_conc':C_conc,
                             'C_sheet':C_sheet, 'H_c':H_coercive, 'Slope':Slope,
                             'H_c1':H_c1, 'H_c2':H_c2, 'Rem':Rem
                            })
    return data_all
