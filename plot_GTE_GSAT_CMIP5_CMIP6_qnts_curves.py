#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 12:13:40 2020

This script plots annual mean ensemble quantiles of GTE and GSAT time series 
from CMIP5 and CMIP6 data. 

Input parameters:
    cmip6_path        path containing CMIP6 GSAT & GTE ensemble NetCDFs
    cmip5_path        path containing CMIP5 GSAT & GTE ensemble NetCDFs

@author: Tim Hermans
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

cmip6_path = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs'
cmip5_path = '/Users/thermans/Documents/Data/CMIP6/GMSL_MonteCarlo_JG/input/CMIP5/'

RCPtoSSP = { #dictionary linking ssps and rcps
	"rcp26": "ssp126",
	"rcp45": "ssp245",
    "rcp85": "ssp585",
    }        

#read cmip6 data
cmip6_gsat = xr.open_dataset(os.path.join(cmip6_path,'tas_CMIP6_n17_1986_2005ref_1850_2100_am.nc'),decode_times=True)
cmip6_gte = xr.open_dataset(os.path.join(cmip6_path,'zostoga_CMIP6_n17_1986_2005ref_qdedr_1850_2100_am.nc'),decode_times=True)

#calculate quantiles
cmip6_gsat_qnts = cmip6_gsat.tas.quantile([.05,.5,.95],dim='model',skipna=True)  
cmip6_gte_qnts = cmip6_gte.zostoga.quantile([.05,.5,.95],dim='model',skipna=True)  

#initiate plot
fig,axs = plt.subplots(3,2)
           
colors = ['C1','C0','C2']

for r,rcp in enumerate(reversed(list(RCPtoSSP.keys()))): #from high to low scenario
    #open CMIP5 GTE (stored per scenario)
    cmip5_gte = xr.open_dataset(os.path.join(cmip5_path,rcp+'_zostoga.nc'),decode_times=True)
    cmip5_gte_qnts = cmip5_gte.global_average_thermosteric_sea_level_change.quantile([.05,.5,.95],dim='region',skipna=True) 

    #plot CMIP5 GTE
    axs[r,0].plot(cmip5_gte['time.year'],cmip5_gte_qnts.isel(quantile=0),color='black',linestyle='dashed',linewidth=1)
    axs[r,0].plot(cmip5_gte['time.year'],cmip5_gte_qnts.isel(quantile=1),color='black',label='CMIP5 (n=21)')
    axs[r,0].plot(cmip5_gte['time.year'],cmip5_gte_qnts.isel(quantile=2),color='black',linestyle='dashed',linewidth=1)
    
    #plot CMIP6 GTE
    cmip6_gte_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=1).plot(ax=axs[r,0],label='CMIP6 (n='+str(len(cmip6_gte.model))+')',color=colors[r]) #median
    axs[r,0].fill_between(cmip6_gte_qnts.year, cmip6_gte_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=0), cmip6_gte_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=2),color=colors[r] ,alpha=0.3) #uncertainty range

    #open CMIP5 GSAT (stored per scenario)
    cmip5_gsat = xr.open_dataset(os.path.join(cmip5_path,rcp+'_tas.nc'),decode_times=True)
    cmip5_gsat_qnts = cmip5_gsat.air_temperature.quantile([.05,.5,.95],dim='region',skipna=True) 

    #plot CMIP5 GSAT
    axs[r,1].plot(cmip5_gsat['time.year'],cmip5_gsat_qnts.isel(quantile=0),color='black',linestyle='dashed',linewidth=1)
    axs[r,1].plot(cmip5_gsat['time.year'],cmip5_gsat_qnts.isel(quantile=1),color='black',label='CMIP5 (n=21)')
    axs[r,1].plot(cmip5_gsat['time.year'],cmip5_gsat_qnts.isel(quantile=2),color='black',linestyle='dashed',linewidth=1)
    
    #plot CMIP6 GSAT
    cmip6_gsat_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=1).plot(ax=axs[r,1],label='CMIP6 (n='+str(len(cmip6_gte.model))+')',color=colors[r]) #median
    axs[r,1].fill_between(cmip6_gsat_qnts.year, cmip6_gsat_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=0), cmip6_gsat_qnts.sel(scen=RCPtoSSP[rcp]).isel(quantile=2),color=colors[r] ,alpha=0.3) #uncertainty range
    
    #set labels & limits
    axs[r,0].set_ylabel('GTE [m]') 
    axs[r,0].set_xlabel('')
    axs[r,0].set_xlim([2005,2100])
    axs[r,0].set_ylim([0,.45])
    axs[r,0].grid()
    axs[r,0].legend(loc='upper left')
    axs[r,0].set_title('('+'abc'[r]+') SSP'+RCPtoSSP[rcp][-3]+'-'+'RCP'+rcp[-2]+'.'+rcp[-1])
    
    axs[r,1].set_ylabel('GSAT [K]')
    axs[r,1].set_xlim([2005,2100])
    axs[r,1].set_ylim([0,7])
    axs[r,1].grid(); 
    #axs[r,1].legend(loc='upper left')
    axs[r,1].set_title('('+'def'[r]+') SSP'+RCPtoSSP[rcp][-3]+'-'+'RCP'+rcp[-2]+'.'+rcp[-1])

    if r==2:
        axs[r,0].set_xlabel('year')
        axs[r,1].set_xlabel('year')
    else:
        axs[r,0].set_xlabel('')
        axs[r,1].set_xlabel('')
        
        axs[r,0].set_xticklabels('')
        axs[r,1].set_xticklabels('')

    
    