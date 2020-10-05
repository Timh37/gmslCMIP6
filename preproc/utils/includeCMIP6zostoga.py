#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function incorporates dedrifted 'zostoga' data of CMIP6 models in one Xarray Dataset.
It looks for the first variant and grid member of each model for which both
scenario and historical run are available. The monthly means are converted to 
annual means and the historical and future period are combined along the time axis.

input: 
targ_years          timeseries years wished to incorporate
scenarios           ssps wished to incorporate
zostoga_path        path to dedrifted zostoga data stored per variable per model

output:
zostoga_ds_am       Dataset with zostoga as function of year, model and scenario
    zostoga             -Zostoga data of model for this scenario (nan if model not incorporated)
    variant_label       -Variant label of model for this scenario (none if model not incorporated)
    grid_label          -Grid label of model for this scenario (none if model not incorporated)
    incorporated        -If model data is available for this scenario (True or False)
    
Created on Fri May 15 15:08:40 2020

@author: thermans
"""
import xarray as xr
import numpy as np
import os
import fnmatch

def includeCMIP6zostoga(zostoga_path,scenarios,targ_years):

    models = [model for model in os.listdir(zostoga_path) if not '.' in model] #get list of model directories
    
    #initialize empty data arrays
    data = np.empty((len(targ_years), len(models),len(scenarios)))
    data[:] = np.nan
    elements = [[None]*len(scenarios)]*len(models)  
        
    zostoga_am = xr.DataArray(data, coords=[targ_years, models, scenarios], dims=['year', 'model', 'scen'])
    variant_label = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    grid_label = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    incorporated = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    
    for model in models: #loop over models
    
        for scenario in scenarios: #loop over scenarios
            
            incorporate = True
                
            ssp_filenames = fnmatch.filter(os.listdir(zostoga_path+model),"*"+scenario+"*nc") #find SSP filenames for current scenario
            
            if not ssp_filenames: #if none available
                incorporate = False
                
            for ssp_filename in ssp_filenames: #loop over SSP filenames
                
                ssp_ds = xr.open_dataset(os.path.join(zostoga_path,model,ssp_filename),decode_times=True) #open SSP file
                
                #find historical with that matches variant and grid of current SSP
                hist_filename = fnmatch.filter(os.listdir(zostoga_path+model), '*historical*'+ssp_ds.parent_variant_label+'*'+ssp_ds.grid_label+'*.nc')
                
                if hist_filename: #as soon as matching historical run is found, use this SSP
                    incorporate = True
                    break
                else: #go to next SSP variant
                    incorporate = False
    
            incorporated.loc[model,scenario] = incorporate
    
            if(incorporate): #if matching historical and SSP combination found:
                hist_ds = xr.open_dataset(os.path.join(zostoga_path,model,hist_filename[0]),decode_times=True)
                
                #convert to annual means
                hist_ds_am = hist_ds.groupby('time.year').mean('time')
                ssp_ds_am = ssp_ds.groupby('time.year').mean('time')
                
                #suture models whose reference state of each experiment is the first year of the experiment
                if model=='MRI-ESM2-0':
                    diffA = (ssp_ds_am['zostoga'].isel(year=0).values-hist_ds_am['zostoga'].isel(year=-1).values) 
                    diffB = (ssp_ds_am['zostoga'].isel(year=1).values-ssp_ds_am['zostoga'].isel(year=0).values)
                    ssp_ds_am['zostoga'] = ssp_ds_am['zostoga'] - (diffA-diffB)

                #append ssp and historical, fill missing values if any
                hist_ssp_ds_am = ssp_ds_am.combine_first(hist_ds_am).interp(year=targ_years,kwargs={'fill_value': 'extrapolate'}) #.sel(year=targ_years) #combine historical and SSP and select target years
                
                #fill data arrays
                zostoga_am.loc[:,model,scenario] = np.squeeze(hist_ssp_ds_am.zostoga) #squeeze (zostoga is a 3D array for some models)
                variant_label.loc[model,scenario] = ssp_ds.variant_label
                grid_label.loc[model,scenario] = ssp_ds.grid_label
                
    #store data arrays in dataset            
    zostoga_am_ds = zostoga_am.to_dataset(name = 'zostoga')
    zostoga_am_ds['variant_label'] = variant_label
    zostoga_am_ds['grid_label'] = grid_label
    zostoga_am_ds['incorporated'] = incorporated #
    #alternatively drop models where incorporated is false for all scenarios:
    #zostoga_am_ds = zostoga_am_ds.where(zostoga_am_ds.isel(year=0).incorporated.sum(dim='scen')>0,drop=True) 
    #zostoga_am_ds = zostoga_am_ds.drop('incorporated')
    
    
    return zostoga_am_ds

#test function:
'''
targ_years = np.arange(1850,2101)
base_years = np.arange(1986,2006)
scenarios = ['ssp126','ssp245','ssp585']
zostoga_path = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted/zostoga/'

my_zostoga_ds = includeCMIP6zostoga(zostoga_path,scenarios,targ_years)

zostoga_ds = my_zostoga_ds.where(my_zostoga_ds.incorporated==True,drop=True).dropna(dim='model')
zostoga_ds = reference_to_baseyears(my_zostoga_ds,base_years,'zostoga')
zostoga_ds['incorporated'] = zostoga_ds.incorporated.astype(int)

zostoga_ds.year.attrs['standard_name'] = 'year'
zostoga_ds.year.attrs['comment'] = 'annual means from Omon'
zostoga_ds.variant_label.attrs['standard_name'] = 'variant_label'
zostoga_ds.grid_label.attrs['standard_name'] = 'grid_label'
zostoga_ds.incorporated.attrs['standard_name'] = 'flag whether model is dedrifted and incorporated in the ensemble'
zostoga_ds.model.attrs['standard_name'] = 'global_climate_model_name'
zostoga_ds.scen.attrs['standard_name'] = 'shared_socioeconomic_pathway'
zostoga_ds.zostoga.attrs['standard_name'] = 'global_average_thermosteric_sea_level_change'
zostoga_ds.zostoga.attrs['long_name'] = 'Global Average Thermosteric Sea Level Change'
zostoga_ds.zostoga.attrs['units'] = 'm'
zostoga_ds.zostoga.attrs['comment'] = 'Dedrifted with parallel part of quadratic fit of full PiControl, reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])
zostoga_ds.attrs['comment'] = 'unfiltered global mean thermosteric sea-level change ensemble for CMIP6 GCMs'
zostoga_ds.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
zostoga_ds.attrs['download_date'] = '21-05-2020'
zostoga_ds.attrs['contact'] = 'user@adress.com'


zostoga_ds.to_netcdf(os.path.join('/Volumes/Naamloos/PhD_Data/CMIP6/zostoga_tas_ens','zostoga_CMIP6_n'+str(len(zostoga_ds.model))+'_unfiltered_'+str(base_years[0])+'_'+str(base_years[-1])+'ref_qdedr_'+str(targ_years[0])+'_'+str(targ_years[-1])+'_am.nc'),
                                                     mode='w',
                                                     encoding={'zostoga':{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
'''