#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function incorporates area-weighted mean 'tas' data of CMIP6 models in one Xarray Dataset.
It looks for the first variant and grid member of each model for which both
scenario and historical run are available. The monthly means are converted to 
annual means and the historical and future period are combined along the time axis.

input: 
targ_years          The years to select from historical to ssp (no interpolation implemented yet)
scenarios           (List of) desired SSPs
tas_path            Path to area-weighted mean tas stored per variable per model

output:
tas_ds_am           Dataset with area-weighted mean tas as function of year, model and scenario
    -tas                 -area-weighted tas data of model for this scenario (nan if model not incorporated)
    -variant_label       -Variant label of model for this scenario (none if model not incorporated)
    -grid_label          -Grid label of model for this scenario (none if model not incorporated)
    -incorporated        -If model data is available for this scenario (True or False)


Created on Fri May 15 15:08:40 2020

@author: thermans
"""
import xarray as xr
import numpy as np
import os
import fnmatch

def includeCMIP6tas(tas_path,scenarios,targ_years):

    models = [model for model in os.listdir(tas_path) if not '.' in model] #get list of model directories
    
    #initialize empty data arrays
    data = np.empty((len(targ_years), len(models),len(scenarios)))
    data[:] = np.nan
    elements = [[None]*len(scenarios)]*len(models)  
        
    tas_am = xr.DataArray(data, coords=[targ_years, models, scenarios], dims=['year', 'model', 'scen'])
    variant_label = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    grid_label = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    incorporated = xr.DataArray(elements, coords=[models,scenarios], dims=['model', 'scen'])
    
    for model in models: #loop over models
    
        for scenario in scenarios: #loop over scenarios
            
            incorporate = True
                
            ssp_filenames = fnmatch.filter(os.listdir(tas_path+model),"*"+scenario+"*nc") #find SSP filenames for current scenario
            
            if not ssp_filenames: #if none available
                incorporate = False
                
            for ssp_filename in ssp_filenames: #loop over SSP filenames
                
                ssp_ds = xr.open_dataset(os.path.join(tas_path,model,ssp_filename),decode_times=True) #open SSP file
                
                #find historical with that matches variant and grid of current SSP
                hist_filename = fnmatch.filter(os.listdir(tas_path+model), '*historical*'+ssp_ds.parent_variant_label+'*.nc')
                
                if hist_filename: #as soon as matching historical run is found, use this SSP
                    incorporate = True
                    break
                else: #go to next SSP file
                    incorporate = False
    
            incorporated.loc[model,scenario] = incorporate
    
            if(incorporate): #if matching historical and SSP combination found:
                hist_ds = xr.open_dataset(os.path.join(tas_path,model,hist_filename[0]),decode_times=True) #open historical file
                
                #convert to annual means
                hist_ds_am = hist_ds.groupby('time.year').mean('time')
                ssp_ds_am = ssp_ds.groupby('time.year').mean('time')

                #append ssp and historical, fill missing values if any
                hist_ssp_ds_am = ssp_ds_am.combine_first(hist_ds_am).interp(year=targ_years,kwargs={'fill_value': 'extrapolate'}) #.sel(year=targ_years) #combine historical and SSP and select target years
                
                #fill data arrays
                tas_am.loc[:,model,scenario] = np.squeeze(hist_ssp_ds_am.tas) #squeeze (the time series is a 3D array for some models)
                variant_label.loc[model,scenario] = ssp_ds.variant_label 
                grid_label.loc[model,scenario] = ssp_ds.grid_label
                
    #store arrays in dataset           
    tas_am_ds = tas_am.to_dataset(name = 'tas')
    tas_am_ds['variant_label'] = variant_label
    tas_am_ds['grid_label'] = grid_label
    tas_am_ds['incorporated'] = incorporated
    #alternatively drop models where incorporated is false for all scenarios:
    #tas_am_ds = tas_am_ds.where(tas_am_ds.isel(year=0).incorporated.sum(dim='scen')>0,drop=True) 
    #tas_am_ds = tas_am_ds.drop('incorporated')
    return tas_am_ds

#test:
'''
targ_years = np.arange(1850,2101)
base_years = np.arange(1986,2006)
scenarios = ['ssp126','ssp245','ssp585']
tas_path = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/tas/'

my_tas_ds = includeCMIP6tas(tas_path,scenarios,targ_years)
tas_ds = reference_to_baseyears(my_tas_ds,base_years,'tas')
tas_ds['incorporated'] = tas_ds.incorporated.astype(int)

tas_ds.year.attrs['standard_name'] = 'year'
tas_ds.variant_label.attrs['standard_name'] = 'variant_label'
tas_ds.grid_label.attrs['standard_name'] = 'grid_label'
tas_ds.incorporated.attrs['standard_name'] = 'flag whether model scenario incorporated in the ensemble'
tas_ds.model.attrs['standard_name'] = 'global_climate_model_name'
tas_ds.scen.attrs['standard_name'] = 'shared_socioeconomic_pathway'
tas_ds.tas.attrs['standard_name'] = 'global_average_surface_air_temperature'
tas_ds.tas.attrs['long_name'] = 'Global average surface air temperature'
tas_ds.tas.attrs['units'] = 'K'
tas_ds.tas.attrs['comment'] = 'The reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])
tas_ds.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
tas_ds.attrs['download_date'] = '20-05-2020'
tas_ds.attrs['contact'] = 'user@adress.com'

tas_ds.to_netcdf(os.path.join('/Volumes/Naamloos/PhD_Data/CMIP6/zostoga_tas_ens','tas_CMIP6_n'+str(len(tas_ds.model))+'_unfiltered_'+str(base_years[0])+'_'+str(base_years[-1])+'ref_'+str(targ_years[0])+'_'+str(targ_years[-1])+'_am.nc'),
                                                     mode='w',
                                                     encoding={'tas':{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
'''