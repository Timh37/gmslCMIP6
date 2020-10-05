
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function generates a NetCDF file with complete ensemble of dedrifted GTE
and GSAT data for CMIP6 GCMs, referenced to base years. 
It also saves the ensemble mean and standard deviation.

Input:
target_years:   timeseries years wished to incorporate
base_years:     years to reference timeseries to
scenarios:      scenarios wished to incorporate
zostoga_path:   path to zostoga data stored per variable per model
tas_path:       path to GSAT data stored per variable per model
output_path:    path to store the ensemble netcdfs
Output:
NetCDF files of ensemble GTE and GSAT time series

@author: thermans
"""
import numpy as np
from utils.includeCMIP6zostoga import includeCMIP6zostoga
from utils.includeCMIP6tas import includeCMIP6tas
from utils.drop_incomplete import drop_incomplete
from utils.reference_var_to_baseyears import reference_var_to_baseyears
import datetime
import os

targ_years = np.arange(1850,2101)
base_years = np.arange(1986,2006)
scenarios = ['ssp126','ssp245','ssp585']
zostoga_path = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_quadratic/zostoga/'
tas_path = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/tas/'
output_path = '/Volumes/Naamloos/PhD_Data/CMIP6/zostoga_tas_ens/'

#pool zostoga and tas data
zostoga_am_ds = includeCMIP6zostoga(zostoga_path,scenarios,targ_years)
tas_am_ds = includeCMIP6tas(tas_path,scenarios,targ_years)

#GISS-E2-1-G is discarded, since it simulates a historical sl fall of 3m
zostoga_am_ds = zostoga_am_ds.drop_sel(model=['GISS-E2-1-G'])

#drop incomplete models (we require model to have all SSPs for both gte and gsat)
zostoga_ds_complete, tas_ds_complete = drop_incomplete(zostoga_am_ds,tas_am_ds)

#reference to base years
zostoga_ds_complete = reference_var_to_baseyears(zostoga_ds_complete,base_years,'zostoga')
tas_ds_complete = reference_var_to_baseyears(tas_ds_complete,base_years,'tas')

#add ensemble mean and standard deviation
tas_ds_complete['tas_ens_mean'] = tas_ds_complete.tas.mean(dim='model')
tas_ds_complete.tas_ens_mean.attrs['standard_name'] = 'ens_mean_global_average_surface_air_temperature'
tas_ds_complete.tas_ens_mean.attrs['long_name'] = 'Ensemble mean global average surface air temperature'
tas_ds_complete.tas_ens_mean.attrs['units'] = 'K'
tas_ds_complete.tas_ens_mean.attrs['comment'] = 'The reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])

tas_ds_complete['tas_ens_std'] = tas_ds_complete.tas.std(dim='model')
tas_ds_complete.tas_ens_std.attrs['standard_name'] = 'ens_std_global_average_surface_air_temperature'
tas_ds_complete.tas_ens_std.attrs['long_name'] = 'Ensemble standard deviation global average surface air temperature'
tas_ds_complete.tas_ens_std.attrs['units'] = 'K'
tas_ds_complete.tas_ens_std.attrs['comment'] = 'The reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])

zostoga_ds_complete['zostoga_ens_mean'] = zostoga_ds_complete.zostoga.mean(dim='model')
zostoga_ds_complete.zostoga_ens_mean.attrs['standard_name'] = 'ens_mean_global_average_thermosteric_sea_level_change'
zostoga_ds_complete.zostoga_ens_mean.attrs['long_name'] = 'Ensemble mean global average thermosteric sea level change'
zostoga_ds_complete.zostoga_ens_mean.attrs['units'] = 'm'
zostoga_ds_complete.zostoga_ens_mean.attrs['comment'] = 'Dedrifted with parallel part of quadratic fit of full PiControl, reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])

zostoga_ds_complete['zostoga_ens_std'] = zostoga_ds_complete.zostoga.std(dim='model')
zostoga_ds_complete.zostoga_ens_std.attrs['standard_name'] = 'ens_std_global_average_thermosteric_sea_level_change'
zostoga_ds_complete.zostoga_ens_std.attrs['long_name'] = 'Ensemble standard deviation global average thermosteric sea level change'
zostoga_ds_complete.zostoga_ens_std.attrs['units'] = 'm'
zostoga_ds_complete.zostoga_ens_std.attrs['comment'] = 'Dedrifted with parallel part of quadratic fit of full PiControl, reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])

#add variable and global attributes
tas_ds_complete.year.attrs['standard_name'] = 'year'
tas_ds_complete.year.attrs['comment'] = 'annual means from Amon'
tas_ds_complete.variant_label.attrs['standard_name'] = 'variant_label'
tas_ds_complete.grid_label.attrs['standard_name'] = 'grid_label'
tas_ds_complete.model.attrs['standard_name'] = 'global_climate_model_name'
tas_ds_complete.scen.attrs['standard_name'] = 'shared_socioeconomic_pathway'
tas_ds_complete.tas.attrs['standard_name'] = 'global_average_surface_air_temperature'
tas_ds_complete.tas.attrs['long_name'] = 'Global average surface air temperature'
tas_ds_complete.tas.attrs['units'] = 'K'
tas_ds_complete.tas.attrs['comment'] = 'The reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])
tas_ds_complete.attrs['comment'] = 'global mean surface air temperature ensemble for CMIP6 GCMs for which tas and dedrifted zostoga is available for both the historical and '+str(scenarios)+', derived from Chris Jones tas data'
tas_ds_complete.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
tas_ds_complete.attrs['download_date'] = '20-05-2020'
tas_ds_complete.attrs['contact'] = 'tim.hermans@nioz.nl'

zostoga_ds_complete.year.attrs['standard_name'] = 'year'
zostoga_ds_complete.year.attrs['comment'] = 'annual means from Omon'
zostoga_ds_complete.variant_label.attrs['standard_name'] = 'variant_label'
zostoga_ds_complete.grid_label.attrs['standard_name'] = 'grid_label'
zostoga_ds_complete.model.attrs['standard_name'] = 'global_climate_model_name'
zostoga_ds_complete.scen.attrs['standard_name'] = 'shared_socioeconomic_pathway'
zostoga_ds_complete.zostoga.attrs['standard_name'] = 'global_average_thermosteric_sea_level_change'
zostoga_ds_complete.zostoga.attrs['long_name'] = 'Global Average Thermosteric Sea Level Change'
zostoga_ds_complete.zostoga.attrs['units'] = 'm'
zostoga_ds_complete.zostoga.attrs['comment'] = 'Dedrifted with parallel part of quadratic fit of full PiControl, reference state is set to the mean of '+str(base_years[0])+'-'+str(base_years[-1])
zostoga_ds_complete.attrs['comment'] = 'global mean thermosteric sea-level change ensemble for CMIP6 GCMs for which tas and dedrifted zostoga is available for both the historical and '+str(scenarios)
zostoga_ds_complete.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
zostoga_ds_complete.attrs['download_date'] = '21-05-2020'
zostoga_ds_complete.attrs['contact'] = 'tim.hermans@nioz.nl'

#store to netcdf
zostoga_ds_complete.to_netcdf(os.path.join(output_path,'zostoga_CMIP6_n'+str(len(zostoga_ds_complete.model))+'_'+str(base_years[0])+'_'+str(base_years[-1])+'ref_qdedr_'+str(targ_years[0])+'_'+str(targ_years[-1])+'_am.nc'),
                                                     mode='w',
                                                     encoding={'zostoga':{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
tas_ds_complete.to_netcdf(os.path.join(output_path,'tas_CMIP6_n'+str(len(zostoga_ds_complete.model))+'_'+str(base_years[0])+'_'+str(base_years[-1])+'ref_awm_'+str(targ_years[0])+'_'+str(targ_years[-1])+'_am.nc'),
                                                     mode='w',
                                                     encoding={'tas':{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
