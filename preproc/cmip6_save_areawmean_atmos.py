'''
Saves area weighted mean atmosphere variable using its grid cell's areas.


Input parameters:
    variable    atmosphere variable of interest
    in_dir      path to NetCDF files with ocean variable, ordered by model
    out_dir     path to store the corrected files
    area_dir    path to areacella (atmosphere grid cell area)

Output:
    awmean   area weighted mean of input variable
    
Created by: Tim Hermans, 29-07-20
'''
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import fnmatch

def cmip6_save_areawmean_atmos(in_dir,out_dir,area_dir,variable):
    for i,model in enumerate((model for model in os.listdir(in_dir) if model!=".DS_Store")): #loop over models
        model_dir = os.path.join(in_dir,model) #model input directory
        print(model)
        
        if model not in os.listdir(area_dir): #if no directory for model in path to areaweights
            hasAreaweights = False
            print(model+' has no areaweights, moving on') 
            continue
        else:
            hasAreaweights = True
            
        for file in sorted(os.listdir(model_dir)): #loop over all experiments for current model
    
            hasAreaweights = True
            
            model_area_dir = os.path.join(area_dir,model) #path to areaweights for model
            areafile=[]
     
            var_ds = xr.open_dataset(os.path.join(model_dir,file),decode_times=False) #load variable file
            coordinates = list(k for k in var_ds[variable].dims if 'time' not in k) #find lon/lat coordinate names
            
            areafile = fnmatch.filter(os.listdir(model_area_dir),"*"+var_ds.grid_label+"*.nc") #look for matching areacella (same grid label)
            
            if not areafile: #if none available
                hasAreaweights=False
                print(model+', '+var_ds.grid_label+' has no areaweights, moving on')
                continue
                
            else:
                areaweights = xr.open_dataset(os.path.join(model_area_dir,areafile[0])) #load first available areaweights with matching grid
                
            if ((np.sum(areaweights.areacella)<5e14) or (np.sum(areaweights.areacella)>5.2e14)): #if no valid total atmosphere surface area, skip file
                print('invalid total surface area for: ' + file)
                continue
            
            #calculate areaweighted mean
            with xr.set_options(keep_attrs=True):
                awmean = ((var_ds[variable] * np.repeat(areaweights.areacella.values[np.newaxis,:,:],len(var_ds.time),axis=0)).sum(dim=coordinates,skipna=True)) / np.sum(areaweights.areacella)
            var_ds[variable] = awmean #replace 'tas' with awmean 'tas'
            
            #save area-weighted mean variable file to netcdf
            if not os.path.exists(os.path.join(out_dir,model)): #make path
                os.mkdir(os.path.join(out_dir,model))
                
            var_ds.to_netcdf(os.path.join(out_dir,model,file[0:-3]+'_awmean.nc'),
                                                         mode='w',
                                                         encoding={variable:{'zlib': True,'complevel': 1,'dtype':'float32'}}
                                                         ) #(over)write a new NetCDF file
            
'''
if __name__ == '__main__':
    variable = 'tas'
    in_dir = 'insert_path_here/tas/' #+'/'+model+'/'
    out_dir = 'insert_path_here/gsat/'
    area_dir = 'insert_path_here/areacella/' #areaweights
    cmip6_save_areawmean_atmos(in_dir,out_dir,area_dir,variable)
'''