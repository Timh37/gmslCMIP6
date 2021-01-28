'''
This script linearly dedrifts CMIP6 model 'zostoga' data using piControl. The drift
is estimated by fitting a quadratic polynomial to the full piControl run and subtracting
the overlapping part of the polynomial from the historical and scenario runs. 
The output is saved to NetCDFs, separately for experiment, variant and grid.

Input parameters:
    in_dir          input directory with time-merged cmip6 variable data (any frequency), organized by model directory
    out_dir         output directory to store dedrifted cmip6 variable data 

Output:
    NetCDF files with dedrifted timeseries
    
Created on Thu Apr 30 14:56:33 2020

@author: thermans
'''

import xarray as xr
import numpy as np
import os
import fnmatch

def fit_quadratic(data, time): #fit quadratic polynomial
    pfit = np.polyfit(time, data,2) 
    return np.array([pfit[0],pfit[1],pfit[2]]) #return coefficients, highest power first!

def evaluate_drift(pfit, time): #evaluate polynomial
    drift = np.polyval(pfit,time)
    return drift

def make_dedrifted_ds(ds,dedrifted_data): #generate dedrifted dataset from raw dataset and add drift coefficients
    dedrifted_data.attrs = ds.zostoga.attrs #copy attributes from original data
    dedrifted_ds = ds.copy(deep=True) #copy original dataset
    
    dedrifted_ds['zostoga'] = dedrifted_data #replace raw zostoga with dedrifted zostoga
    dedrifted_ds.zostoga.attrs['dedrifted'] = 'Dedrifted using overlapping quadratic polynomial fitted to full piControl'
    return dedrifted_ds

def cmip6_dedrift_zostoga_quadratic(in_dir,out_dir):
    for model in (model for model in os.listdir(in_dir) if model!=".DS_Store"):
        print('Current model: '+model)
        model_dir = os.path.join(in_dir,model)
        
        pic_fns = fnmatch.filter(os.listdir(model_dir), "*piControl*nc") #search for piControl runs for current model
    
        if not pic_fns: #skip model if no piControl available
            print('Model has no piControl')
            continue
        
        #loop over piControl runs for current model
        for pic_fn in pic_fns: 
            
            #open piControl run (decode false, xarray can't handle piC calendar years)
            pic = xr.open_dataset(os.path.join(model_dir,pic_fn),decode_times=False) 
            
            print('Current piControl variant: '+pic.variant_label+', grid: '+pic.grid_label)
             
            #search for historical filenames with this grid
            hist_fns = fnmatch.filter(os.listdir(model_dir), "*historical*"+pic.grid_label+"*nc")
            ssp_fns = fnmatch.filter(os.listdir(model_dir), "*ssp*"+pic.grid_label+"*nc")
            
            if not hist_fns: #if no historical file, we cannot dedrift, so move on to next piControl run
                print('Current piControl has no historical')
                continue #go to next piControl
                
            #apply quadratic polynomial to full piControl
            control_pfit = xr.apply_ufunc(
                    fit_quadratic, pic.zostoga , pic.time,
                    input_core_dims=[["time"], ["time"]], #core dimension time, loop over the others
                    output_core_dims=[["coefs"]], #outputs 2nd degree, 1st degree and intercept
                    vectorize=True, 
                    dask='allowed', #allow calculating in chunks (dask='parallelized' doesn't work)
                    output_dtypes=[float],
                    output_sizes={"coefs": 3}, #output must be numpy array
                    )
            
            print('Fitted quadratic polynomial to piControl')
            
            # dedrift historical        
            for hist_fn in hist_fns: #loop over historical filenames with the same grid as piControl
                hist = xr.open_dataset(os.path.join(model_dir,hist_fn),decode_times=False)
    
                overlaps = 0 #identifies if piControl time and historical time overlap
                
                #if historical is  branched from current piControl
                if (((hist.parent_experiment_id=='piControl') or (hist.parent_experiment_id=='p i C o n t r o l')) and (hist.parent_variant_label == pic.variant_label)) : #if historical is  branched from piControl
        
                    #determine parallel parts using branch time in parent
                    if model.startswith('BCC') or model.startswith('CAMS'): #CAMS and BCC give branch time in parent as year instead of day of branching
                        #Convert branch year to branch days since reference time (not completely correct, since not all years have numdays of year 1)
                        hist.attrs['branch_time_in_parent'] = (hist.branch_time_in_parent - int(pic.time.units[11:15]))*np.sum(np.diff(pic.time)[0:12])
                        hist.attrs['branch_time_comment'] = 'branch_time_in_parent has been converted to days since piControl reference time'
                    elif model=='EC-Earth3':
                        hist.attrs['branch_time_in_parent'] = 149749 #from Thomas Reerink @ KNMI
                    
                    #CASE 1: historical time referenced to 1850 and starts with 15-16 days since 1850 (month 1)
                    if (hist.time.units[11:15] == '1850') & (hist.time[0]<17): 
                        #historical parallel to piControl at hist_zostoga.time+branch time
                        time_hist_in_pic = hist.time + hist.branch_time_in_parent 
                    
                    #CASE 2: historical time referenced to earlier than 1850, so historical time starts at larger than ~15-16 days since reference time
                    elif (int(hist.time.units[11:15]) < 1850) & (hist.time[0]>17):
                        #historical parallel to piControl at hist_zostoga.time+branch time, minus the extra days introduced because of the different historical reference time
                        time_hist_in_pic = (hist.time - (hist.time[0]-15.5)) + hist.branch_time_in_parent 
                        
                    else:
                        raise Exception('Error: historical branch_time_in_parent and reference time not recognized.')   
                
                    idx_hinc = np.where((pic.time>time_hist_in_pic[0]-15) & (pic.time<time_hist_in_pic[-1]+15)) #find where historical time is approximately contained in piControl
                    overlap_hist_pic, idx_hisc, idx_cish = np.intersect1d(pic.time,time_hist_in_pic,return_indices=True) #find exact overlap between piControl and historical times
                    
                    if np.array_equal(time_hist_in_pic,overlap_hist_pic): #if control and historical times overlap exactly      
                        overlaps = 1
                    
                    elif len(pic.time[idx_hinc])==len(time_hist_in_pic): #if historical time is contained in control time, but no precise overlap
                        idx_hisc = np.squeeze(idx_hinc)
                        overlaps = 1
                        
                    else:
                        print('piControl and historical do not overlap, cannot dedrift historical')
                            
                    if overlaps == 1: #if piControl and historical overlap, we can dedrift
                        #evaluate polynomial of piControl for historical
                        hist_drift = xr.apply_ufunc(
                            evaluate_drift,control_pfit,time_hist_in_pic, 
                            input_core_dims=[["coefs"], ["time"]], 
                            output_core_dims=[["time"]], 
                            vectorize=True,
                            dask='allowed',
                            output_dtypes=[hist.zostoga.dtype],
                            )
                        hist_drift['time'] = hist.time
                        
                        hist_drift = hist_drift - hist_drift.isel(time=0) #we're interested in subtracting the slope, not absolute values
                        hist_dedrifted = hist.zostoga - hist_drift #subtract drift
                    
                        hist_dedrifted_ds = make_dedrifted_ds(hist,hist_dedrifted)
                        
                        if not os.path.exists(out_dir+model): #make path
                            os.mkdir(out_dir+model)
                    
                        #save to netcdf, use overwrite and compress 'zostoga' variable
                        hist_dedrifted_ds.to_netcdf(os.path.join(out_dir,model,hist_fn[0:-3]+'_dedrifted.nc'),
                                                 mode='w',
                                                 encoding={'zostoga':{'zlib': True,'complevel': 1}}
                                                 ) #(over)write a new NetCDF file
                        
                        print('Dedrifted '+hist.experiment_id+' ('+hist.variant_label+')')  
                        
                        hist_dedrifted_ds.close()
                        
                else: #current historical not branched from current piControl    
                    print('Current ' +hist.experiment_id+' ('+hist.variant_label+') not branched from current piControl')
                    continue #move on to next historical (we cannot dedrift SSPs either in this case)
                
                # dedrift SSPs
                for ssp_fn in ssp_fns:
                    ssp = xr.open_dataset(os.path.join(model_dir,ssp_fn),decode_times=False)
                
                    overlaps = 0
                    
                    #if current ssp is branched from current historical
                    if (((ssp.parent_experiment_id=='historical') or (ssp.parent_experiment_id=='h i s t o r i c a l')) and (ssp.parent_variant_label == hist.variant_label)): #if ssp is branched from historical
                        
                        #determine parallel parts using branch time in parent
                        if model.startswith('BCC') or model.startswith('CAMS'): #CAMS and BCC give branch time in parent as year of branching
                            #Convert branch year to branch days since reference time
                            ssp.attrs['branch_time_in_parent'] = (ssp.branch_time_in_parent - int(hist.time.units[11:15]))*np.sum(np.diff(pic.time)[0:12])
                            ssp.attrs['branch_time_comment'] = 'branch_time_in_parent has been converted to days since historical reference time'
                        elif model=='EC-Earth3':
                            ssp.attrs['branch_time_in_parent'] = 60265 #from Thomas Reerink @ KNMI
                        
                        #CASE 1: SSP and historical are both referenced w.r.t. 1850
                        if (ssp.time.units[11:15] == '1850') & (ssp.time.units == hist.time.units):
                            #branch time of ssp in historical already included in ssp_time if referenced w.r.t. 1850, so parallel to piControl at ssp_zostoga.time + historical branch time
                            time_ssp_in_pic = ssp.time + hist.branch_time_in_parent
                        
                        #CASE 2: SSP is referenced w.r.t. 2015, so ssp time starts with 15-16 days since 2015    
                        elif (ssp.time.units[11:15] == '2015') & (ssp.time[0]<17): 
                            #parallel to piControl at ssp.time + historical branch time in piControl + ssp branch time in historical
                            time_ssp_in_pic = ssp.time + hist.branch_time_in_parent + ssp.branch_time_in_parent
                            
                        #CASE 3: SSP time and historical time are both referenced to earlier than 1850
                        elif (int(ssp.time.units[11:15])<1850) & (ssp.time.units == hist.time.units): #(NorESM for instance)
                            #parallel at ssp_zostoga.time + historical branch time, minus the extra days introduced because of the different historical and ssp reference time
                            time_ssp_in_pic = ssp.time + hist.branch_time_in_parent - (hist.time[0]-15.5)
                        
                        else:
                            raise Exception('Error: ssp branch_time_in_parent and reference time not recognized.')  
                
                        idx_sspinc = np.where((pic.time>=time_ssp_in_pic[0]-15) & (pic.time<=time_ssp_in_pic[-1]+15)) #find where ssp time is approximately contained in piControl
                        overlap_ssp_pic, idx_sspisc, idx_cisssp = np.intersect1d(pic.time,time_ssp_in_pic,return_indices=True) #find exact overlap between piControl and ssp times
                 
                        if np.array_equal(time_ssp_in_pic,overlap_ssp_pic): #if control and ssp times overlap exactly  
                            overlaps = 1
                
                        #CASE 2: if ssp time is contained in control time but for some reason shifted by less than half a month
                        elif len(pic.time[idx_sspinc])==len(time_ssp_in_pic): #if ssp time is contained in control time, but no precise overlap
                            idx_sspisc = np.squeeze(idx_sspinc)
                            overlaps = 1
                        
                        else:
                            print('piControl and '+ssp.experiment_id+' do not overlap, cannot dedrift '+ssp.experiment_id)
                            
                        if overlaps == 1: #if piControl and ssp overlap, we can dedrift
                            ssp_drift = xr.apply_ufunc(  
                                evaluate_drift,control_pfit,time_ssp_in_pic,
                                input_core_dims=[["coefs"], ["time"]], 
                                output_core_dims=[["time"]], 
                                vectorize=True,
                                dask='allowed',
                                output_dtypes=[hist.zostoga.dtype],
                                )
                            ssp_drift['time'] = ssp.time
                            #subtract parallel piControl drift from raw ssp zostoga (note we subtract the fit, but not the initial value at hist_zostoga.time 0)
                            ssp_dedrifted = ssp.zostoga - ssp_drift + hist_drift.isel(time=0)
                        
                            ssp_dedrifted_ds = make_dedrifted_ds(ssp,ssp_dedrifted)
            
                            if not os.path.exists(os.path.join(out_dir,model)): #make path
                                    os.mkdir(os.path.join(out_dir,model))
                            
                            #save to netcdf, use overwrite and compress 'zostoga' variable
                            ssp_dedrifted_ds.to_netcdf(os.path.join(out_dir+model,ssp_fn[0:-3]+'_dedrifted.nc'),
                                                     mode='w',
                                                     encoding={'zostoga':{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
                            
                            print('Dedrifted '+ssp.experiment_id+' ('+ssp.variant_label+')')  
                            
                            ssp_dedrifted_ds.close()
                    else:
                        print('Current ' +ssp.experiment_id+' ('+ssp.variant_label+') not branched from current historical')
                        continue
                            
                    ssp.close()
                    
                hist.close()
            pic.close()
#for testing
'''
if __name__ == '__main__':
    in_dir = 'insert_path_here' #set raw, time-merged zostoga directory
    out_dir = 'insert_path_here' #set output directory

    cmip6_dedrift_zostoga_quadratic(in_dir,out_dir)
'''
