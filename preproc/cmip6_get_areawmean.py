'''
This script calculates the area weighted average of a CMIP6 2D field using
its grid cell's areas.

Input parameters:
    data        xarray DataArray containing the variable
    cell_area   xarray DataArray containing cell areas (areacella or areacello depending on realm)
    realm       variable in ocean or atmos realm 

Output:
    areawmean   area weighted mean of input variable
    
Created by: Tim Hermans, 29-07-20
'''
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def getAreaWeightedMean(data, cell_area, realm): #remove area weighted mean from variables
    coordinates = list(k for k in data.dims if 'time' not in k) #get coordinate dimensions
    
    if realm == 'atmos':
        total_area = cell_area.sum(skipna=True)   
    elif realm =='ocean':
        cell_area = cell_area*np.isfinite(data.isel(time=0)) #mask with land values from data
        total_area = cell_area.sum(skipna=True)   
    else:
        raise TypeError("Expected 'ocean' or 'atmos' for realm, got " + realm)

    areawmean = (data * cell_area).sum(dim=coordinates,skipna=True) / total_area
        
    return areawmean

'''EXAMPLE
path = '/Volumes/Naamloos/PhD_Data/CMIP6/downloading/07_2020/areacello/CanESM5/'

area = xr.open_dataset(path+'areacello_Ofx_CanESM5_historical_r1i1p1f1_gn.nc') #apply zos mask in case areacello is not masked on land
ssp585 = xr.open_dataset('/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/zos/CanESM5/zos_Omon_CanESM5_ssp585_r1i1p1f1_gn_201501-218012.nc')

test= getAreaWeightedMean(ssp585.zos,area.areacello,'atmos')
plt.figure()
test.plot()
'''