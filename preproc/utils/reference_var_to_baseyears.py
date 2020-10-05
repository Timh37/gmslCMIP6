
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function references a variable from a dataset to the queried reference
year or reference period.

Input:
    am_ds       annual mean dataset with time dimension 'year'
    baseyears   reference period
    var         variable to reference
Output:
    am_ds       annual mean dataset with variable referenced to baseyars

@author: thermans
"""
import numpy as np

def reference_var_to_baseyears(am_ds,baseyears,var):
        
    if np.isscalar(baseyears):
        am_ds[var] = am_ds[var] - am_ds[var].sel(year=baseyears)
    else:
        am_ds[var] = am_ds[var] - am_ds[var].sel(year=baseyears).mean(dim='year')
        
    return am_ds