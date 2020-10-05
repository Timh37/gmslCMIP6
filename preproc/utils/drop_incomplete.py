
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function finds the overlap of models for which GTE and GSAT
are available for all SSPs, and drops models for which data is missing.

@author: thermans
"""
import numpy as np

def drop_incomplete(zostoga_ds,tas_ds):
    #only include models for GTE and gSAT that have incorporated=True for all ssps
    zostoga_ds_incorporated = zostoga_ds.where(zostoga_ds.incorporated==True,drop=True).dropna(dim='model')
    tas_ds_incorporated = tas_ds.where(tas_ds.incorporated==True,drop=True).dropna(dim='model')
    
    #overlapping models GTE and GSAT
    complete_models = np.intersect1d(zostoga_ds_incorporated.model.values,tas_ds_incorporated.model.values)
    
    #select models with complete GTE and GSAT data and drop the incorporate variable
    zostoga_ds_complete = zostoga_ds_incorporated.sel(model=complete_models).drop_vars('incorporated')
    tas_ds_complete = tas_ds_incorporated.sel(model=complete_models).drop_vars('incorporated')
    
    return zostoga_ds_complete, tas_ds_complete