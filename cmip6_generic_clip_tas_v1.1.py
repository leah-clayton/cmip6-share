#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:22:20 2024

@author: leahclayton
"""

import xarray as xr
import os

# set base input location for data, ending in /AMES/NEX/GDDP-CMIP6 (no / at end)
input_loc = '/file/path/AMES/NEX/GDDP-CMIP6'

# set output location (folder doesn't have to exist)
output_loc = '/file/path/save_folder'

# define data time period/ssp
ssps = ['ssp245']

# define variable
variable = 'tas'

# define time period
start_year = 2015
final_year = 2100

# define spatial extent
lon_min, lon_max = 235.370, 257.875
lat_min, lat_max = 31.125, 49.125

#%% don't edit below this line

# final range year
final_range = final_year + 1

#%%
# all models
model = "TaiESM1"

variant = 'r1i1p1f1' # TaiESM1

# two letter component of file name by model
code = 'gn' # TaiESM1

#%%
def create_dir_if_not_exists(save_loc):
    if not os.path.exists(save_loc):
        os.makedirs(save_loc)
    else:
        print(f"Folder '{save_loc}' already exists.")
            
#%%
for ssp in ssps:
    folder_path = input_loc + f'/{model}/{ssp}/{variant}'
    
    for year in range(start_year, final_range):
        save_loc = output_loc + f'/{year}/{model}'
        create_dir_if_not_exists(save_loc)
    
        file_path = folder_path + f'/{variable}/{variable}_day_{model}_{ssp}_{variant}_{code}_{year}_v1.1.nc'
        save_path = save_loc + f'/cmip6_{model}_{variable}_{ssp}_{year}.nc'
        ds = xr.open_dataset(file_path)
        
        clip_ds_lon = ds.where((ds['lon'] >= lon_min) & (ds['lon'] <= lon_max), drop=True)
        clip_ds = clip_ds_lon.where((ds['lat'] >= lat_min) & (ds['lat'] <= lat_max), drop=True)
        print(clip_ds)
        
        clip_ds.to_netcdf(save_path)
        
        ds.close()
        clip_ds.close()