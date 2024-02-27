#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:15:07 2024

@author: leahclayton
"""

import xarray as xr
import os

# set base input location for data, ending in /AMES/NEX/GDDP-CMIP6 (no / at end)
input_loc = '/file/path/AMES/NEX/GDDP-CMIP6'

# set output location (folder doesn't have to exist)
output_loc = '/file/path/save_folder'

# define data time period/ssp ('historical', 'ssp245', 'ssp585')
ssp = 'ssp245'

# define variables (list any variables being used)
variables = ['pr', 'rsds', 'tas']

# define time period
start_year = 2015
final_year = 2100

# define spatial extent for clipping using NASA NEX GDDP CMIP6 convention 
# (given is for Western US, Eastern boarder of CO and west)
lon_min, lon_max = 235.370, 257.875
lat_min, lat_max = 31.125, 49.125

#%% comment out what models you aren't using and coresponding variants and codes

# all models
models = ["ACCESS-CM2",
    "ACCESS-ESM1-5",
    "CanESM5",
    "CMCC-CM2-SR5", # tas, tasmin, and tasmax retracted but REPLACED
    "CMCC-ESM2", # no ssp126 ssp370
    "CNRM-CM6-1", # no ssp126 ssp370
    "CNRM-ESM2-1",
    "EC-Earth3",
    "EC-Earth3-Veg-LR",
    "FGOALS-g3",
    "GFDL-CM4", # no ssp126 ssp370
    "GFDL-ESM4",
    "GISS-E2-1-G",
    "HadGEM3-GC31-LL", # no ssp370
    "INM-CM4-8",
    "INM-CM5-0",
    "KACE-1-0-G",
    "MIROC-ES2L",
    "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR",
    "MRI-ESM2-0",
    "NorESM2-LM",
    "NorESM2-MM",
    "TaiESM1",
    "UKESM1-0-LL"
]

variants = ['r1i1p1f1', # ACCESS-CM2
    'r1i1p1f1', # ACCESS-ESM1-5
    'r1i1p1f1', # CanESM5
    'r1i1p1f1', # CMCC-CM2-SR5
    'r1i1p1f1', # CMCC-ESM2
    'r1i1p1f2', # CNRM-CM6-1
    'r1i1p1f2', # CNRM-ESM2-1
    'r1i1p1f1', # EC-Earth3
    'r1i1p1f1', # EC-Earth3-Veg-LR
    'r3i1p1f1', # FGOALS-g3
    'r1i1p1f1', # GFDL-CM4
    'r1i1p1f1', # GFDL-ESM4
    'r1i1p1f2', # GISS-E2-1-G
    'r1i1p1f3', # HadGEM3-GC31-LL
    'r1i1p1f1', # INM-CM4-8 
    'r1i1p1f1', # INM-CM5-0
    'r1i1p1f1', # KACE-1-0-G
    'r1i1p1f2', # MIROC-ES2L
    'r1i1p1f1', # MPI-ESM1-2-HR
    'r1i1p1f1', # MPI-ESM1-2-LR
    'r1i1p1f1', # MPI-ESM2-0
    'r1i1p1f1', # NorESM2-LM
    'r1i1p1f1', # NorESM2-MM
    'r1i1p1f1', # TaiESM1
    'r1i1p1f2' #UKESM1-0-LL
    ]

# two letter component of file name by model
codes= ['gn', # ACCESS-CM2
    'gn', # ACCESS-ESM1-5
    'gn', # CanESM5
    'gn', # CMCC-CM2-SR5
    'gn', # CMCC-ESM2
    'gr', # CNRM-CM6-1
    'gr', # CNRM-ESM2-1
    'gr', # EC-Earth3
    'gr', # EC-Earth3-Veg-LR
    'gn', # FGOALS-g3
    'gr1', # GFDL-CM4
    'gr1', # GFDL-ESM4
    'gn', # GISS-E2-1-G
    'gn', # HadGEM3-GC31-LL
    'gr1', # INM-CM4-8 
    'gr1', # INM-CM5-0
    'gr', # KACE-1-0-G
    'gn', # MIROC-ES2L
    'gn', # MPI-ESM1-2-HR
    'gn', # MPI-ESM1-2-LR
    'gn', # MPI-ESM2-0
    'gn', # NorESM2-LM
    'gn', # NorESM2-MM
    'gn', # TaiESM1
    'gn' #UKESM1-0-LL
    ]
#%% don't edit below this line

# final range year
final_range = final_year + 1

#%%
def create_dir_if_not_exists(save_loc):
    if not os.path.exists(save_loc):
        os.makedirs(save_loc)
    else:
        print(f"Folder '{save_loc}' already exists.")
            
#%%
for i in range(0, 26):
    model = models[i]
    print(model)
    
    variant = variants[i]
    code = codes[i]
    folder_path = input_loc + f'/{model}/{ssp}/{variant}'
    
    for year in range(start_year, final_range):
        save_loc = output_loc + f'/{year}/{model}'
        create_dir_if_not_exists(save_loc)
        
        for variable in variables:
            file_path = folder_path + f'/{variable}/{variable}_day_{model}_{ssp}_{variant}_{code}_{year}.nc'
            save_path = save_loc + f'/cmip6_{model}_{variable}_{ssp}_{year}.nc'
            ds = xr.open_dataset(file_path)
            
            clip_ds_lon = ds.where((ds['lon'] >= lon_min) & (ds['lon'] <= lon_max), drop=True)
            clip_ds = clip_ds_lon.where((ds['lat'] >= lat_min) & (ds['lat'] <= lat_max), drop=True)
            print(clip_ds)
            
            clip_ds.to_netcdf(save_path)
            
            ds.close()
            clip_ds.close()