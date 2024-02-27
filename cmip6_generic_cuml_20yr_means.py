#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:21:36 2024

@author: leahclayton
"""

import xarray as xr
import numpy as np
import os
import rasterio
from rasterio.transform import from_bounds
import numpy as np

"""
User Inputs -- this script designed for 'app_water' input
"""
# Enter the year start and end dates as integers
start_year = 2021
time_1_end = 2040
time_2_end = 2060
time_3_end = 2080
end_year = 2100

# Set SSPs/time periods ('historical', 'ssp245', 'ssp585')
ssps = ['ssp245','ssp585']

# model options, ones with comments not used in ensemble as of 2/15/24
models = ["ACCESS-CM2",
    "ACCESS-ESM1-5",
    # "BCC-CSM2-MR", # no hurs
    # "CESM2", # no tasmax or tasmin
    # "CESM2-WACCM", # no tasmax or tasmin # no ssp126 ssp370
    "CanESM5",
    "CMCC-CM2-SR5", # tas, tasmin, and tasmax retracted but REPLACED
    "CMCC-ESM2", # no ssp126 ssp370
    "CNRM-CM6-1", # no ssp126 ssp370
    "CNRM-ESM2-1",
    "EC-Earth3",
    "EC-Earth3-Veg-LR",
    "FGOALS-g3",
    "GFDL-CM4", # no ssp126 ssp370
    # "GFDL-CM4_gr2", # no ssp126 ssp370
    "GFDL-ESM4",
    "GISS-E2-1-G",
    "HadGEM3-GC31-LL", # no ssp370
    # "HadGEM3-GC31-MM", # no ssp370 or ssp245
    # "IITM-ESM", # no tasmax or tasmin
    "INM-CM4-8",
    "INM-CM5-0",
    # "IPSL-CM6A-LR",# no huss
    "KACE-1-0-G",
    # "KIOST-ESM", # no ssp370 # hurs ssp245 missing 2058 # hurs ssp126 missing 2023
    "MIROC-ES2L",
    # "MIROC6",
    "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR",
    "MRI-ESM2-0",
    # "NESM3", # no hurs or huss # no ssp370
    "NorESM2-LM",
    "NorESM2-MM",
    "TaiESM1",
    "UKESM1-0-LL"]

"""
------ Do not edit any code below this line except input and output paths ------------
Edit input and output paths in lines 85 and 86 to match use
"""
#%% initialization
time_2_start = time_1_end + 1
time_3_start = time_2_end + 1
time_4_start = time_3_end + 1
end_range = end_year + 1

# check the existence of the output directory
def create_dir_if_not_exists(output_loc):
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)
    else:
        print(f"Folder '{output_loc}' already exists.")
        
#%% Start calculations
for ssp in ssps:
    input_path = f'/home/cmip6_{ssp}_app_water'
    output_path = f'/home/lkc33/cmip6_{ssp}_app_water_analysis'
    create_dir_if_not_exists(output_path)

    #%% Time period 1: 2021-2040    
    models_1 = []
    for model in models: 
        annual_cuml_years = []
        # open each year and calculate annual cumlative applied water
        for year in range(start_year, time_2_start):
            file_path = input_path + f'/{year}/cmip6_rddr_app_water_{ssp}_{model}_{year}.nc'
            ds_year = xr.open_dataset(file_path)
            annual_cuml = ds_year.sum(dim='time')
            annual_cuml_years.append(annual_cuml)
            ds_year.close()
        
        # combine all annual cumulative applied water netCDFs and find the mean for 20 years
        ds_years = xr.concat(annual_cuml_years, dim='year')
        mean_20yrs = ds_years.mean(dim='year')
        models_1.append(mean_20yrs)
        ds_years.close()
    
    # calculate the ensemble mean from all models and save output netCDF
    time_1_ensemble_mems = xr.concat(models_1, dim='ensembles')
    time_1_ensemble = time_1_ensemble_mems.mean(dim='ensembles')
    print(time_1_ensemble)
    
    # save as netCDF
    time_1_save = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{start_year}_{time_1_end}.nc'
    time_1_ensemble.to_netcdf(time_1_save)
    
    # save as raster
    time_1_array = time_1_ensemble['app_water']
    print(time_1_array)
    time_1_raster = np.flip(time_1_array,0)
    print(time_1_raster)

    # Get spatial information
    lat = time_1_raster.coords['lat']
    lon = time_1_raster.coords['lon']
    lat_min, lat_max = lat.min().values, lat.max().values
    lon_min, lon_max = lon.min().values, lon.max().values
    rows, cols = len(lat), len(lon)

    transform = from_bounds(lon_min, lat_min, lon_max, lat_max, cols, rows)

    # Define metadata for the output GeoTIFF
    profile = {
        'driver': 'GTiff',
        'height': rows,
        'width': cols,
        'count': 1,
        'dtype': time_1_raster.dtype,
        'crs': 'EPSG:4326',
        'transform': transform
    }

    # Output GeoTIFF file name
    time_1_tif = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{start_year}_{time_1_end}.tif'

    # Write the data slice to a GeoTIFF file
    with rasterio.open(time_1_tif, 'w', **profile) as dst:
        dst.write(time_1_raster.values, 1)
    
    time_1_ensemble_mems.close()
    time_1_ensemble.close()

    #%% Time period 2: 2041-2060
    models_2 = []
    for model in models: 
        annual_cuml_years = []
        for year in range(time_2_start, time_3_start):
            file_path = input_path + f'/{year}/cmip6_rddr_app_water_{ssp}_{model}_{year}.nc'
            ds_year = xr.open_dataset(file_path)
            annual_cuml = ds_year.sum(dim='time')
            annual_cuml_years.append(annual_cuml)
            ds_year.close()
        
        ds_years = xr.concat(annual_cuml_years, dim='year')
        mean_20yrs = ds_years.mean(dim='year')
        models_2.append(mean_20yrs)
        ds_years.close()
    
    time_2_ensemble_mems = xr.concat(models_2, dim='ensembles')
    time_2_ensemble = time_2_ensemble_mems.mean(dim='ensembles')
    
    time_2_save = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_2_start}_{time_2_end}.nc'
    time_2_ensemble.to_netcdf(time_2_save)
    
    time_2_array = time_2_ensemble['app_water']
    time_2_raster = np.flip(time_2_array,0)

    lat = time_2_raster.coords['lat']
    lon = time_2_raster.coords['lon']
    lat_min, lat_max = lat.min().values, lat.max().values
    lon_min, lon_max = lon.min().values, lon.max().values
    rows, cols = len(lat), len(lon)

    transform = from_bounds(lon_min, lat_min, lon_max, lat_max, cols, rows)

    profile = {
        'driver': 'GTiff',
        'height': rows,
        'width': cols,
        'count': 1,
        'dtype': time_2_raster.dtype,
        'crs': 'EPSG:4326',
        'transform': transform
    }

    time_2_tif = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_2_start}_{time_2_end}.tif'

    with rasterio.open(time_2_tif, 'w', **profile) as dst:
        dst.write(time_2_raster.values, 1)
    
    time_2_ensemble_mems.close()
    time_2_ensemble.close()
    
    #%% Time period 3: 2061-2080
    models_3 = []
    for model in models: 
        annual_cuml_years = []
        for year in range(time_3_start, time_4_start):
            file_path = input_path + f'/{year}/cmip6_rddr_app_water_{ssp}_{model}_{year}.nc'
            ds_year = xr.open_dataset(file_path)
            annual_cuml = ds_year.sum(dim='time')
            annual_cuml_years.append(annual_cuml)
            ds_year.close()
        
        ds_years = xr.concat(annual_cuml_years, dim='year')
        mean_30yrs = ds_years.mean(dim='year')
        models_3.append(mean_30yrs)
        ds_years.close()
    
    time_3_ensemble_mems = xr.concat(models_3, dim='ensembles')
    time_3_ensemble = time_3_ensemble_mems.mean(dim='ensembles')
    
    time_3_save = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_3_start}_{time_3_end}.nc'
    time_3_ensemble.to_netcdf(time_3_save)
    
    time_3_array = time_3_ensemble['app_water']
    time_3_raster = np.flip(time_3_array,0)

    lat = time_3_raster.coords['lat']
    lon = time_3_raster.coords['lon']
    lat_min, lat_max = lat.min().values, lat.max().values
    lon_min, lon_max = lon.min().values, lon.max().values
    rows, cols = len(lat), len(lon)

    transform = from_bounds(lon_min, lat_min, lon_max, lat_max, cols, rows)

    profile = {
        'driver': 'GTiff',
        'height': rows,
        'width': cols,
        'count': 1,
        'dtype': time_3_raster.dtype,
        'crs': 'EPSG:4326',
        'transform': transform
    }

    time_3_tif = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_3_start}_{time_3_end}.tif'

    with rasterio.open(time_3_tif, 'w', **profile) as dst:
        dst.write(time_3_raster.values, 1)
    
    time_3_ensemble_mems.close()
    time_3_ensemble.close()

    #%% Time period 4: 2081-2100
    models_4 = []
    for model in models: 
        annual_cuml_years = []
        for year in range(time_4_start, end_range):
            file_path = input_path + f'/{year}/cmip6_rddr_app_water_{ssp}_{model}_{year}.nc'
            ds_year = xr.open_dataset(file_path)
            annual_cuml = ds_year.sum(dim='time')
            annual_cuml_years.append(annual_cuml)
            ds_year.close()
        
        ds_years = xr.concat(annual_cuml_years, dim='year')
        mean_40yrs = ds_years.mean(dim='year')
        models_4.append(mean_40yrs)
        ds_years.close()
    
    time_4_ensemble_mems = xr.concat(models_4, dim='ensembles')
    time_4_ensemble = time_4_ensemble_mems.mean(dim='ensembles')
    
    time_4_save = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_4_start}_{end_year}.nc'
    time_4_ensemble.to_netcdf(time_4_save)
    
    time_4_array = time_4_ensemble['app_water']
    time_4_raster = np.flip(time_4_array,0)

    lat = time_4_raster.coords['lat']
    lon = time_4_raster.coords['lon']
    lat_min, lat_max = lat.min().values, lat.max().values
    lon_min, lon_max = lon.min().values, lon.max().values
    rows, cols = len(lat), len(lon)

    transform = from_bounds(lon_min, lat_min, lon_max, lat_max, cols, rows)

    profile = {
        'driver': 'GTiff',
        'height': rows,
        'width': cols,
        'count': 1,
        'dtype': time_4_raster.dtype,
        'crs': 'EPSG:4326',
        'transform': transform
    }

    time_4_tif = output_path + f'/cmip6_{ssp}_mean_ann_cuml_app_water_{time_4_start}_{end_year}.tif'

    with rasterio.open(time_4_tif, 'w', **profile) as dst:
        dst.write(time_4_raster.values, 1)
    
    time_4_ensemble_mems.close()
    time_4_ensemble.close()