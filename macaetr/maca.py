# -*- coding: utf-8 -*-
"""
Python tools for downloading and using MACA downscaled and bias corrected climate data including the calculation of ASCE reference evapotranspiration. 
"""

import numpy as np
import pandas as pd
import refet
from pathlib import Path



class MACA(object):
    """
    """
    server_prefix = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/'
    # MACAv1 having errors on server
    product_info = {
        'macav2':{
            'server_name':'agg_macav2metdata',
            'resolution':'4-km (1/24-deg)',
            'variables':('tmax','tmin','rel_hum_max','rel_hum_min','precip',
                         'rs','wind_east','wind_north','spec_hum',
                         'vapor_pres_def')
        },
        'livneh':{
            'server_name':'macav2livneh',
            'resolution':'~6-km (1/16-deg)',
            'variables':('tmax','tmin','precip','rs','wind_mean','spec_hum')
        }
    }

    var_info = {
        'tmax': {
            'name':'tmax',
            'internal_name':'air_temperature',
            'server_name':'tasmax',
            'units':'kelvin'
        },
        'tmin': {
            'name':'tmin',
            'internal_name':'air_temperature',
            'server_name':'tasmin',
            'units':'kelvin'
        },
        'spec_hum':{
            'name':'specific_humidity',
            'internal_name':'specific_humidity',
            'server_name':'huss',
            'units':'kg/kg'
        },
        'precip':{
            'name':'precipitation',
            'internal_name':'precipitation',
            'server_name':'pr',
            'units':'mm'
        },
        # wind at 10 m
        'wind_north':{
            'name':'northward_wind',
            'internal_name':'northward_wind',
            'server_name':'vas',
            'units':'m/s'
        },
        # wind at 10 m
        'wind_east':{
            'name':'eastward_wind',
            'internal_name':'eastward_wind',
            'server_name':'uas',
            'units':'m/s'
        },
        'wind_mean':{
            'name':'wind_mean',
            'internal_name':'wind_speed',
            'server_name':'was',
            'units':'m/s'
        },
        'rel_hum_max':{
            'name':'rel_hum_max',
            'internal_name':'relative_humidity',
            'server_name':'rhsmax',
            'units':'%'
        },
        'rel_hum_min':{
            'name':'rel_hum_min',
            'internal_name':'relative_humidity',
            'server_name':'rhsmin',
            'units':'%'
        },
        'vapor_pres_def':{
            'name':'vapor_pres_def',
            'internal_name':'vpd',
            'server_name':'vpd',
            'units':'%'
        },
        'rs':{
            'name':'rs',
            'internal_name':'surface_downwelling_shortwave_flux_in_air',
            'server_name':'rsds',
            'units':'w/m2'
        }

    }

    models = (
        'bcc-csm1-1',
        'bcc-csm1-1-m',
        'BNU-ESM',
        'CanESM2',
        'CCSM4',
        'CNRM-CM5',
        'CSIRO-Mk3-6-0',
        'GFDL-ESM2G',
        'GFDL-ESM2M',
        'HadGEM2-CC365',
        'HadGEM2-ES365',
        'inmcm4',
        'IPSL-CM5A-MR',
        'IPSL-CM5A-LR',
        'IPSL-CM5B-LR',
        'MIROC5',
        'MIROC-ESM',
        'MIROC-ESM-CHEM',
        'MRI-CGCM3',
        'NorESM1-M'
    )
    scenarios = ('rcp45', 'rcp85', 'historical')
    date_limits = {
        'rcp45':(2006,2099),
        'rcp85':(2006,2099),
        'historical':(1950,2005)
    }

    def __init__(self):
        pass


    def _check_arguments(self,variable, date_start, date_end, product, model, 
            scenario):
        """Check if MACA selection options are valid"""
        if model not in MACA.models:
            print(
                f'ERROR: {model} is not a valid climate model option, '
                f'pick from:\n{MACA.models}'
            )
            return False
        if product not in MACA.product_info.keys():
            print(
                f'ERROR: {product} is not a valid MACA product option, '
                f'pick from:\n{tuple(MACA.product_info.keys())}'
            )
            return False
        if variable not in MACA.product_info.get(product,{}).get('variables'):
            print(
                f'ERROR: {variable} is not a valid climate variable for the {product} dataset,'
                f' pick from:\n{MACA.product_info.get(product,dict()).get("variables")}'
            )
            return False
        if scenario not in MACA.scenarios:
            print(
                f'ERROR: {scenario} is not a valid emission scenario, pick from:\n'
                f' {MACA.scenarios}'
            )
            return False
        # check time period
        start = pd.to_datetime(date_start).year
        valid_start = MACA.date_limits.get(scenario)[0]
        end = pd.to_datetime(date_end).year
        valid_end = MACA.date_limits.get(scenario)[1]

        if start < valid_start:
            print(
                f'ERROR: date {date_start} preceeds the {scenario} time period, time period'
                f' begins on {valid_start}'
            )
            return False

        if end > valid_end:
            print(
                f'ERROR: date {date_end} exceeds the {scenario} time period, time period'
                f' ends on {valid_end}'
            )
            return False
        else: # passes all checks
            return True

    def _make_server_url(self, variable, product, model, scen):
        """
        Construct url string for netCDF for OpenDaP server for given selection.
        Does NOT check if selection options are valid.
        """
        pref = MACA.server_prefix
        prod = self.product_info.get(product).get('server_name')
        variable = self.var_info.get(variable).get('server_name')
        model = f'{model}_r1i1p1'
        start = MACA.date_limits.get(scen)[0]
        end = MACA.date_limits.get(scen)[1]
        suff = 'CONUS_daily.nc' if product == 'macav2' else 'CONUS_daily_aggregated.nc'

        return f'{pref}{prod}_{variable}_{model}_{scen}_{start}_{end}_{suff}'

    def download(self, lat, lon, start_date, end_date, variable, 
            model='GFDL-ESM2G', product='macav2', scenario='rcp85',
            ret='np_array'):
        """
        """
        ret_options = ('np_array','xarray','pd_series','pd_dataframe')
        valid = self._check_arguments(variable, start_date, end_date, 
                product, model, scenario)
        if not valid:
            print('Aborting download!')
            return
        server_url = self._make_server_url(variable, product, model, scenario)

        lon += 360 # lats are positive in MACA datasets
        ds = xarray.open_dataset(server_url).sel(
            time=pd.date_range(start_date,end_date),
            lon=lon,
            lat=lat,
            method='nearest'
        ).drop('crs')

        self.data = None
        self.metadata = ds.attrs
        self.coordinates = ds.coords
        self.var_name = self.var_info.get(variable).get('name')
        self.var_units = self.var_info.get(variable).get('units')
        internal_name = f'{self.var_info.get(variable).get("internal_name")}'
        if ret == 'np_array':
            self.data = ds[internal_name].data
        elif ret == 'xarray':
            self.data = ds
        elif ret == 'pd_series':
            self.data = pd.Series(
                index=pd.date_range(start_date,end_date),
                data=ds[internal_name].data
            )
            self.data.name = self.var_name
            self.data.index.name = 'date'
        elif ret == 'pd_dataframe':
            self.data = ds.drop(['lat','lon']).to_dataframe().rename(
                columns={internal_name:self.var_name}
            )
            self.data.index.name = 'date'
        else:
            print(
                f'WARNING: {ret} is an invalid return type option, returning the default'
                f' numpy array, pick from:\n{ret_options}'
            )
            self.data = ds[f'{self.var_info.get(variable).get("internal_name")}'].data


        return self.data

    def daily_ETr(self, reference="tall"):
        """
        """
        pass

