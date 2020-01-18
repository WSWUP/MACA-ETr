# -*- coding: utf-8 -*-
"""
Python tools for downloading and using MACA downscaled and bias corrected climate data including the calculation of ASCE reference evapotranspiration. 
"""
import pkg_resources
import numpy as np
import pandas as pd
import refet
import xarray
from scipy import spatial

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
            'units':'kPa'
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
        self.data = None
        self.coordinates = None


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
        model = f'{model}_r1i1p1' if not model == 'CCSM4' else f'{model}_r6i1p1'
        start = MACA.date_limits.get(scen)[0]
        end = MACA.date_limits.get(scen)[1]
        suff = 'CONUS_daily.nc#fillmismatch' if product == 'macav2' else 'CONUS_daily_aggregated.nc#fillmismatch'

        return f'{pref}{prod}_{variable}_{model}_{scen}_{start}_{end}_{suff}'

    def download(self, lat, lon, start_date, end_date, variable, 
            model='GFDL-ESM2G', product='macav2', scenario='rcp85',
            ret='pd_dataframe'):
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

        self.scenario = scenario
        self.model = model
        self.product = product

        # attrs will get overwritten with subsequent downloads
        self.metadata = ds.attrs
        if self.coordinates:
            t1=self.coordinates.get('lon').values != ds.coords.get('lon').values
            t2=self.coordinates.get('lat').values != ds.coords.get('lat').values
            if t1 or t2 and ret == 'pd_dataframe':
                print('Downloading a new location, creating a new dataframe.')
                self.data = None
        self.coordinates = ds.coords
        self.centroid_lat = float(ds.coords.get('lat').values)
        self.centroid_lon = float(ds.coords.get('lon').values) - 360 # dec. deg.
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
            if isinstance(self.data, pd.DataFrame):
                self.data[self.var_name] = ds.drop(
                    ['lat','lon']).to_dataframe().rename(
                        columns={internal_name:self.var_name}
                )
            else:
                self.data = ds.drop(['lat','lon']).to_dataframe().rename(
                    columns={internal_name:self.var_name}
                )
            self.data.index.name = 'date'
        else:
            print(
                f'WARNING: {ret} is an invalid return type option, returning the default'
                f' numpy array, pick from:\n{ret_options}'
            )
            self.data = ds[
                f'{self.var_info.get(variable).get("internal_name")}'
            ].data


        return self.data

    def _get_elev(self, lat, lon, product):

        inf = f'metadata/{product}_cell_data.csv'
        try:
            if pkg_resources.resource_exists('macaetr', inf):
                meta_path = pkg_resources.resource_filename('macaetr', inf)
        except:
            print(
                f'ERROR: could not find {product} gridcell metadata file '
                'please provide elevation and rerun'
            )
            return
        df = pd.read_csv(meta_path)
        pts = list(zip(df.LAT,df.LON)) # grid centroid locations
        tree = spatial.KDTree(pts)
        index = tree.query(np.array([lat,lon]))[1]
        return df.loc[index, 'ELEV']


    def daily_refet(self, lat, lon, start_date, end_date, anemometer_height=10,
            elev=None, model='GFDL-ESM2G', product='macav2', scenario='rcp85'):
        """
        Download MACAv2 or livneh data required and calculate short and tall
        (abbreviated ETo and ETr) ASCE standardized reference 
        evapotranspiration.

        """

        if product == 'livneh':
            get_vars = ['tmax', 'tmin','rs','wind_mean','spec_hum']
        elif product == 'macav2':
            get_vars = [
                'tmax','tmin','rs','wind_east','wind_north', 'vapor_pres_def'
            ]
        else:
            print(f'ERROR: {product} not a valid dataset product, aborting.')
            return
        # wipe out any existing instance data to avoid issues
        self.data = None

        for v in get_vars:
            self.download(lat, lon, start_date, end_date, v, 
                model=model, product=product, scenario=scenario,
                ret='pd_dataframe'
            )

        self.data['tavg'] = self.data[['tmin','tmax']].mean(1)
        # convert temps to celcius
        self.data[['tmin_c','tmax_c','tavg_c']] =\
            self.data[['tmin','tmax','tavg']] - 273.15
        
        if elev is None:
            # get it from metadata based on product
            elev = self._get_elev(lat, lon, product)

        self.elev = elev
        length = len(self.data)

        if product == 'macav2':
            self.data['wind_mean'] = self.data[
                ['eastward_wind','northward_wind']
            ].abs().mean(1) # mean of absolute wind velocities
            # estimate sat. vapor press
            es = 0.6108*np.exp(17.27*self.data.tavg_c/(self.data.tavg_c+237.3))
            self.data['ea_kpa'] = self.data.vapor_pres_def + es

        elif product == 'livneh':
            self.data['pa'] =  np.full(length, refet.calcs._air_pressure(elev))
            self.data['ea_kpa'] = refet.calcs._actual_vapor_pressure(
                self.data.specific_humidity, self.data.pa
            )

        tmin = self.data.tmin_c
        tmax = self.data.tmax_c
        rs = self.data.rs
        ea = self.data.ea_kpa
        uz = self.data.wind_mean
        zw = anemometer_height 
        lats = np.full(length, self.centroid_lat)
        doy = tmin.index.dayofyear
        elevs = np.full(length, elev)

        # refet converts rs units, everything else in target units
        input_units = {'rs': 'w/m2'}
        REF = refet.Daily(
            tmin, tmax, ea, rs, uz, zw, elevs, lats, doy, method='asce',
            input_units=input_units
        )

        self.data['ETr_mm'] = REF.etr() 
        self.data['ETo_mm'] = REF.eto() 




