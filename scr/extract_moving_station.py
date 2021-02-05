#!/usr/bin/python3

import os
import sys
import argparse
#import csv
import requests
import pandas as pd
from io import StringIO
import xarray as xr
import matplotlib.pyplot as plt
#import datetime as dt
import json
import yaml
import logging
from logging.handlers import TimedRotatingFileHandler
from datetime import datetime
import pprint

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c","--cfg",dest="cfgfile",
            help="Configuration file", required=True)
    parser.add_argument("-s","--startday",dest="startday",
            help="Start day in the form YYYY-MM-DD", required=True)
    parser.add_argument("-e","--endday",dest="endday",
            help="End day in the form YYYY-MM-DD", required=True)
    args = parser.parse_args()

    try:
        datetime.strptime(args.startday,'%Y-%m-%d')
    except ValueError:
        raise ValueError
    try:
        datetime.strptime(args.endday,'%Y-%m-%d')
    except ValueError:
        raise ValueError

    if args.cfgfile is None:
        parser.print_help()
        parser.exit()

    return args

def parse_cfg(cfgfile):
    # Read config file
    print("Reading", cfgfile)
    with open(cfgfile, 'r') as ymlfile:
        cfgstr = yaml.full_load(ymlfile)

    return cfgstr

def initialise_logger(outputfile = './log'):
    # Check that logfile exists
    logdir = os.path.dirname(outputfile)
    if not os.path.exists(logdir):
        try:
            os.makedirs(logdir)
        except:
            raise IOError
    # Set up logging
    mylog = logging.getLogger()
    mylog.setLevel(logging.INFO)
    #logging.basicConfig(level=logging.INFO, 
    #        format='%(asctime)s - %(levelname)s - %(message)s')
    myformat = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(myformat)
    mylog.addHandler(console_handler)
    file_handler = logging.handlers.TimedRotatingFileHandler(
            outputfile,
            when='w0',
            interval=1,
            backupCount=7)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(myformat)
    mylog.addHandler(file_handler)

    return(mylog)

def getships(frostcfg):

    myrequest = {'municipality': 'skip'}

    # Connect and read information
    try:
        r = requests.get(frostcfg['endpointsources'],
                myrequest,
                auth=(frostcfg['client_id'],""))
    except:
        mylog.error('Something went wrong extracting metadata.')
        raise
    # Check if the request worked, print out any errors
    if not r.ok:
        mylog.error('Returned status code was %s', r.status_code)
        print(r.text)
        raise
    mylist = r.json()
    return(mylist)

def extractdata(frostcfg,station,stmd,output):

    # Connect and read metadata for the station
    mylog.info('Retrieving metadata for station: %s', station)
    myrequest = {'ids': station}
    r = requests.get(frostcfg['endpointsources'],
            myrequest,
            auth=(frostcfg['client_id'],""))
    # Check if the request worked, print out any errors
    if not r.ok:
        mylog.error('Returned status code was %s with response %s', r.status_code, r.text)
        raise
    metadata = json.loads(r.text)

    # Connect and read information on available parameters at the station
    mylog.info('Retrieving parameters for station: %s', station)
    myrequest = {'sources': station}
    r = requests.get(frostcfg['endpointparameters'],
            myrequest,
            auth=(frostcfg['client_id'],""))
    # Check if the request worked, print out any errors
    if not r.ok:
        mylog.error('Returned status code was %s, message was:\n%s', r.status_code, r.text)
        raise
    parameterlist = json.loads(r.text)
    pprint.pprint(parameterlist)
    sys.exit()
    # Create list of parameters to extract, latitude and olongitude are automatically added. Only checking for observations, i.e. data within hours and minutes, not daily or monthly aggregates.
    myparameters = ''
    i = 0
    for item in parameterlist['data']:
        if 'PT' in item['timeResolution']:
            print(item['timeResolution']+' - '+item['elementId'])
            if i == 0:
                myparameters = item['elementId']
            else:
                myparameters += ','+item['elementId']
            i+=1 

    # Create request for observations
    mylog.info('Retrieving data for station: %s', station)
    myrequest = {
            'sources': station,
            'elements': myparameters,
            'fields': ','.join(frostcfg['fields']),
            'referencetime': '/'.join([args.startday,args.endday])
            }
    # Connect and read observations
    try:
        r = requests.get(frostcfg['endpointobs'],
                myrequest,
                auth=(frostcfg['client_id'],""))
    except:
        mylog.error('Something went wrong extracting data.')
        raise
    # Check if the request worked, print out any errors
    if r.status_code == 412:
        mylog.error('Information returned indicates that no data is available for this time period for station %s', station)
        return
    if not r.status_code == 200:
        mylog.error('Returned status code was %s\nmessage:\n%s', r.status_code, r.text)
        raise
    # Read into  Pandas DataFrame, assuming - is used for missing values.
    df = pd.read_csv(StringIO(r.text),header=0,
        mangle_dupe_cols=True, parse_dates=['referenceTime'],
        index_col=False,na_values=['-'])

    columnnames = list(df.columns)
    columnnames = [s.replace('(-)','') for s in list(df.columns)]
    columnnames[columnnames.index('referenceTime')] = 'time'
    df.columns = columnnames

    timemin = min(df['time'])
    timemax = max(df['time'])
    datasetstart = timemin.strftime('%Y-%m-%dT%H:%M:%SZ')
    datasetend = timemax.strftime('%Y-%m-%dT%H:%M:%SZ')
    datasetstart4filename = timemin.strftime('%Y%m%d')
    datasetend4filename = timemax.strftime('%Y%m%d')
    mytimes = (pd.to_datetime(df['time'], utc=True)-pd.Timestamp("1970-01-01", tz='UTC')) // pd.Timedelta('1s')
    # Make sure time is dimension and coordinate variable
    ds_timeseries = df.set_index(['time']).to_xarray()
##    da_timeseries.name = 'soil_temperature'
##    da_timeseries.attrs['name'] = 'temperature'
##    da_timeseries.attrs['standard_name'] = 'soil_temperature'
##    da_timeseries.attrs['units'] = 'degree_Celsius'
##    da_timeseries.depth.attrs['name'] = 'depth'
##    da_timeseries.depth.attrs['standard_name'] = 'depth'
##    da_timeseries.depth.attrs['long_name'] = 'depth below surface in centimeters'
##    da_timeseries.depth.attrs['units'] = 'centimeter'
##    da_timeseries.depth.attrs['positive'] = 'down'
##    da_timeseries.time.attrs['standard_name'] = 'time'
##    da_timeseries.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00+0'

    # Need to convert from dataarray to dataset in order to add global attributes
    ds_timeseries.attrs['featureType'] = 'timeSeries'
    ds_timeseries.attrs['title'] = 'Weather station information from ship '+stmd['name']
    ds_timeseries.attrs['summary'] = output['abstract']
    ds_timeseries.attrs['license'] = metadata['license']
    ds_timeseries.attrs['time_coverage_start'] = datasetstart
    ds_timeseries.attrs['time_coverage_end'] = datasetend
    ds_timeseries.attrs['geospatial_lat_min'] = min(ds_timeseries.data_vars['latitude'].values)
    ds_timeseries.attrs['geospatial_lat_max'] = max(ds_timeseries.data_vars['latitude'].values)
    ds_timeseries.attrs['geospatial_lon_min'] = min(ds_timeseries.data_vars['longitude'].values)
    ds_timeseries.attrs['geospatial_lon_max'] = max(ds_timeseries.data_vars['longitude'].values)
    ds_timeseries.attrs['creator_name'] = stmd['PrincipalInvestigator'] 
    ds_timeseries.attrs['creator_email'] = stmd['PrincipalInvestigatorEmail']
    ds_timeseries.attrs['creator_url'] = stmd['PrincipalInvestigatorOrganisationURL']
    ds_timeseries.attrs['creator_institution'] = stmd['PrincipalInvestigatorOrganisation']
    # Remember to fix these FIXME
    ds_timeseries.attrs['keywords'] = 'Earth Science > Cryosphere > Frozen Ground > Permafrost > Permafrost Temperature,Earth Science > Land Surface > Soils > Soil temperature'
    ds_timeseries.attrs['keywords_vocabulary'] = 'GCMD'
    ds_timeseries.attrs['publisher_name'] = ''
    ds_timeseries.attrs['publisher_email'] = 'adc@met.no'
    ds_timeseries.attrs['publisher_url'] = 'https://adc.met.no/'
    ds_timeseries.attrs['publisher_institution'] = 'Norwegian Meteorlogical Institute'
    ds_timeseries.attrs['Conventions'] = 'ACDD, CF-1.8'
    ds_timeseries.attrs['date_created'] = metadata['createdAt']
    ds_timeseries.attrs['history'] = metadata['createdAt']+': Data extracted from the MET Observation Database through Frost and stored as NetCDF-CF'
    ds_timeseries.attrs['source'] = 'Soil temperature from permafrost boreholes'
    ds_timeseries.attrs['wigosId'] = metadata['data'][0]['wigosId']
    ds_timeseries.attrs['METNOId'] =  station
    ds_timeseries.attrs['project'] = stmd['Project']

    ds_timeseries.encoding['unlimited_dims'] = 'time'

    pprint.pprint(metadata['data'])
    print(metadata['data'][0]['name'])
    print(metadata['data'][0]['stationHolders'])
    print(ds_timeseries)

    # Dump to Netcdf
    outputfile = output['destdir']+'/ship-'+metadata['data'][0]['wigosId']+'_'+datasetstart4filename+'-'+datasetend4filename+'.nc'
    mylog.info('Dumping data to NetCDF-CF:\n%s', outputfile)
    try:
        ds_timeseries.to_netcdf(outputfile)
    except:
        mylog.error('Creation of NetCDF file didn\'t work properly\n%s', sys.exc_info()[0])
        raise
##            encoding={'depth': {'dtype':'int32'},
##                'time': {'dtype': 'int32'},
##                'soil_temperature': {'dtype': 'float32'}
##                })
    return

if __name__ == '__main__':
    
    # Parse command line arguments
    try:
        args = parse_arguments()
    except:
        raise SystemExit('Command line arguments didn\'t parse correctly.')

    # Parse configuration file
    cfgstr = parse_cfg(args.cfgfile)

    # Initialise logging
    mylog = initialise_logger(cfgstr['output']['logfile'])
    mylog.info('Configuration of logging is finished.')

    # Find all ships available
    mylog.info('Retrieve all ships available in the data storage.')
    try:
        ships = getships(cfgstr['frostcfg'])
    except:
        mylog.warn('Couldn\'t get the list of ships in data storage.')
        raise SystemExit()
    pprint.pprint(ships)

    # Loop through ships
    mylog.info('Processing ships from the list.')
    for item in ships['data']:
        mylog.info('Extracting information for %s - %s',item['name'], item['id'])
        try:
            extractdata(cfgstr['frostcfg'], item['id'], cfgstr['stations']['SN15270'], cfgstr['output'])
        except:
            mylog.error('Something went horrible wrong here.')
            raise SystemExit()
