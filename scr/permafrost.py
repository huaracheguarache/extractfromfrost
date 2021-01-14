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

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c","--cfg",dest="cfgfile",
            help="Configuration file", required=True)
    parser.add_argument("-s","--startday",dest="startday",
            help="Start day in the form YYYY-MM-DD", required=True)
    parser.add_argument("-e","--endday",dest="endday",
            help="End day in the form YYYY-MM-DD", required=True)
    args = parser.parse_args()
    print(args)
    print(args.cfgfile)
    print(args.startday)
    print(args.endday)

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

def extractdata(frostcfg,station,stmd,dumplocation,abstract):

    # Create request for observations
    myrequest = 'ids='+station
    # Connect and read metadata
    r = requests.get(frostcfg['endpointmeta'],
            myrequest,
            auth=(frostcfg['client_id'],""))
    # Check if the request worked, print out any errors
    if r.status_code != 200:
        print('Error! Returned status code %s' % r.status_code)
        print(r.text)
        sys.exit()

    print(r.text)
    metadata = json.loads(r.text)

    # Create request for observations
    myrequest = ('sources='+station+'&elements='
        +'.'.join(frostcfg['elements'])
        +'&fields='+','.join(frostcfg['fields'])
        +'&referencetime='+'/'.join([args.startday,args.endday]))
    # Connect and read observations
    r = requests.get(frostcfg['endpointobs'],
            myrequest,
            auth=(frostcfg['client_id'],""))
    # Check if the request worked, print out any errors
    if r.status_code != 200:
        print('Error! Returned status code %s' % r.status_code)
        print(r.text)
        sys.exit()
    # Read into  Pandas DataFrame
    df = pd.read_csv(StringIO(r.text),header=0,
        mangle_dupe_cols=True, parse_dates=['referenceTime'],
        index_col=False,na_values=['-'])

    # Create data frame for each time step (i.e. a profile)
    #print(list(df.columns))
    ntime = 0
    mytimes = []
    myprofiles = list()
    mydepths = []
    for i in df.loc[:,'referenceTime']:
        #print(type(i))
        mytimes.append(i)
        mydata = {
                'depth':df.filter(like='depth_below_surface',axis='columns').iloc[ntime].values,
                'temperature':df.filter(like='soil_temperature',axis='columns').iloc[ntime].values
                }
        mytmpdata = pd.DataFrame(mydata).sort_values(by='depth')
        myprofiles.append(mytmpdata.temperature)
        mydepths = mytmpdata.depth
        ntime += 1

    datasetstart = min(mytimes).strftime('%Y-%m-%dT%H:%M:%SZ')
    datasetend = max(mytimes).strftime('%Y-%m-%dT%H:%M:%SZ')
    mytimes = (pd.to_datetime(mytimes,
        utc=True)-pd.Timestamp("1970-01-01", tz='UTC')) // pd.Timedelta('1s')
    da_profile = xr.DataArray(myprofiles, 
            dims=['time','depth'],
            coords={
                'time':mytimes,
                'depth':mydepths})
    da_profile.name = 'soil_temperature'
    da_profile.attrs['name'] = 'temperature'
    da_profile.attrs['standard_name'] = 'soil_temperature'
    da_profile.attrs['units'] = 'degree_Celsius'
    da_profile.depth.attrs['name'] = 'depth'
    da_profile.depth.attrs['standard_name'] = 'depth'
    da_profile.depth.attrs['long_name'] = 'depth below surface in centimeters'
    da_profile.depth.attrs['units'] = 'centimeter'
    da_profile.depth.attrs['positive'] = 'down'
    da_profile.time.attrs['standard_name'] = 'time'
    da_profile.time.attrs['units'] = 'seconds since 1970-01-01 00:00:00+0'

    # Need to convert from dataarray to dataset in order to add global
    # attributes
    ds_profile = da_profile.to_dataset()
    ds_profile.attrs['title'] = 'Permafrost borehole measurements at '+stmd['name']
    ds_profile.attrs['summary'] = abstract
    ds_profile.attrs['license'] = metadata['license']
    ds_profile.attrs['time_coverage_start'] = datasetstart
    ds_profile.attrs['time_coverage_end'] = datasetend
    ds_profile.attrs['geospatial_lat_min'] = metadata['data'][0]['geometry']['coordinates'][1]
    ds_profile.attrs['geospatial_lat_max'] = metadata['data'][0]['geometry']['coordinates'][1]
    ds_profile.attrs['geospatial_lon_min'] = metadata['data'][0]['geometry']['coordinates'][0]
    ds_profile.attrs['geospatial_lon_max'] = metadata['data'][0]['geometry']['coordinates'][0]
    ds_profile.attrs['creator_name'] = stmd['PrincipalInvestigator'] 
    ds_profile.attrs['creator_email'] = stmd['PrincipalInvestigatorEmail']
    ds_profile.attrs['creator_url'] = stmd['PrincipalInvestigatorOrganisationURL']
    ds_profile.attrs['creator_institution'] = stmd['PrincipalInvestigatorOrganisation']
    ds_profile.attrs['keywords'] = 'Earth Science > Cryosphere > Frozen Ground > Permafrost > Permafrost Temperature'
    ds_profile.attrs['keywords_vocabulary'] = 'GCMD'
    ds_profile.attrs['publisher_name'] = ''
    ds_profile.attrs['publisher_email'] = 'adc@met.no'
    ds_profile.attrs['publisher_url'] = 'https://adc.met.no/'
    ds_profile.attrs['publisher_institution'] = 'Norwegian Meteorlogical Institute'
    ds_profile.attrs['Conventions'] = 'ACDD, CF-1.8'
    ds_profile.attrs['date_created'] = metadata['createdAt']
    ds_profile.attrs['history'] = metadata['createdAt']+': Data extracted from the MET Observation Database through Frost and stored as NetCDF-CF'
    ds_profile.attrs['source'] = 'Soil temperature from permafrost boreholes'
    ds_profile.attrs['wigosId'] = metadata['data'][0]['wigosId']
    ds_profile.attrs['project'] = stmd['Project']

    # Plotting  works best on DataArray, not sure how to do on DataSet, i.e.
    # extract DataArray from Dataset first
    #da_profile.plot(aspect=0.80,size=6,y='depth',yincrease=False)
    #plt.savefig('myfig.png')
    #plt.show()

    # Dump to Netcdf
    #print(ds_profile)
    ds_profile.to_netcdf(dumplocation,
            encoding={'depth': {'dtype':'int32'},
                'time': {'dtype': 'int32'},
                'soil_temperature': {'dtype': 'float32'}
                })
    return

if __name__ == '__main__':
    
    # Parse command line arguments
    args = parse_arguments()
    print(args.cfgfile)

    # Parse configuration file
    cfgstr = parse_cfg(args.cfgfile)
    #print(cfgstr)
    #print(cfgstr['frostcfg']['client_id'])

    # Loop through stations
    for station,content in cfgstr['stations'].items():
        if station in ['SN99380','SN99927','SN99879']:
            continue
        print('Requesting data for', station)
        outputfile = cfgstr['output']['destdir']+'/'+content['filename']+'.nc'
        print(outputfile)
        extractdata(cfgstr['frostcfg'], station, content,
                outputfile,cfgstr['output']['abstract'])
        sys.exit()
