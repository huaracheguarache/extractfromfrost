#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
import logging.handlers
import lxml.etree as ET
from netCDF4 import Dataset
import numpy

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--dest",dest="destination",
            help="Destination which to add NCML to", required=True)
    parser.add_argument("-l","--log",dest="logdir",
            help="Destination where to put logfiles", required=True)
    parser.add_argument("-o","--overwrite",action='store_true',
            help="Overwrite if NCML is existing")
    args = parser.parse_args()

    """
    if args.cfgfile is False:
        parser.print_help()
        parser.exit()
    """

    return args

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

def traverse_structure(myfolder):

    for item in os.listdir(myfolder):
        mydir = '/'.join([myfolder,item])
        if not os.path.isdir(mydir):
            continue
        if not item.startswith('SN'):
            mylog.warn('Apparently this is not a station folder.')
            continue
        mylog.info('Processing folder: %s', mydir)
        myncmlfile = '/'.join([mydir,item])+'-aggregated.ncml'
        try:
            create_ncml(myncmlfile, mydir)
        except Exception as e:
            mylog.error('Something failed.')
            raise

def create_ncml(myncmlfile, aggdir):
    # Check if NCML already exist
    if os.path.isfile(myncmlfile) and not args.overwrite:
        mylog.warning("%s already exists, won't do anything", myncmlfile)
        return

    # Create NCML
    ns_map = {None:'http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2'}
    root = ET.Element(ET.QName('netcdf'), nsmap = ns_map)
    aggel = ET.SubElement(root, ET.QName('aggregation'))
    aggel.set('dimName','time')
    aggel.set('type','joinExisting')
    aggel.set('recheckEvery','1 day')
    """
    Removed in favor of listing fo files
    scanel = ET.SubElement(aggel, ET.QName('scan'))
    scanel.set('location', aggdir)
    scanel.set('suffix','.nc')
    """
    # Set up specific files to include, assuming data stored in years
    for item in os.listdir(aggdir):
        curdir = '/'.join([aggdir,item])
        # Check content of yearly folder
        if os.path.isdir(curdir):
            # Process files
            for item2 in sorted(os.listdir(curdir)):
                myfile = '/'.join([curdir,item2])
                # Open NetCDF and check content
                if myfile.endswith('.nc'):
                    myncds = Dataset(myfile)
                    tmp = myncds.variables['time'][:]
                    tmpstring = ' '.join(str(num) for num in tmp)
                    netcdf = ET.SubElement(aggel, ET.QName('netcdf'))
                    netcdf.set('location',myfile)
                    netcdf.set('coordValue', tmpstring)

    # Dump NCML file
    et = ET.ElementTree(root)
    et.write(myncmlfile, xml_declaration=True, encoding='UTF-8', pretty_print=True)

if __name__ == '__main__':
    
    # Parse command line arguments
    try:
        args = parse_arguments()
    except:
        raise SystemExit('Command line arguments didn\'t parse correctly.')

    # Parse configuration file
    #cfgstr = parse_cfg(args.cfgfile)

    # Initialise logging
    #output_dir = cfgstr['output']
    mylog = initialise_logger(args.logdir)
    mylog.info('Configuration of logging is finished.')

    try:
        traverse_structure(args.destination)
    except Exception as e:
        mylog.error('Something failed %s', e)
