#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Checking CF-NetCDF files in a structure to ensure that all files have the same number of variables and allows aggregation in time.
"""
__author__ = "Øystein Godøy"
#__copyright__ = "Copyright Info"
__credits__ = ["Øystein Godøy",]
__license__ = """
    This file is part of extractfromfrost.

    This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This file is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
    """
__version__ = "0.0.0"
__maintainer__ = "Øystein Godøy"
__email__ = "steingod@met.no"
#__status__ = "status of the file"

import os
import sys
import argparse
import logging
import logging.handlers
from netCDF4 import Dataset
import numpy
import pytz
from datetime import datetime

def parse_arguments():
    """
    Set up the command line interface.
    """
    parser = argparse.ArgumentParser(description="Checker for CF-NetCDF files to ensure consistency over files allowing aggregation in time.",epilog="")
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
    """
    Set up the logging for the application.
    """
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
    """
    Assuming data from one station is organised in a single folder and with sub folders for each year. This function loops through all stations.
    """

    for item in os.listdir(myfolder):
        mydir = '/'.join([myfolder,item])
        if not os.path.isdir(mydir):
            continue
        if not item.startswith('SN'):
            mylog.warn('Apparently this is not a station folder.')
            continue
        mylog.info('Processing folder: %s', mydir)
        try:
            check_netcdf(mydir)
        except Exception as e:
            mylog.error('Something failed.')
            raise

def check_netcdf(stdir):
    """
    Check the individual files in the folder for each station. If some files miss variables that have been added later, these are added to the respective files and set to all missing values. This function works backwards under the assumption that there is a larger probability for variables to be added than removed.
    """

    # Loop through folder
    for item in sorted(os.listdir(stdir), reverse=True):
        # Check content of yearly folder
        curdir = '/'.join([stdir,item])
        myvariables = list()
        if os.path.isdir(curdir):
            # Process files for each year and extract list of variables
            for item2 in sorted(os.listdir(curdir), reverse=True):
                myfile = '/'.join([curdir,item2])
                # Only process NetCDF files and check content
                if myfile.endswith('.nc'):
                    mylog.info('Processing file: %s', myfile)
                    myncds = Dataset(myfile)
                    tmpvars = list(myncds.variables.keys())
                    if len(myvariables) == 0:
                        myvariables = tmpvars
                    if not bool(set(myvariables).intersection(tmpvars)):
                        mylog.warning('This file has different variables than others\nChecker: %s\nFile: %s',myvariables,tmpvars)
                        myvariables = tmpvars
                        sys.exit()
                    else:
                        mylog.info('This file has the same variables as other files, continuing.')
                    myncds.close()

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
    mylog.info('Starting new check of CF-NetCDF files consistency...')

    try:
        traverse_structure(args.destination)
    except Exception as e:
        mylog.error('Something failed %s', e)
