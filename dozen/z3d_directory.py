#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 20:39:17 2018

@author: andymcaliley

Functions to quickly parse z3d-related metadata from all files in a directory
"""
# %gui tk
# must run above line IN CONSOLE first if in Spyder
import os
from fnmatch import fnmatch
from math import pi

import numpy as np
import pandas as pd
from pyproj import Proj

from . import z3dio
from . import timeio
from . import util

datetimeformat = "%Y-%m-%d %H:%M:%S"

initialdir = './'
cal_dir = './'


def directory_info(metadata_csv='metadata.csv',use_csv=True,initialdir=initialdir,ask_dir=True):
    '''
    read metadata about z3d files in directory
    read from csv file if available
    otherwise, read from directory
    '''
    if use_csv:
        try:
            dirinfo = load_directory_info(metadata_csv)
        except:
            dirinfo = get_z3d_directory_info(initialdir=initialdir,ask_dir=ask_dir)
            dirinfo.to_csv(metadata_csv)
    else:
        dirinfo = get_z3d_directory_info(initialdir=initialdir,ask_dir=ask_dir)
    return dirinfo


def load_directory_info(filepath):
    '''
    Read metadata summary (rx or tx)
    '''
    df = pd.read_csv(filepath,index_col=0)
    # get start and end as dates
    no_records = (df['num_records']==0)
    if no_records.sum()>0:
        df.loc[~no_records,'start'] = pd.to_datetime(df.loc[~no_records,'start'],utc=True).dt.tz_convert('US/Mountain')
        df.loc[~no_records,'end'] = pd.to_datetime(df.loc[~no_records,'end'],utc=True).dt.tz_convert('US/Mountain')
        df.loc[no_records,'start'] = -1
        df.loc[no_records,'end'] = -1
        df.loc[no_records,'valid'] = False
    else:
        df.loc[:,'start'] = pd.to_datetime(df.start,utc=True).dt.tz_convert('US/Mountain')
        df.loc[:,'end'] = pd.to_datetime(df.end,utc=True).dt.tz_convert('US/Mountain')
    return df


def get_z3d_directory_info(initialdir=initialdir,ask_dir=True):
    '''
    Get select metadata about each z3d file in a directory and its subdirectories
    '''
    if ask_dir:
        rootdir = util.getdir(initialdir=initialdir,title='Select Z3D directory',mustexist=True)
    else:
        rootdir = initialdir

    # read metadata from all z3d files in root directory structure
    # TODO: implement summary-fileinfo keyword pairs, based on ZenAcq version
    # utm13n = Proj(init='epsg:4326',proj="utm",zone=13)
    utm13n = Proj(proj="utm",zone=13,ellps='WGS84')
    pattern = "*.z3d"
    summaries = []
    for path, subdirs, files in os.walk(rootdir):
        for name in files:
            if fnmatch(name.lower(), pattern):
                filepath = os.path.join(path, name)
                summary = {}
                summary['valid']=True
                summary['error']=''
                summary['filename']=name
                summary['filepath']=path
                summary['fullpath']=os.path.join(path,name)
                try:
                    file_info = z3dio.get_file_info(filepath)
                    metadata = file_info['metadata']
                    summary['sampling_rate'] = metadata.get('A/D Rate')
                    lat = metadata.get('Lat')*180/pi
                    lon = metadata.get('Long')*180/pi
                    summary['easting'],summary['northing'] = utm13n(lon,lat)
                    # summary['easting'],summary['northing'] = lon,lat
                    summary['num_records'] = file_info['num_records']
                    # schedule_info = {}
                    # if 'Tx.Freq' in metadata[0]:
                    if 'TX.SENSE' in metadata:
                        # NOTE: this only ids channel 1 TX files
                        summary['type'] = 'TX'
                        summary['tx_sn'] = metadata.get('TX.SN')
                    else:
                        summary['type'] = 'RX'
                        summary['tx_sn'] = None
                    # These items were in schedule_info for TX
                    summary['txfreq'] = metadata.get('Tx.Freq')
                    summary['txduty'] = metadata.get('Tx.Duty')
                    # These items were in schedule_info for RX
                    summary['rx_station'] = metadata.get('RX.STN')
                    summary['job_number'] = metadata.get('JOB.NUMBER')
                    summary['channel_station'] = metadata.get('CH.STN')
                    summary['component'] = metadata.get('CH.CMP')
                    summary['box_number'] = metadata.get('Box number')
                    summary['card_number'] = metadata.get('ChannelSerial')
                    # summary['schedule_info'] = schedule_info
                    if 'CH.OFFSET.XYZ1' in metadata:
                        offset1 = metadata.get('CH.OFFSET.XYZ1').split(':')
                    else:
                        offset1 = metadata.get('CH.XYZ1').split(':')
                    dx1_keys = ['dx1','dy1','dz1']
                    dx2_keys = ['dx2','dy2','dz2']
                    o1_floats = [float(o1) for o1 in offset1]
                    for idx in range(len(o1_floats)):
                        summary[dx1_keys[idx]] = o1_floats[idx]
                    if metadata.get('CH.CMP')[0] == 'H':
                        # mag field
                        if 'CH.ANTSN' in metadata:
                            summary['antenna_number'] = metadata.get('CH.ANTSN')
                        else:
                            summary['antenna_number'] = metadata.get('CH.NUMBER')
                        summary['azimuth'] = metadata.get('CH.AZIMUTH')
                        summary['inclination'] = metadata.get('CH.INCL')
                        summary['dx2'] = None
                        summary['dy2'] = None
                        summary['dz2'] = None
                    elif metadata.get('CH.CMP')[0] == 'E':
                        # electric field
                        summary['antenna_number'] = None
                        summary['azimuth'] = None
                        summary['inclination'] = None
                        if 'CH.OFFSET.XYZ2' in metadata:
                            offset2 = metadata.get('CH.OFFSET.XYZ2').split(':')
                        else:
                            offset2 = metadata.get('CH.XYZ2').split(':')
                        o2_floats = [float(o2) for o2 in offset2]
                        for idx in range(len(o2_floats)):
                            summary[dx2_keys[idx]] = o2_floats[idx]
                    elif metadata.get('CH.CMP')[0] == 'R':
                        # electric field, transmitter ref channel
                        summary['antenna_number'] = None
                        summary['azimuth'] = None
                        summary['inclination'] = None
                        if 'CH.OFFSET.XYZ2' in metadata:
                            offset2 = metadata.get('CH.OFFSET.XYZ2').split(':')
                        else:
                            offset2 = metadata.get('CH.XYZ2').split(':')
                        o2_floats = [float(o2) for o2 in offset2]
                        for idx in range(len(o2_floats)):
                            summary[dx2_keys[idx]] = o2_floats[idx]
                    else:
                        raise ValueError('CH.CMP not recognized as E, H, or Ref component')
                    if file_info['num_records']==0:
                        summary['start']=-1
                        summary['end']=-1
                        summary['valid']=False
                    else:
                        [dt_start,dt_end] = timeio.get_start_and_end_times_mountain(file_info)
                        summary['start'] = dt_start
                        summary['end'] = dt_end
                except Exception as e:
                    summary['valid'] = False
                    summary['error'] = repr(e)
                summaries.append(summary)
    df = pd.DataFrame(summaries)
    return df


def check_z3d_times(initialdir=initialdir,ask_dir=True):
    '''
    Get select metadata about each z3d file in a directory and its subdirectories
    '''
    if ask_dir:
        rootdir = util.getdir(initialdir=initialdir,title='Select Z3D directory',mustexist=True)
    else:
        rootdir = initialdir
    # check all z3d files in root directory structure
    pattern = "*.z3d"
    summaries = []
    for path, subdirs, files in os.walk(rootdir):
        for name in files:
            if fnmatch(name.lower(), pattern):
                filepath = os.path.join(path, name)
                summary = {}
                summary['valid']=True
                summary['gps_check']=False
                summary['error']=''
                summary['fullpath']=os.path.join(path,name)
                try:
                    gps_times = z3dio.read_gps_times(filename)
                    record_timedeltas = np.diff(gps_times)
                    if all(record_timedeltas == 1):
                        summary['gps_check'] = True
                    else:
                        summary['gps_check'] = False
                except Exception as e:
                    summary['valid'] = False
                    summary['error'] = repr(e)
                summaries.append(summary)
    df = pd.DataFrame(summaries)
    return df


def get_nearest_station(rx,location_file='RXlocations.csv',search_radius=200,campaign_increment=0):
    '''
    Populate nearest_station column in rx
    '''
    # get station number from location
    column_names = ['station_number','easting','northing']
    rx_true_loc = pd.read_csv(location_file,header=None,names=column_names,index_col='station_number')
    # (num_rows,num_columns) = rx.shape
    # nearest_station=np.empty(num_rows,np.nan)
    # rx.insert(num_columns,'nearest_station',nearest_station)
    rx['nearest_station']=np.nan
    for true_station,true_loc in rx_true_loc.iterrows():
        distance = np.sqrt((rx.easting-true_loc.easting)**2 + (rx.northing-true_loc.northing)**2)
        is_close = distance < search_radius
        unassigned = np.isnan(rx.nearest_station[is_close])
        if any(~unassigned):
            assigned_station = rx.nearest_station[is_close][~unassigned].iloc[0]
            raise ValueError('Assigning station '+str(true_station+campaign_increment)+', already assigned as '+str(assigned_station))
        rx.loc[is_close,'nearest_station']=true_station+campaign_increment


def find_overlaps(tx,rx,overlaps_csv,use_csv=True):
    '''
    Find TX-RX time series overlaps

    tx: tx metadata dataframe
    rx: rx metadata dataframe
    overlaps: overlaps dataframe, as returned by find_overlaps,
        with indeces corresponding to tx and rx

    Returns a dataframe with overlaps and tx/rx indeces
    '''
    # from dateutil.tz import tzutc
    # from dateutil.parser import parse
    # import Plot_txrx
    # def date_utc(s):
    #     return parse(s,tzinfos=tzutc)

    if use_csv:
        try:
            overlaps = pd.read_csv(overlaps_csv,index_col=0,parse_dates=['start','end'])
        except:
            overlaps = _find_overlaps(tx,rx)
            overlaps.to_csv(overlaps_csv)
    else:
        overlaps = _find_overlaps(tx,rx)

    overlaps.start = overlaps.start.dt.tz_localize('UTC')
    overlaps.start = overlaps.start.dt.tz_convert('US/Mountain')
    overlaps.end = overlaps.end.dt.tz_localize('UTC')
    overlaps.end = overlaps.end.dt.tz_convert('US/Mountain')
    return overlaps


def _find_overlaps(tx,rx):
    '''
    Find tx-rx overlaps, without csv or timezone conversions
    '''
    tx_inds = np.array([],dtype=int)
    rx_inds = np.array([],dtype=int)
    all_starts = np.array([],dtype='datetime64')
    all_ends = np.array([],dtype='datetime64')
    for irx, rxrow in rx.iterrows():
        starts = tx.start.map(lambda dt: max(rxrow.start,dt))
        ends = tx.end.map(lambda dt: min(rxrow.end, dt))
        is_overlap = starts<ends
        tx_inds = np.append(tx_inds,tx.index[is_overlap].values)
        rx_inds = np.append(rx_inds,[irx]*sum(is_overlap))
        all_starts = np.append(all_starts,starts[is_overlap].values)
        all_ends = np.append(all_ends,ends[is_overlap].values)
        if irx == 0:
            print('Organizing tx/rx overlaps...')
        if irx%100 == 0:
            print('Row ID '+str(irx)+', '+str(len(rx))+' rows.')
    overlaps = pd.DataFrame(data={'tx_ind':tx_inds,
                                  'rx_ind':rx_inds,
                                  'start':all_starts,
                                  'end':all_ends})
    return overlaps


def full_overlaps(tx,rx,overlaps):
    '''
    Get tx-rx overlap metadata

    tx: tx metadata dataframe
    rx: rx metadata dataframe
    overlaps: overlaps dataframe, as returned by find_overlaps,
        with indeces corresponding to tx and rx

    Returns a dataframe with metadata for every overlap
    '''
    orx = rx.loc[overlaps.rx_ind].reset_index()
    otx = tx.loc[overlaps.tx_ind].reset_index()
    # make ok columns
    orx['ok'] = orx.valid & orx.original & orx.in_range
    otx['ok'] = otx.valid & otx.original & otx.in_range
    overlaps_all = pd.concat([orx,otx,overlaps],
                             keys=['rx','tx','overlap'],
                             axis=1,
                             sort=False)
    # flatten multiindex
    overlaps_all.columns=pd.Index([e[0]+'_'+e[1] for e in overlaps_all.columns.tolist()])
    return overlaps_all


def read_zen_cals(directory=cal_dir,ask_dir=True):
    '''
    Read Zen board calibration files (zen###.cal) in the directory specified
    '''
    if ask_dir:
        rootdir = util.getdir(initialdir=directory,title='Select calibration directory',mustexist=True)
    else:
        rootdir = directory
    # antcal_pattern = 'antenna.cal'
    zencal_pattern = 'zen*.cal'
    zencals = []
    zennums = []
    for path, subdirs, files in os.walk(rootdir):
        for name in files:
            # if fnmatch(name,antcal_pattern):
            #     filepath = os.path.join(path,name)
            #     ant_cal = z3dio.read_antcal(filepath)
            # elif fnmatch(name,zencal_pattern):
            if fnmatch(name,zencal_pattern):
                filepath = os.path.join(path,name)
                # zencal = {}
                # zencal['filename']=name
                # zencal['cal']= z3dio.read_syscal(filepath)
                # zencal['cal_head'] = z3dio.read_syscal_header(filepath)
                # zen_no = zencal['GDP.SN']
                # zencals.append(zencal)

                syscal = z3dio.read_syscal(filepath)
                syscal_head = z3dio.read_syscal_header(filepath)
                if 'GDP.SN' in syscal_head:
                    syscal_zen_num = syscal_head['GDP.SN']
                else:
                    syscal_zen_num = name[3:-4]

                zencals.append(syscal)
                zennums.append(syscal_zen_num)
        # don't walk beyond top directory
        break
    cals = pd.concat(zencals,keys=zennums)
    # return (cals,antcal)
    return cals

