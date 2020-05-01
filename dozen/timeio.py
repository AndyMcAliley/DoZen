import datetime
from datetime import timezone

import numpy as np
import pandas as pd

datetimeformat = "%Y-%m-%d %H:%M:%S"
mst = timezone(datetime.timedelta(hours=-7),'MST')
mdt = timezone(datetime.timedelta(hours=-6),'MDT')
time_types = {'timestamp','datetime','datetime64'}

def year_1_days_to_dt(days):
    '''
    Convert days since year 1 (x axis on mpl plots) to datetime
    '''
    days_1970 = 719163
    dt = datetime.datetime.utcfromtimestamp(0)+datetime.timedelta(days=days-days_1970)
    return dt

def week_seconds_to_utc(gpsweek,gpsseconds,leapseconds=(19-37),astype='datetime'):
    '''Convert gps weeks/seconds to utc
    https://gist.github.com/jeremiahajohnson/eca97484db88bcf6b124#file-gpstoutc-py
    '''
    assert astype in time_types, 'astype unrecognized'
    if astype == 'datetime64':
        # timezone-agnostic numpy datetime64 (more precise)
        epoch = np.datetime64("1980-01-06 00:00:00")
        elapsed = np.timedelta64(int(gpsweek*7),'D') + np.timedelta64(int((gpsseconds+leapseconds)*1e9),'ns')
        dt = (epoch + elapsed)
        return dt
    if astype == 'datetime':
        # timezone-aware datetime (less precise)
        epoch = datetime.datetime.strptime("1980-01-06 00:00:00",datetimeformat)
        elapsed = datetime.timedelta(days=(gpsweek*7),seconds=(gpsseconds+leapseconds))
        # return datetime.datetime.strftime(epoch + elapsed,datetimeformat)
        dt = (epoch + elapsed).replace(tzinfo=timezone.utc)
        # to convert to timezone: dt.astimezone(mst)
        return dt
    if astype == 'timestamp':
        epoch = pd.to_datetime("1980-01-06 00:00:00")
        elapsed = pd.to_timedelta(int(gpsweek*7),'D') + pd.to_timedelta(int((gpsseconds+leapseconds)*1e9),'ns')
        dt = (epoch + elapsed).tz_localize('utc')
        return dt

def gps_timestamp_to_utc(timestamp,gpsweek,leapseconds=(19-37),astype='datetime'):
    assert astype in time_types, 'astype unrecognized'
    gps_seconds = timestamp/1024.
    return week_seconds_to_utc(gpsweek,gps_seconds,leapseconds,astype)

def aware_dt_to_mountain(dt):
    ''' Convert an aware datetime object to mountain time '''
    dst_dts = [datetime.datetime(2017,3,12,2,0,0,0,tzinfo=mst),
           datetime.datetime(2017,11,5,2,0,0,0,tzinfo=mdt),
           datetime.datetime(2018,3,11,2,0,0,0,tzinfo=mst),
           datetime.datetime(2018,11,4,2,0,0,0,tzinfo=mdt),
           datetime.datetime(2019,3,10,2,0,0,0,tzinfo=mst),
           datetime.datetime(2019,11,3,2,0,0,0,tzinfo=mdt),
           datetime.datetime(2020,3,8,2,0,0,0,tzinfo=mst),
           datetime.datetime(2020,11,1,2,0,0,0,tzinfo=mdt)
          ]
    i=0
    found = False
    for i in range(0,len(dst_dts)):
        if dt > dst_dts[i] and dt <= dst_dts[i+1]:
            found=True
            break
    if found:
        if i%2==0:
            return dt.astimezone(mdt)
        else:
            return dt.astimezone(mst)
    else:
        raise ValueError('time is outside (admittedly narrow) range of allowable timestamps')

# def week_rollovers(gps_time):
    ''' Return a list of indices where '''

def gps_week_seconds_to_mountain(gps_week_start,gps_seconds,leapseconds=(19-37),astype='datetime'):
    '''
    Convert an array of gps_seconds to mountain datetimes,
    given gps week at start of array
    '''
    #TODO: allow for numpy datetime64 and pandas timestamp
    assert astype in time_types, 'astype unrecognized'
    num_seconds = len(gps_seconds)
    # First, find all week rollovers
    # EITHER
    # rollovers = np.array([False]*len(gps_seconds))
    # rollovers[1:] = gps_seconds[:-1] > gps_seconds[1:]
    # OR use indices not booleans
    # robustify by looking for a rollover where the seconds roll back by at least most of a week
    # 600000 seconds is one week minus 80 minutes
    rollovers = np.where(gps_seconds[:-1] > gps_seconds[1:]+600000)[0]
    rollovers = rollovers+1
    # Create vector of gps_week
    gps_week = np.array([gps_week_start]*num_seconds)
    for rollover in rollovers:
        # increment all gps_week values after rollover
        gps_week[rollover:] = gps_week[rollover:]+1
    # Next, convert, as vectorized as possible
    if astype == 'datetime':
        epoch = np.datetime64("1980-01-06 00:00:00")
        elapsed_weeks = [np.timedelta64(int(gw*7),'D') for gw in gps_week]
        elapsed_seconds = [np.timedelta64(int(gs)+leapseconds,'s') for gs in gps_seconds]
        # return epoch + elapsed_weeks+elapsed_seconds
        dt = (epoch + elapsed_weeks + elapsed_seconds).astype(datetime.datetime) 
        dt = [d.replace(tzinfo=timezone.utc).astimezone(mst) for d in dt]
        dt = [aware_dt_to_mountain(d) for d in dt]
    elif astype == 'datetime64':
        # Not sure how to handle timezone conversion gracefully here since numpy implements no tz tools
        raise ValueError('datetime64 not implemented')
        epoch = np.datetime64("1980-01-06 00:00:00")
        elapsed_weeks = [np.timedelta64(int(gw*7),'D') for gw in gps_week]
        elapsed_seconds = [np.timedelta64(int(gs)+leapseconds,'s') for gs in gps_seconds]
        # return epoch + elapsed_weeks+elapsed_seconds
        dt = (epoch + elapsed_weeks + elapsed_seconds)
        # convert timezone
        # dt = [d.replace(tzinfo=timezone.utc).astimezone(mst) for d in dt]
        # dt = [aware_dt_to_mountain(d) for d in dt]
    elif astype == 'timestamp':
        epoch = pd.Timestamp("1980-01-06 00:00:00").tz_localize('utc')
        elapsed_weeks = [pd.Timedelta(int(gw*7),'D') for gw in gps_week]
        elapsed_seconds = [pd.Timedelta(gs+leapseconds,'s') for gs in gps_seconds]
        # return epoch + elapsed_weeks+elapsed_seconds
        dt = epoch + pd.Series(elapsed_weeks) + pd.Series(elapsed_seconds)
        dt = dt.dt.tz_convert('US/Mountain')
    return dt



def get_start_and_end_times_mountain(file_info,include_final_second=False,astype='datetime'):
    '''
    Get beginning and end times of a z3d file from the metadata 
    as extracted by z3dio.get_file_info
    if include_final_second, give end (actually, start?) time for last sample in last record
    '''
    # TODO: remove datetime usages from DoZen,
    #       convert all code to astype='timestamp'
    #       change default to 'timestamp'
    assert astype in time_types, 'astype unrecognized'
    start_gps_week = file_info['metadata']['GpsWeek']
    start_gps_timestamp = file_info['first_record']['GPS_timestamp']
    end_gps_timestamp = file_info['last_record']['GPS_timestamp']
    # check for GPS week rollover
    if file_info['GPS_week_rollover']:
        end_gps_week = start_gps_week + 1
    else:
        end_gps_week = start_gps_week
    dt_start_utc = gps_timestamp_to_utc(start_gps_timestamp,start_gps_week,astype=astype)
    if include_final_second:
        final_second_fraction = (file_info['metadata']['A/D Rate']-1)/file_info['metadata']['A/D Rate']
        dt_end_utc = gps_timestamp_to_utc(end_gps_timestamp+final_second_fraction*1024,end_gps_week,astype=astype)
    else:
        dt_end_utc = gps_timestamp_to_utc(end_gps_timestamp,end_gps_week,astype=astype)
    if astype=='datetime64':
        # convert to timestamp since it is time zone aware
        dt_start_utc_ts = pd.Timestamp(dt_start_utc).tz_localize('utc')
        dt_start_mountain_ts = dt_start_utc_ts.tz_convert('US/Mountain')
        dt_end_utc_ts = pd.Timestamp(dt_end_utc).tz_localize('utc')
        dt_end_mountain_ts = dt_end_utc_ts.tz_convert('US/Mountain')
        # remove time zone, convert to datetime64
        dt_start_mountain = dt_start_mountain_ts.tz_localize(None).to_datetime64()
        dt_end_mountain = dt_end_mountain_ts.tz_localize(None).to_datetime64()
    if astype=='datetime':
        dt_start_mountain = aware_dt_to_mountain(dt_start_utc)
        dt_end_mountain = aware_dt_to_mountain(dt_end_utc)
    if astype=='timestamp':
        dt_start_mountain = aware_dt_to_mountain(dt_start_utc)
        dt_end_mountain = aware_dt_to_mountain(dt_end_utc)
    return [dt_start_mountain,dt_end_mountain]


def get_times_mountain(z3d,time_range=[None,None]):
    '''
    Return an array of datetimes for all data in z3d file,
    using metadata as extracted by z3dio.read_z3d
    '''
    [dt_start,dt_end] = get_start_and_end_times_mountain(z3d, 
                                                         include_final_second=True, 
                                                         astype='timestamp')
    dt_times = pd.date_range(start=dt_start, end=dt_end,
                             periods=len(z3d['data']), tz=None,
                             normalize=False, name=None, closed=None)
    if time_range[0]:
        dt_start = time_range[0]
    if time_range[1]:
        dt_end = time_range[1]
    window_slice = slice_to_window_series(dt_times,dt_start,dt_end)
    # curve_data = z3d['data'][window_slice]
    dt_times = dt_times[window_slice]
    dt_times = dt_times.tz_localize(None).astype('datetime64[ns]')
    # hv_time_range = (dt_start.tz_localize(None), dt_end.tz_localize(None))
    # dt_start = np.datetime64(dt_start.tz_localize(None))
    # dt_end = np.datetime64(dt_end.tz_localize(None))
    return dt_times


def slice_to_window_series(x,start,stop):
    '''
    Return a slice of x between start and stop
    x is a pandas date_range
    Assumes x is ordered
    x, start, and stop must have the same time-zone awareness
    '''

    first = x.searchsorted(start)
    last = x.searchsorted(stop,side='right')
    return slice(first,last)
    # after_start = np.where([xe>=start for xe in x])[0]
    # if len(after_start)>0:
    #     first = after_start[0]
    #     after_last = np.where([xe>stop for xe in x])[0]
    #     if len(after_last)>0:
    #         last = after_last[0]
    #     else:
    #         last = None
    # else:
    #     first = 0
    #     last = 0
    # return slice(first,last)





