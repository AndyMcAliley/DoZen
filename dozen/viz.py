'''
Read a z3d
Plot the time series
Use datashader to handle large time series data
'''

# from glob import glob
from datetime import timedelta

import numpy as np
import pandas as pd
import holoviews as hv
from holoviews import renderer,Curve
from holoviews.operation.datashader import datashade, shade, dynspread
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import tools as btools
from bokeh.models import Button
from bokeh.models.widgets import Select, PreText
from bokeh.layouts import layout
# import panel as pn

from .z3dio import read_z3d
from . import timeio
# import .process

# plot settings
# brenderer = renderer('bokeh')
hv.config.image_rtol = 0.1
dynspread.max_px=20
dynspread.threshold=0.5
shade.cmap="#30a2da" # to match HV Bokeh default
shade.cmap="#3062aa"
# shade.cmap=hv.plotting.util.process_cmap('Blues')

# limit doc
__all__ = ['verbose_formatter','no_logo','plot_ts']

def verbose_formatter():
    '''
    Format bokeh time axis verbosely
    '''
    vf = DatetimeTickFormatter()
    vf.microseconds = ['%f us']
    vf.milliseconds = ['%S.%3N s']
    # vf.milliseconds = ['%H:%M:%S.%3N']
    vf.seconds = ['%H:%M:%S']
    vf.minsec = ['%H:%M:%S']
    vf.minutes = ['%m/%d %H:%M']
    vf.hourmin = ['%m/%d %H:%M']
    vf.hours = ['%m/%d %H:%M']
    vf.days = ['%m/%d', '%a%d']
    vf.months = ['%m/%Y', '%b %Y']
    vf.years = ['%Y']
    return vf

def no_logo(plot, element):
    plot.state.toolbar.logo = None


def plot_ts(filename,time_range=[None,None]):
    '''
    Plot the time series in z3d file filename
    return a holoviews dynspread object
    time_range must be tz-aware
    '''
    # print('reading...')
    z3d = read_z3d(filename)
    sampling_rate = z3d['metadata']['A/D Rate']
    if 'TX.SENSE' in z3d['metadata']:
        # NOTE: this only ids channel 1 TX files
        station_type = 'TX'
    else:
        station_type = 'RX'
    # These items were in schedule_info for RX
    try:
        station = z3d['metadata']['CH.STN']
    except KeyError:
        try:
            station = z3d['metadata']['RX.STN']
        except KeyError:
            station = ''
    try:
        station = str(int(station))
    except ValueError:
        station = str(station)
    channel = z3d['metadata']['Channel']
    component = z3d['metadata']['CH.CMP']
    if z3d['metadata']['CH.CMP'][0] == 'H':
        # mag field
        # antenna_number = z3d['metadata']['CH.ANTSN']
        pass
    elif z3d['metadata']['CH.CMP'][0] == 'E':
        # electric field
        # antenna_number = None
        pass
    elif z3d['metadata']['CH.CMP'][0] == 'R':
        # electric field, transmitter ref channel
        # antenna_number = None
        pass
    else:
        raise ValueError('CH.CMP not recognized as either E or H component')
    if z3d['num_records']==0:
        valid = False
    else:
        [dt_start,dt_end] = timeio.get_start_and_end_times_mountain(z3d,include_final_second=True,astype='timestamp')

    # print('getting times...')
    # dt_start_naive = dt_start.replace(tzinfo=None)
    # dt_end_naive = dt_end.replace(tzinfo=None)
    # dt_start_naive = pd.Timestamp(dt_start.replace(tzinfo=None))
    # dt_end_naive = pd.Timestamp(dt_end.replace(tzinfo=None))
    # period = timedelta(seconds=1)/sampling_rate
    # sampling_increment = timedelta(seconds=1./sampling_rate)
    # dt_times = pd.date_range(start=dt_start_naive, periods=len(z3d['data']),
    #                          freq=sampling_increment, tz=None,
    #                          normalize=False, name=None, closed=None)
    dt_times = pd.date_range(start=dt_start, end=dt_end,
                          periods=len(z3d['data']), tz=None,
                          normalize=False, name=None, closed=None)#.values

    # print('windowing data...')
    # print(type(dt_start))
    # print(type(time_range[0]))
    # window data to time_range
    if time_range[0]:
        dt_start = time_range[0]
    if time_range[1]:
        dt_end = time_range[1]
    window_slice = timeio.slice_to_window_series(dt_times,dt_start,dt_end)
    curve_data = z3d['data'][window_slice]
    dt_times = dt_times[window_slice]
    # ... or don't
    # curve_data = z3d['data']
    # print('windowed!')

    # holoviews seems to want timezone-naive type datetime64[ns]. Throws warning otherwise
    try:
        dt_times = dt_times.tz_localize(None).astype('datetime64[ns]')
        # time_range[0] = time_range[0].tz_localize(None)
        # time_range[1] = time_range[1].tz_localize(None)
        # hv_time_range = (time_range[0].tz_localize(None), time_range[1].tz_localize(None))
        hv_time_range = (dt_start.tz_localize(None), dt_end.tz_localize(None))
        # print(type(dt_times))
        # print(type(dt_times[0]))
        # print(dt_times[0])
        # print(type(dt_start))
        # print(dt_start)
        dt_start = np.datetime64(dt_start.tz_localize(None))
        # print(type(dt_start))
        # print(dt_start)
        dt_end = np.datetime64(dt_end.tz_localize(None))
    except Exception as e:
        print(e)
        traceback.print_exc()

    # print(len(curve_data))
    # print(len(dt_times))
    # print(hv_time_range)

    # print('plotting...')
    title_str = '{} Station: {}, '+\
        'Channel {:.0f} ({}), {} to {} Mountain'
    title_str = title_str.format(station_type, station,
                                 channel, component,
                                 np.datetime_as_string(dt_start, unit='s'),
                                 np.datetime_as_string(dt_end, unit='s'))
    ts_tools = ['save','pan','xwheel_zoom','box_zoom','undo','reset']
    signal = hv.Dimension('signal',label='Signal',unit='V')
    time = hv.Dimension('time',label='Time',range=hv_time_range)
    curve = Curve((dt_times, curve_data),
                  time,
                  signal,
                  label=title_str
                 )#.redim.range(time=time_range)
                     # kdims=['Time'],vdims=['Signal, V']
    # dscurve = datashade(curve, cmap=["blue"]).opts(width=800)
    # dscurve = datashade(curve,name=title_str).opts(width=800)
    dscurve = dynspread(datashade(curve).opts(#framewise=True,
                                              xformatter=verbose_formatter(),
                                              default_tools=ts_tools))
    dsfinal=dscurve.opts(plot=dict(hooks=[no_logo]))
    # print('done!')
    return dsfinal


class QuickPlot():
    # plot some time series!
    start = 0
    end = None
    
    ts_tools = ['save','pan','xwheel_zoom','box_zoom','undo','reset']
    
    def ds(self,curve):
        # apply datashading and dynspread to a curve
        if self.end:
            dscurve = hv.operation.datashader.dynspread(hv.operation.datashader.datashade(curve).opts(
                xformatter=verbose_formatter(),
                default_tools=self.ts_tools))
        else:
            dscurve = hv.operation.datashader.dynspread(hv.operation.datashader.datashade(curve).opts(
                default_tools=self.ts_tools))
        dsfinal=dscurve.opts(plot=dict(hooks=[no_logo]))
        return dsfinal
    
    def times(self,ts):
        if self.end:
            dt_times = pd.date_range(start=self.start, end=self.end,
                          periods=len(ts), tz=None,
                          normalize=False, name=None, closed=None)#.values
            dt_times = dt_times.tz_localize(None).astype('datetime64[ns]')
            hv_time_range = (self.start.tz_localize(None), self.end.tz_localize(None))
            #dt_start = np.datetime64(dt_start.tz_localize(None))
            #dt_end = np.datetime64(dt_end.tz_localize(None))
            x_dim = hv.Dimension('time',label='Time',range=hv_time_range)
            
        else:
            dt_times = np.arange(len(ts))
            x_dim = hv.Dimension('sample',label='Sample number',range=(0,len(ts)))
        return (dt_times,x_dim)

    def p1(self,ts,x_arr=None):
        signal = hv.Dimension('signal',label='Signal')
        
        if x_arr is not None:
            # User-specified x axis
            dt_times = x_arr
            x_dim = hv.Dimension('x',label='x')
        # assume 4096 Hz
        elif self.end:
            dt_times = pd.date_range(start=self.start, end=self.end,
                          periods=len(ts), tz=None,
                          normalize=False, name=None, closed=None)#.values
            dt_times = dt_times.tz_localize(None).astype('datetime64[ns]')
            hv_time_range = (self.start.tz_localize(None), self.end.tz_localize(None))
            #dt_start = np.datetime64(dt_start.tz_localize(None))
            #dt_end = np.datetime64(dt_end.tz_localize(None))
            x_dim = hv.Dimension('time',label='Time',range=hv_time_range)
            
        else:
            dt_times = np.arange(len(ts))
            x_dim = hv.Dimension('sample',label='Sample number',range=(0,len(ts)))
        
        curve = hv.Curve((dt_times, ts),x_dim,signal)
        p = self.ds(curve).opts(height=250,width=800)
        return p
    
    def p2(self,tx_ts,rx_ts,xlim=None,rxlim=None,txlim=None):
        assert len(tx_ts)==len(rx_ts), 'tx and rx time series must be the same length'
        
        tx_signal = hv.Dimension('tx_signal',label='TX Signal',unit='A')
        rx_signal = hv.Dimension('rx_signal',label='RX Signal',unit='V')
        
        # assume 4096 Hz
        if self.end:
            dt_times = pd.date_range(start=self.start, end=self.end,
                          periods=len(tx_ts), tz=None,
                          normalize=False, name=None, closed=None)#.values
            dt_times = dt_times.tz_localize(None).astype('datetime64[ns]')
            hv_time_range = (self.start.tz_localize(None), self.end.tz_localize(None))
            #dt_start = np.datetime64(dt_start.tz_localize(None))
            #dt_end = np.datetime64(dt_end.tz_localize(None))
            x_dim = hv.Dimension('time',label='Time',range=hv_time_range)
            
        else:
            dt_times = np.arange(len(tx_ts))
            x_dim = hv.Dimension('sample',label='Sample number')
        
        tx_curve = hv.Curve((dt_times, tx_ts),x_dim,tx_signal)
        rx_curve = hv.Curve((dt_times, rx_ts),x_dim,rx_signal)
        txp = self.ds(tx_curve).opts(height=200,width=800)
        rxp = self.ds(rx_curve).opts(height=200,width=800)
        if xlim:
            txp = txp.opts(xlim=xlim)
            rxp = rxp.opts(xlim=xlim)
        if rxlim:
            rxp = rxp.opts(ylim=rxlim)
        if txlim:
            txp = txp.opts(ylim=txlim)
        return (txp+rxp).cols(1)


def _old_plot_series(ts,ts2=None,sampling_rate=4096):
    '''
    Time series plots
    Plot one or two time series on the same x-axis
    '''
    fig,ax1=plt.subplots(figsize=(9,6))
    c1 = 'tab:blue'
    ax1.set_xlabel('seconds')
    time_in_ms = np.arange(len(ts))/float(sampling_rate)
    ax1.plot(time_in_ms,ts,color=c1)
    ax1.tick_params(axis='y',labelcolor=c1)
    ax1.set_ylabel('Signal 1',color=c1)
    if ts2 is not None:
        c2 = 'tab:red'
        ax2 = ax1.twinx()
        ax2.plot(time_in_ms,ts2,color=c2)
        ax2.set_ylabel('Signal 2',color=c2)
        ax2.tick_params(axis='y',labelcolor=c2)
    fig.tight_layout()
    # plt.show()
    if ts2 is not None:
        return(fig,ax1,ax2)
    else:
        return(fig,ax1)



def _old_plot_ts(txfilename,rxfilename,rx_station):
    tx=z3dio.read_z3d(txfilename)
    rx=z3dio.read_z3d(rxfilename)

    sampling_rate = round(rx['metadata']['A/D Rate'])
    tx_mountain_times = timeio.gps_week_seconds_to_mountain(tx['metadata']['GpsWeek'],tx['gps_times'])
    rx_mountain_times = timeio.gps_week_seconds_to_mountain(rx['metadata']['GpsWeek'],rx['gps_times'])
    [tx_start,tx_end] = timeio.get_start_and_end_times_mountain(tx)
    [rx_start,rx_end] = timeio.get_start_and_end_times_mountain(rx)

    rx_overlap_start_second = rx_mountain_times.index(tx_start)
    rx_start_index = rx_overlap_start_second*sampling_rate
    rx_end_index = rx_start_index + len(tx['data'])

    tx_start_naive = tx_start.replace(tzinfo=None)

    # form array of datetimes without pandas
    # TODO: stop assuming total overlap of tx and rx
    # TODO: stop assuming same sampling rate
    period = timedelta(seconds=1)/sampling_rate
    tx_times = tx_start_naive+np.arange(len(tx['data']))*period
    # plot with matplotlib
    fig,ax1=plt.subplots(figsize=(9,6))
    mpl_times = matplotlib.dates.date2num(tx_times)
    c1 = 'tab:blue'
    c2 = 'tab:red'
    ax1.set_xlabel('Date/Time, Mountain Time')
    ax1.plot_date(mpl_times,tx['data'],fmt='-',color=c1)
    ax1.tick_params(axis='y',labelcolor=c1)
    if 'TX.SENSE' in rx['metadata']:
        ax1.set_ylabel('Transmitter Current, A',color=c1)
    else:
        ax1.set_ylabel('Transmitter Signal, V',color=c1)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Receiver Signal, V',color=c2)
    ax2.plot_date(mpl_times,rx['data'][rx_start_index:rx_end_index],fmt='-',color=c2)
    ax2.tick_params(axis='y',labelcolor=c2)
    # datetimeformat = "%Y-%m-%d %H:%M:%S"


    # show it!
    # plt.title(ntpath.basename(filename))
    title_string = 'TX: '+str(int(tx['metadata']['RX.STN']))+', RX: '+\
                   str(int(rx_station))+', '+str(tx['metadata']['Tx.Freq'])+\
                   ' Hz, '+rx['metadata']['CH.CMP']
    plt.title(title_string)
    fig.tight_layout()
    plt.show()

