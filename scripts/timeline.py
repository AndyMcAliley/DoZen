'''
Interactive timeline, with time series plotting
'''

from functools import partial
from time import sleep
import os, sys
from pathlib import Path
from shutil import copyfile

import numpy as np
import pandas as pd
# plotting
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import panel as pn
import holoviews as hv
from holoviews.operation.datashader import datashade, dynspread
from holoviews.streams import Params
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models.tools import HoverTool, TapTool
from bokeh.models import CustomJS, ColumnDataSource
from bokeh.io import show
from bokeh.plotting import figure
import param

# add parent directory to path
# import os, sys
# module_path = os.path.abspath(os.path.join('..'))
# if module_path not in sys.path:
    # sys.path.append(module_path)

# import DoZen
# add parent directory to path
try:
    # works when called externally via panel serve
    script_path = os.path.abspath(os.path.join(__file__,'..'))
    module_path = os.path.abspath(os.path.join(__file__,'..','..'))
except NameError:
    # works in notebook
    script_path = str(Path().resolve())
    module_path = str(Path().resolve().parent)
if module_path not in sys.path:
    sys.path.append(module_path)
from dozen import z3dio, timeio, viz, z3d_directory, util

# hv.notebook_extension('bokeh','matplotlib')
hv.extension('bokeh','matplotlib')
pn.extension()

settings_file = 'timeline_settings.csv'
try:
    settings = pd.read_csv(settings_file,index_col=None)
except FileNotFoundError:
    print('Settings file not found locally. Copying default into local folder.')
    old_settings_file = settings_file
    settings_file = os.path.join(script_path,'timeline_settings.csv')
    copyfile(settings_file,old_settings_file)
    settings = pd.read_csv(settings_file,index_col=None)

print('Settings loaded from '+settings_file)

def rx_select(attr,old,new,rxsource,txsources):
    
    #selection_text.object = 'rx_select'
    #print('rx_select method') 
    #print(type(rxsource))
    #print(type(rxsource.data))
    #print(rxsource.data.keys())
    #print("Those are rxsource.data's keys")
    
    rxfile=[]
    rxpath=[]
    txfile=[]
    txpath=[]
    
    if len(new)>0:
        rxpath = rxsource.data['fullpath'][new]
        rxfile = rxsource.data['filename'][new]
        
        for txsource in txsources:
            txsource_indices = txsource.selected.indices
            if txsource_indices:
                #print('Found a TX while in rx_select!')
                txpath = txsource.data['fullpath'][txsource_indices]
                txfile = txsource.data['filename'][txsource_indices]
        if rxpath and not txpath:
            #print(rxpath[0])
            selection_text.object = 'RX file: '+rxfile[0]
            #ts.filenames = [txpath,rxpath]
        else:
            if rxpath and txpath:
                selection_text.object = 'TX file: '+txfile[0]+'<br /> RX file: '+rxfile[0]
            #print("Don't plot yet! TX is also selected.")
    else:
        selection_text.object = 'Nothing selected'
        #print('len of new is <= 0')
        #print(new)
    tl.rxfile = rxfile
    tl.rxpath = rxpath
    tl.txfile = txfile
    tl.txpath = txpath


def tx_select(attr,old,new,txsource,rxsource):
    #print('tx_select method')
    rxfile=[]
    rxpath=[]
    txfile=[]
    txpath=[]

    if len(new)>0:
        #print('new tx_select')
        txpath = txsource.data['fullpath'][new]
        txfile = txsource.data['filename'][new]
        rxpath = txsource.data['rxfullpath'][new]
        rxfile = txsource.data['rxfilename'][new]

        selection_text.object = 'TX file: '+txfile[0]+'<br /> RX file: '+rxfile[0]
        #ts.filenames = [txpath,rxpath]
    else:
        rxsource_indices = rxsource.selected.indices
        if len(rxsource_indices)>0:
            #print('rx tx_select')
            rxpath = rxsource.data['fullpath'][rxsource_indices]
            rxfile = rxsource.data['filename'][rxsource_indices]
            #ts.filenames = [txpath,rxpath]
            selection_text.object = 'RX file: '+rxfile[0]
            #print(rxpath)
        else:
            #print('nothing tx_select')
            #ts.filenames = [txpath,rxpath]
            selection_text.object = 'Nothing selected'
    tl.rxfile = rxfile
    tl.rxpath = rxpath
    tl.txfile = txfile
    tl.txpath = txpath


def plot_selected(events):
    #print(type(events))
    #print(tl.txfile)
    #print(tl.rxfile)
    #if tl.txfile and tl.rxfile:
    #    selection_text.object = 'TX file: '+tl.txfile[0]+'<br /> RX file: '+tl.rxfile[0]
    #elif tl.rxfile:
    #    selection_text.object = 'RX file: '+rxfile[0]
    #else:
    #    selection_text.object = 'Nothing selected'
    #print('Changing ts.filenames in 3... 2... 1...')
    ts.filenames = [tl.txpath,tl.rxpath]
    #print('Calling ts.time_series in 3...2...1...')
    #ts.time_series()
    #print('End of plot_selected')


class Timeline(param.Parameterized):
    campaign = param.ObjectSelector(default=settings.name.values[0],objects=settings.name.values.tolist())
    rx_component = param.ObjectSelector(default='Ex', objects=['Ex','Ey','Hx','Hy','Hz'])
    
    rxfile = []
    rxpath = []
    txfile = []
    txpath = []

    @param.depends('rx_component','campaign', watch=True)
    def timeline(self):
        #print('tl.timeline triggered!')
        [self.tx_source,self.rx_source] = self.sources[self.campaign]
        p = self.plot_timeline(self.rx_component)
        return p
        
    
    def __init__(self):
        self.sources={}
        for campaign_name in settings.name.values:
            # Check for duplicate names
            assert sum(settings.name==campaign_name)==1, 'Unique entry for '+campaign_name+' not found'
            self.sources[campaign_name]=self.load_campaign(campaign_name)  
        [self.tx_source,self.rx_source]=self.sources[self.campaign]
    
    
    def load_campaign(self,campaign_name):

        cs = settings.loc[settings.name==campaign_name]

        print('Loading '+campaign_name+' from '+cs.rx_dir.values[0])

        rx = z3d_directory.get_z3d_directory_info(initialdir=cs.rx_dir.values[0],ask_dir=False)
        tx = z3d_directory.get_z3d_directory_info(initialdir=cs.tx_dir.values[0],ask_dir=False)


        # figure out station number
        # look for 3 or 4 digit number plus letters
        parsed_station_folder = rx.filepath.str.extract(r'/(?:BC|sub|Sub|SUB|LMA)(\d+[a-zA-Z]+)/')[0]

        # get station number from location
        z3d_directory.get_nearest_station(rx,
                                          location_file=cs.rx_station_location_file.values[0],
                                          search_radius=cs.search_radius.values[0],
                                          campaign_increment=cs.campaign_increment.values[0])

        # remove duplicates, keeping only the first
        tx.drop_duplicates(subset='filename',keep='first',inplace=True)
        rx.drop_duplicates(subset='filename',keep='first',inplace=True)
        # drop tx Ex files
        tx.drop(tx[tx.type=='RX'].index,inplace=True)
        # drop invalid files
        rx.drop(rx[~rx.valid].index,inplace=True)
        tx.drop(tx[~tx.valid].index,inplace=True)
        # drop unassigned stations
        rx.dropna(subset=['nearest_station'],inplace=True)

        # overlaps = z3d_directory.find_overlaps(tx,rx,overlaps_csv=os.path.join(script_path,cs.overlaps_csv.values[0]))
        overlaps = z3d_directory.find_overlaps(tx,rx,overlaps_csv=cs.overlaps_csv.values[0])

        rxstn_overlaps = rx.loc[overlaps.rx_ind].nearest_station.values
        unique_rx_stations,station_indices = np.unique(rx.nearest_station.values,return_inverse=True)
        try:
            rxstn_overlaps_indices = np.asarray([np.asarray(unique_rx_stations==r).nonzero()[0] for r in rxstn_overlaps])
        except:
            print('Station mismatch. Try deleting ' + 
                  cs.overlaps_csv.values[0] + 
                  ' and run again. If that fails, check ' +
                  cs.rx_station_location_file.values[0] +
                  ' and the timeline settings file.')
            raise

        rx_source = rx.copy()
        rx_source['station'] = rx_source['nearest_station']
        rx_source['station_index'] = station_indices
        rx_source['txfreq']=0
        rx_source.start = rx_source.start.map(lambda dt: dt.replace(tzinfo=None))
        rx_source.end = rx_source.end.map(lambda dt: dt.replace(tzinfo=None))

        tx_source = overlaps.copy()
        tx_source['station'] = tx.loc[tx_source.tx_ind].rx_station.values
        tx_source['type'] = tx.loc[tx_source.tx_ind].type.values
        tx_source['filename'] = tx.loc[tx_source.tx_ind].filename.values
        tx_source['fullpath'] = tx.loc[tx_source.tx_ind].fullpath.values
        tx_source['sampling_rate'] = tx.loc[tx_source.tx_ind].sampling_rate.values
        tx_source['txfreq'] = tx.loc[tx_source.tx_ind].txfreq.values
        tx_source['component'] = rx.loc[tx_source.rx_ind].component.values
        tx_source['rxfullpath'] = rx.loc[tx_source.rx_ind].fullpath.values
        tx_source['rxfilename'] = rx.loc[tx_source.rx_ind].filename.values
        tx_source['station_index'] = rxstn_overlaps_indices
        tx_source.start = tx_source.start.map(lambda dt: dt.replace(tzinfo=None))
        tx_source.end = tx_source.end.map(lambda dt: dt.replace(tzinfo=None))

        return(tx_source,rx_source)
    
    
    def plot_timeline(self,rx_component):
        # plot with bokeh
        unique_rx_stations,_ = np.unique(self.rx_source.station.values,return_inverse=True)
        tx_stations = self.tx_source.station.unique()
        tx_stations.sort()
        tx_colors = ['#FFCCCC','#FF00FF','#DDDD00','#FF8800','#00FFFF','#FF0000','#0000FF','#00FF00']
        tx_colors = tx_colors[:len(tx_stations)]
        TOOLS = "save,pan,xwheel_zoom,box_zoom,tap,reset"
        bar_height = 0.5

        min_date = self.rx_source.start.min()
        max_date = self.rx_source.end.max()
        date_range = max_date-min_date
        xmin = min_date - 0.05*date_range
        xmax = max_date + 0.05*date_range

        p = figure(x_axis_type="datetime", tools=TOOLS, 
                   plot_height=300, plot_width=800, 
                   title = 'Active receivers',
                   x_range=(xmin,xmax),
                   y_range=(-bar_height,len(unique_rx_stations)-bar_height),
                   active_scroll="xwheel_zoom"
                  )
        p.add_tools(HoverTool(tooltips = [('station', '@type @station'),
                                          ('start', '@start{%F %T}'),
                                          ('end', '@end{%F %T}'),
                                          ('TX frequency', '@txfreq{0.[00000]} Hz'),
                                          ('sampling rate', '@sampling_rate Hz'),
                                          ('filename', '@filename')],
                              formatters={'start': 'datetime',
                                          'end': 'datetime'}))

        #rx_component = self.rx_component
        rcsource = ColumnDataSource(self.rx_source[self.rx_source.component==rx_component])
        tcsource = self.tx_source[self.tx_source.component==rx_component]
        
        p.grid.grid_line_alpha=0.3
        p.hbar('station_index',bar_height,'start','end',source=rcsource,
               fill_color="#DDDDDD", line_color="black")
        # plot tx
        tcsources=[]
        for i_station,tx_station,tx_color in zip(np.arange(len(tx_stations)),tx_stations,tx_colors):
            tcsources.append(ColumnDataSource(tcsource[tcsource.station==tx_station]))
            tcsources[i_station].selected.on_change('indices',
                                                    partial(tx_select,
                                                            txsource=tcsources[i_station],
                                                            rxsource=rcsource))

            p.hbar('station_index',bar_height,'start','end',source=tcsources[i_station],
                   fill_color=tx_color,line_color=tx_color,legend='{0:.0f}'.format(tx_station))

        # get rx to respond on click
        rcsource.selected.on_change('indices',
                                    partial(rx_select,
                                            rxsource=rcsource,
                                            txsources=tcsources))    
        
        p.xaxis.formatter = viz.verbose_formatter()
        y_tick_values = [v for v in range(len(unique_rx_stations))]
        y_tick_labels = ["%d" % d for d in unique_rx_stations]
        y_labels = dict(zip(y_tick_values,y_tick_labels))
        p.yaxis.ticker = y_tick_values
        p.yaxis.major_label_overrides = y_labels
        p.yaxis.axis_label = "Receiver Station"
        p.legend.location = "top_left"
        p.legend.click_policy="hide"
        return p


class Time_series(param.Parameterized):
    #sources = []
    filenames = param.List(default=[[],[]])
    #filenames = [[],[]]
    
    #@param.depends('filenames', watch=True)
    def time_series(self):
        empty_curve = hv.Curve([]).options(xaxis=None, yaxis=None, show_frame=True,
                                               height=200,width=800)
        print('Plotting time series')
        tx_path = self.filenames[0]
        rx_path = self.filenames[1]
        #rx_path = []
        #print(tx_path,rx_path)
        if tx_path and rx_path:
            #print('tx and rx!')
            # get plot time limits
            fi = z3dio.get_file_info(tx_path[0])
            [tstart,tend] = timeio.get_start_and_end_times_mountain(fi,include_final_second=True,astype='timestamp')
            #tstart=pd.Timestamp(tstart.replace(tzinfo=None))
            #tend=pd.Timestamp(tend.replace(tzinfo=None))

            d1 = viz.plot_ts(tx_path[0],time_range=[tstart,tend])
            d1rd = d1.redim(signal='TX Signal')
            #d2 = viz.plot_ts(rx_path[0])
            d2 = viz.plot_ts(rx_path[0],time_range=[tstart,tend])
            d2rd = d2.redim(signal='RX Signal')
            d1rd.opts(height=200,width=800)
            d2rd.opts(height=200,width=800)
            q=hv.Layout(d1rd+d2rd).cols(1)
            #print('I did it')
        elif rx_path:
            #print('only rx.')
            d2 = viz.plot_ts(rx_path[0])
            d2rd = d2.redim(signal='RX Signal')
            #print('redimensioned!')
            #print(type(d2rd))
            q = hv.Layout((empty_curve,d2rd.opts(height=200,width=800))).cols(1)
            #print('I did it')
        else:
            #print('I got nuthin.')
            q = hv.Layout((empty_curve,empty_curve)).cols(1)
        #print('Shootin back a q atchu')
        #print(type(q))
        return q


tl = Timeline()
ts = Time_series()
selection_text = pn.pane.Markdown('Nothing selected')
plot_button = pn.widgets.Button(name='Plot time series', button_type='primary')
plot_button.param.watch(plot_selected,'clicks')

# serve app
# pn.Column(tl.param,tl.timeline)
pn.Column(tl.param,tl.timeline,selection_text,plot_button,ts.time_series).servable()

# # save app as static html (no plotting)
# panel = pn.Column(tl.param,tl.timeline)
# panel.save('timeline.html',embed=True)

