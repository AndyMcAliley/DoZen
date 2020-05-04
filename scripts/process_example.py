import os
import sys
import traceback

import numpy as np
from scipy import signal
import pandas as pd

from dozen import z3d_directory, z3dio, timeio, process

# Inputs and settings

# survey campaign
rx_file = '../data/preprocess/campaign_rx.csv'
tx_file = '../data/preprocess/campaign_tx.csv'
overlaps_file = '../data/preprocess/overlaps.csv'
# calibration files
cal_dir = '../data/calibrations/'
antcal_file = '../data/calibrations/antenna.cal'
# Results file
results_file = 'DoZen.avg'
# Save odd harmonic fourier coefficients for every window for every time series?
save_coefficients = True
coefficients_dir = 'coeff'
os.mkdir(coefficients_dir)

# mag orientations

# subset data offsets
start_offset = 2
end_offset = 0

# filtering settings
pad_before = 34
pad_after = 50
notch_bandwidth = 0.2

# short time Fourier transform settings
# window_shape = ('kaiser',3*np.pi)
# window_shape = ('tukey',1./3.)
# window_shape = 'tukey_ten_sample_taper'
window_shape = ('hann')
# window_shape = 'boxcar'
cycles_per_window = 2 # cycle defined as one transmitter waveform
strides_per_cycle = 0.25

# Read z3d directory info
print('Reading directory info...')
rx = pd.read_csv(rx_file,index_col=0)
tx = pd.read_csv(tx_file,index_col=0)
# get start and end as dates
for i_row in rx.index:
    file_info = z3dio.get_file_info(rx.loc[i_row,'fullpath'])
    if file_info['num_records']==0:
        rx.at[i_row,'start'] = -1
        rx.at[i_row,'end']= -1
        rx.at[i_row,'valid']=False
    else:
        [start,end] = timeio.get_start_and_end_times_mountain(file_info)
        rx.at[i_row,'start'] = pd.Timestamp(start).tz_convert('US/Mountain')
        rx.at[i_row,'end']= pd.Timestamp(end).tz_convert('US/Mountain')
for i_row in tx.index:
    file_info = z3dio.get_file_info(tx.loc[i_row,'fullpath'])
    if file_info['num_records']==0:
        tx.at[i_row,'start'] = -1
        tx.at[i_row,'end']= -1
        tx.at[i_row,'valid']=False
    else:
        [start,end] = timeio.get_start_and_end_times_mountain(file_info)
        tx.at[i_row,'start'] = pd.Timestamp(start).tz_convert('US/Mountain')
        tx.at[i_row,'end']= pd.Timestamp(end).tz_convert('US/Mountain')

# Fix errors in station numbering, invalid/duplicate files
# check for duplicate files
# remove duplicates, keeping only the first
tx.drop_duplicates(subset='filename',keep='first',inplace=True)
rx.drop_duplicates(subset='filename',keep='first',inplace=True)
# drop tx Ex files
tx.drop(tx[tx.type=='RX'].index,inplace=True)
# drop tx 256 hz files
tx.drop(tx[tx.sampling_rate==256].index,inplace=True)
# drop invalid files
rx.drop(rx[~rx.valid].index,inplace=True)
tx.drop(tx[~tx.valid].index,inplace=True)
# drop aborted tx files (fewer than 30 seconds)
tx.drop(tx[tx.num_records<30].index,inplace=True)
# drop unassigned stations
rx.dropna(subset=['rx_station_qc'],inplace=True)
# TODO: drop bad tx files (user ID'd?)

# find TX-RX overlaps
print('Finding overlaps...')
overlaps = z3d_directory.find_overlaps(tx,rx,overlaps_csv=overlaps_file)
# trim one bad TX signal: 2019-07-24, 11:41:30, 0.5 Hz
# clip_time_209 = pd.Timestamp('2019-07-24 11:45:20').tz_localize('US/Mountain')
# overlaps.loc[overlaps.tx_ind==209,'end'] = clip_time_209

# Read calibration files
# cal_head = z3dio.read_syscal_header(cal_file)
print('Reading calibration files...')
cals = z3d_directory.read_zen_cals(cal_dir,ask_dir=False)
antcal = z3dio.read_antcal(antcal_file)

# store settings to be accessed in get_calibrated_fc
settings = {}
settings['notch_bandwidth'] = notch_bandwidth
settings['cycles_per_window'] = cycles_per_window
settings['strides_per_cycle'] = strides_per_cycle
settings['window_shape'] = window_shape
settings['cals'] = cals
settings['antcal'] = antcal

def get_calibrated_fc(overlap_data,transmitter_frequency,tx_mask,
                      sampling_rate,component,box_number,card_number,
                      antenna_number,settings):
    '''
    return calibrated Fourier coefficients, useful for computing transfer functions

    'settings' is a dictionary with these keywords:
        notch_bandwidth
        cycles_per_window
        strides_per_cycle
        window_shape
        cals
        antcal
    settings holds variables that don't change from overlap to overlap
    '''

    samples_per_cycle = int(round(sampling_rate/transmitter_frequency))

    # filter 60 Hz noise
    print('Filtering...')
    if transmitter_frequency > 16:
        # notch filter
        f0 = 60.0
        bw = settings['notch_bandwidth'] # -3 dB bandwidth
        Q = f0/bw
        w0 = f0*2./sampling_rate
        numerator,denominator = signal.iirnotch(w0,Q)
        filt_data = signal.filtfilt(numerator,denominator,overlap_data)
    else:
        # use mean running period
        filt_data = process.waveform_average_filter(
            overlap_data,sampling_rate=sampling_rate,mask=tx_mask)

    # compute robust running mean
    print('Applying drift correction...')
    taur = process.robust_running_mean(filt_data,samples_per_cycle)
    # apply drift correction
    drift = filt_data-taur

    # compute short time Fourier transform
    print('Computing Fourier coefficients...')
    samples_per_window = int(round(samples_per_cycle*settings['cycles_per_window']))
    stride=int(round(samples_per_cycle/settings['strides_per_cycle']))
    if settings['window_shape'] == 'tukey_ten_sample_taper':
        samples_per_window += 10
        ten_sample_taper = 20./samples_per_window
        window_shape = ('tukey',ten_sample_taper)
    else:
        window_shape = settings['window_shape']

    # stft_str = '{} sampling_rate, {} nperseg, {} noverlap'
    # print(stft_str.format(sampling_rate,samples_per_window,samples_per_window-stride))

    f,t,fc = signal.stft(drift,sampling_rate,window=window_shape,nperseg=samples_per_window,
                          noverlap=samples_per_window-stride,boundary=None,padded=False)
    num_windows = len(t)
    num_freq = len(f)

    # apply calibrations
    print('Applying calibrations...')
    try:
        zen_cal = settings['cals'].loc[box_number]
        fccb = process.apply_board_cal(fc,f,zen_cal,card_number,sampling_rate)
    except KeyError:
        print('Zen {} board calibration not found'.format(box_number))
        fccb = fc
    if component[0] == 'H':
        fcc = process.apply_antenna_cal(fccb,f,settings['antcal'],antenna_number)
    else:
        fcc=fccb
    return (f,t,fcc)


# initialize results arrays
n_pairs = overlaps.shape[0]
processed = [False]*n_pairs
tx_stations = np.empty(n_pairs)
rx_stations = np.empty(n_pairs)
rx_runs = ['']*n_pairs
rx_components = ['']*n_pairs
tx_frequencies = np.empty(n_pairs)
sampling_rates = np.empty(n_pairs)
num_sampless = np.empty(n_pairs)
signal_durations = np.empty(n_pairs)
samples_per_cycles = np.empty(n_pairs)
num_cycless = np.empty(n_pairs)
tx_filenames = ['']*n_pairs
rx_filenames = ['']*n_pairs
dx1s = np.empty(n_pairs)
dy1s = np.empty(n_pairs)
dz1s = np.empty(n_pairs)
dx2s = np.empty(n_pairs)
dy2s = np.empty(n_pairs)
dz2s = np.empty(n_pairs)
azimuths = np.empty(n_pairs)
inclinations = np.empty(n_pairs)
num_nanss = np.empty(n_pairs)
num_clipss = np.empty(n_pairs)
max_clipfree_tx_cycless = np.empty(n_pairs)
min_currents = np.empty(n_pairs)
max_currents = np.empty(n_pairs)
tf1_real_LS = np.empty(n_pairs)
tf1_imag_LS = np.empty(n_pairs)
tf3_real_LS = np.empty(n_pairs)
tf3_imag_LS = np.empty(n_pairs)
tf5_real_LS = np.empty(n_pairs)
tf5_imag_LS = np.empty(n_pairs)
tf7_real_LS = np.empty(n_pairs)
tf7_imag_LS = np.empty(n_pairs)
tf9_real_LS = np.empty(n_pairs)
tf9_imag_LS = np.empty(n_pairs)
tf1_var_LS = np.empty(n_pairs)
tf3_var_LS = np.empty(n_pairs)
tf5_var_LS = np.empty(n_pairs)
tf7_var_LS = np.empty(n_pairs)
tf9_var_LS = np.empty(n_pairs)
tf1_real_RM = np.empty(n_pairs)
tf1_imag_RM = np.empty(n_pairs)
tf3_real_RM = np.empty(n_pairs)
tf3_imag_RM = np.empty(n_pairs)
tf5_real_RM = np.empty(n_pairs)
tf5_imag_RM = np.empty(n_pairs)
tf7_real_RM = np.empty(n_pairs)
tf7_imag_RM = np.empty(n_pairs)
tf9_real_RM = np.empty(n_pairs)
tf9_imag_RM = np.empty(n_pairs)
tf1_var_RM = np.empty(n_pairs)
tf3_var_RM = np.empty(n_pairs)
tf5_var_RM = np.empty(n_pairs)
tf7_var_RM = np.empty(n_pairs)
tf9_var_RM = np.empty(n_pairs)
n_iters = np.empty(n_pairs)
n_goods = np.empty(n_pairs)

# initialize temporary data storage objects
processed_tx_data = pd.DataFrame(index=tx.index,
                                 columns=['processed','num_samples','calibrated_fc','tx_mask'])
processed_tx_data['processed'] = False
ovr = process.overlap_reader(tx,rx,overlaps,read_all_tx=True,
                             start_offset=start_offset,end_offset=end_offset)

# iterate over pairs, or pick a pair
# for i_pair in [16411]:
# for i_pair in np.arange(93,100):
for i_pair in np.arange(n_pairs):
    try:
        print('Processing pair {}...'.format(i_pair))

        # extract relevant rows and info from tx and rx
        i_tx = overlaps.tx_ind.iloc[i_pair]
        i_rx = overlaps.rx_ind.iloc[i_pair]
        tx_row = tx.loc[i_tx]
        rx_row = rx.loc[i_rx]
        # get filenames
        tx_filename = tx_row.fullpath
        tx_filenames[i_pair] = tx_filename
        rx_filename = rx_row.fullpath
        rx_filenames[i_pair] = rx_filename
        # get geometries
        dx1s[i_pair] = rx_row.dx1
        dy1s[i_pair] = rx_row.dy1
        dz1s[i_pair] = rx_row.dz1
        dx2s[i_pair] = rx_row.dx2
        dy2s[i_pair] = rx_row.dy2
        dz2s[i_pair] = rx_row.dz2
        azimuths[i_pair] = rx_row.azimuth
        inclinations[i_pair] = rx_row.inclination
        # get station numbers
        tx_station = tx_row.channel_station
        tx_stations[i_pair] = tx_station
        rx_station = rx_row.rx_station_qc
        rx_stations[i_pair] = rx_station
        rx_run = rx_row.run_qc
        rx_runs[i_pair] = rx_run
        # get component
        rx_component = rx_row.component
        rx_components[i_pair] = rx_component
        # get tx frequency
        transmitter_frequency = tx_row.txfreq
        tx_frequencies[i_pair] = transmitter_frequency
        # get sampling rates
        tx_sampling_rate = round(tx_row.sampling_rate)
        rx_sampling_rate = round(rx_row.sampling_rate)
        assert tx_sampling_rate == rx_sampling_rate, 'tx and rx have different sampling rates'
        sampling_rate = rx_sampling_rate
        sampling_rates[i_pair] = sampling_rate

        # get overlap data and info
        # print('Reading z3d data...')
        # (tx_overlap_data,rx_overlap_data) = process.overlap_data(tx,rx,overlaps,i_pair)
        (tx_overlap_data,rx_overlap_data) = ovr.next_overlap_data(i_pair)
        num_samples = len(tx_overlap_data)
        num_sampless[i_pair] = num_samples
        # signal_time = np.arange(num_samples)/sampling_rate
        signal_duration = num_samples/sampling_rate
        signal_durations[i_pair] = signal_duration
        # assume sampling_rate = n*transmitter_frequency
        samples_per_cycle = int(round(sampling_rate/transmitter_frequency))
        samples_per_cycles[i_pair] = samples_per_cycle
        num_cycles = num_samples//samples_per_cycle
        num_cycless[i_pair] = num_cycles

        # check that TX Sense was applied correctly
        min_currents[i_pair] = min(tx_overlap_data)
        max_currents[i_pair] = max(tx_overlap_data)
        print('Transmitter current max={}, min={}'.format(max(tx_overlap_data),
                                                          min(tx_overlap_data)))

        # check for nans
        tx_nans = np.isnan(tx_overlap_data)
        rx_nans = np.isnan(rx_overlap_data)
        nans = np.where(tx_nans | rx_nans)[0]
        num_nans = len(nans)
        num_nanss[i_pair] = num_nans
        frac_nans = num_nans / num_samples
        max_nanfree_sec = (max(np.diff(np.append(np.append([-1],nans),num_samples)))-1)/sampling_rate
        if len(nans)>0:
            nan_message = 'Total: {} nans, {} of signal is nan, {} seconds nan-free'
            print(nan_message.format(num_nans, frac_nans, max_nanfree_sec))
            print('Interpolating over nans...')
            # interpolate over nans
            ttt = lambda z: z.nonzero()[0]
            tx_overlap_data[tx_nans] = np.interp(ttt(tx_nans),ttt(~tx_nans),
                                                 tx_overlap_data[~tx_nans])
            rx_overlap_data[rx_nans] = np.interp(ttt(rx_nans),ttt(~rx_nans),
                                                 rx_overlap_data[~rx_nans])

        # check for clipping
        tx_clips = np.abs(tx_overlap_data)>=39.8
        rx_clips = np.abs(rx_overlap_data)>=1.99
        clips = np.where(tx_clips | rx_clips)[0]
        # num_tx_clips = len(tx_clips)
        # num_rx_clips = len(rx_clips)
        num_clips = len(clips)
        num_clipss[i_pair] = num_clips
        # frac_tx_clips = num_tx_clips / num_samples
        # frac_rx_clips = num_rx_clips / num_samples
        frac_clips = num_clips / num_samples
        clip_pad = np.append(np.append([-1],clips),num_samples)
        num_clipfree_samples = np.diff(clip_pad)-1
        max_clipfree_index = np.argmax(num_clipfree_samples)
        max_clipfree_samples = num_clipfree_samples[max_clipfree_index]
        max_clipfree_sec = max_clipfree_samples/sampling_rate
        max_clipfree_tx_cycles = max_clipfree_sec*transmitter_frequency
        max_clipfree_tx_cycless[i_pair] = max_clipfree_tx_cycles

        if len(clips)>0:
            clip_message = 'Total: {} clips, {} clipped, {} seconds clip-free'
            print(clip_message.format(num_clips, frac_clips, max_clipfree_sec))
            if max_clipfree_tx_cycles >= 3:
                # trim to what's clip-free
                print('Trimming!')
                rx_overlap_data = rx_overlap_data[clip_pad[max_clipfree_index]+1:clip_pad[max_clipfree_index+1]]
                tx_overlap_data = tx_overlap_data[clip_pad[max_clipfree_index]+1:clip_pad[max_clipfree_index+1]]
            else:
                raise ValueError('Insufficient clip-free data')

        # Is this tx series the one to reuse?
        # check to see if this overlap includes the full tx time series
        full_tx = (tx_row.start==overlaps.start.iloc[i_pair] and
                   tx_row.end==overlaps.end.iloc[i_pair])
        # if RX signal clips, TX needs to be reprocessed
        # rx_is_clipless = len(np.where(rx_clips)[0])==0
        rx_is_clipless = np.sum(rx_clips)==0
        # check if this tx series is the same length as the last time it was processed
        same_num_samples = (processed_tx_data.num_samples[i_tx] == len(tx_overlap_data))
        # has this tx series been processed before
        tx_processed = processed_tx_data.processed[i_tx] 
        if full_tx and rx_is_clipless and tx_processed and not same_num_samples:
            # I don't know how this would ever happen
            print('Error: Unexpectedly different number of tx samples; must reprocess')

        # reuse tx if appropriate
        if (tx_processed and rx_is_clipless and 
            full_tx and same_num_samples):
            # if tx has been processed and no rx clips and 
            # the rx series fully overlaps the tx series and the tx series is the same 
            # length as the last time it was processed
            ft,tt,fctc = processed_tx_data.calibrated_fc[i_tx]
            tx_mask = processed_tx_data.tx_mask[i_tx]
        else:
            # tx must be processed
            tx_mask = process.mask_near_on_time(
                tx_overlap_data,transmitter_frequency,sampling_rate,
                pad_before,pad_after)
            ft,tt,fctc = get_calibrated_fc(tx_overlap_data,transmitter_frequency,
                                           tx_mask,sampling_rate,tx_row.component,
                                           tx_row.box_number,tx_row.card_number,
                                           tx_row.antenna_number,settings)
            if rx_is_clipless and full_tx:
                # only store results if rx has no clips and therefore tx has not been trimmed
                if tx_processed:
                    print('Error: reprocessing previously processed tx')
                processed_tx_data.at[i_tx,'processed'] = True
                processed_tx_data.at[i_tx,'num_samples'] = len(tx_overlap_data)
                processed_tx_data.at[i_tx,'calibrated_fc'] = (ft,tt,fctc)
                processed_tx_data.at[i_tx,'tx_mask'] = tx_mask

        # process rx
        fr,tr,fcrc = get_calibrated_fc(rx_overlap_data,transmitter_frequency,
                                       tx_mask,sampling_rate,rx_component,
                                       rx_row.box_number,rx_row.card_number,
                                       rx_row.antenna_number,settings)

        # mean and covariance
        odd_harmonics = [1,3,5,7,9]
        i_fundamental = np.where(ft==transmitter_frequency)[0][0]
        odd_harmonics = [odd_harmonic*i_fundamental for odd_harmonic in odd_harmonics]
        # remove harmonics if there aren't enough processed frequencies
        odd_harmonics = odd_harmonics[:(len(ft)//2)]
        odd_harmonics = [odd_harmonic for odd_harmonic in odd_harmonics if odd_harmonic<len(ft)]
        num_odd_harmonics = len(odd_harmonics)
        # print(fctc.shape,fcrc.shape,odd_harmonics)
        tf_hat_LS_r, tf_var_LS_r = process.LS(fctc[odd_harmonics],fcrc[odd_harmonics])
        # save results, making sure that results are length 5
        tf_hat_LS = np.full(5,np.nan+1j*np.nan,dtype=complex)
        tf_var_LS = np.full(5,np.nan)
        tf_hat_LS[:num_odd_harmonics] = tf_hat_LS_r
        tf_var_LS[:num_odd_harmonics] = tf_var_LS_r
        tf1_real_LS[i_pair] = tf_hat_LS[0].real
        tf1_imag_LS[i_pair] = tf_hat_LS[0].imag
        tf3_real_LS[i_pair] = tf_hat_LS[1].real
        tf3_imag_LS[i_pair] = tf_hat_LS[1].imag
        tf5_real_LS[i_pair] = tf_hat_LS[2].real
        tf5_imag_LS[i_pair] = tf_hat_LS[2].imag
        tf7_real_LS[i_pair] = tf_hat_LS[3].real
        tf7_imag_LS[i_pair] = tf_hat_LS[3].imag
        tf9_real_LS[i_pair] = tf_hat_LS[4].real
        tf9_imag_LS[i_pair] = tf_hat_LS[4].imag
        tf1_var_LS[i_pair] = tf_var_LS[0]
        tf3_var_LS[i_pair] = tf_var_LS[1]
        tf5_var_LS[i_pair] = tf_var_LS[2]
        tf7_var_LS[i_pair] = tf_var_LS[3]
        tf9_var_LS[i_pair] = tf_var_LS[4]
        print('Least squares: done! Starting regression m-estimate...')

        # regression_m estimates
        tf_hat_r, tf_cov_r, n_iter, n_good = process.regression_m_slow(fctc[odd_harmonics],
                                                                       fcrc[odd_harmonics])
        print('Regression m-estimate done!')
        # make sure that results are length 5
        tf_hat = np.full((5,2),np.nan)
        tf_cov = np.full((5,2,2),np.nan)
        tf_hat[:num_odd_harmonics] = tf_hat_r
        tf_cov[:num_odd_harmonics] = tf_cov_r

        # save all results to data structure
        tf1_real_RM[i_pair] = tf_hat[0,0]
        tf1_imag_RM[i_pair] = tf_hat[0,1]
        tf3_real_RM[i_pair] = tf_hat[1,0]
        tf3_imag_RM[i_pair] = tf_hat[1,1]
        tf5_real_RM[i_pair] = tf_hat[2,0]
        tf5_imag_RM[i_pair] = tf_hat[2,1]
        tf7_real_RM[i_pair] = tf_hat[3,0]
        tf7_imag_RM[i_pair] = tf_hat[3,1]
        tf9_real_RM[i_pair] = tf_hat[4,0]
        tf9_imag_RM[i_pair] = tf_hat[4,1]
        tf1_var_RM[i_pair] = tf_cov[0,0,0]
        tf3_var_RM[i_pair] = tf_cov[1,0,0]
        tf5_var_RM[i_pair] = tf_cov[2,0,0]
        tf7_var_RM[i_pair] = tf_cov[3,0,0]
        tf9_var_RM[i_pair] = tf_cov[4,0,0]
        n_iters[i_pair] = n_iter[0]
        n_goods[i_pair] = n_good[0]
        processed[i_pair] = True

        if save_coefficients:
            # Save Fourier Coefficients
            coeff_filepath = os.path.join(coefficients_dir,str(i_pair))
            np.save(coeff_filepath,np.dstack((fctc[odd_harmonics],fcrc[odd_harmonics])))

        # QC plot of transfer function values: first 5 odd harmonics
        # compute transfer function
        # tf = fcrc/fctc
        # for n_freq in np.arange(1,10,2): # 1,3,5,7,9
        #     fign,axn = plot.plot_short_time_fourier(tf[n_freq*cycles_per_window,:],times=tt,line_width=0.1)
        #     axn.scatter(tf_hat_LS[n_freq].real,tf_hat_LS[n_freq].imag,color='orange')
        #     axn.scatter(tf_hat[n_freq,0],tf_hat[n_freq,1],color='red')
        #     fign.show()
        # TODO: plot covariance ellipses

        # TODO: figure out if this is the L1 norm minimum or not
        # TODO: compare to scipy.optimize
        # TODO: look at # 'good' points, residuals, etc.

    except Exception as e:
        print('Error processing pair {}...'.format(i_pair))
        print(tx_filename, rx_filename)
        print(e)
        traceback.print_exc()

# save results in dataframe
# tups = list(zip(tx_stations,rx_stations,num_clipss,tx_frequencies))
# column_names = ['tx_station','rx_station','num_clips','tx_frequency']
# df = pd.DataFrame(tups,column_names)
df = pd.DataFrame()
df['processed'] = processed
df['tx_station'] = tx_stations
df['rx_station'] = rx_stations
df['rx_run'] = rx_runs
df['rx_component'] = rx_components
df['num_nans'] = num_nanss
df['num_clips'] = num_clipss
df['max_clipfree_tx_cycles'] = max_clipfree_tx_cycless
df['sampling_rate'] = sampling_rates
df['num_samples'] = num_sampless
df['signal_duration'] = signal_durations
df['samples_per_cycle'] = samples_per_cycles
df['num_cycles'] = num_cycless
df['tx_filename'] = tx_filenames
df['rx_filename'] = rx_filenames
df['dx1'] = dx1s
df['dy1'] = dy1s
df['dz1'] = dz1s
df['dx2'] = dx2s
df['dy2'] = dy2s
df['dz2'] = dz2s
df['azimuth'] = azimuths
df['inclination'] = inclinations
df['transmitter_frequency'] = tx_frequencies
df['min_current'] = min_currents
df['max_current'] = max_currents
df['tf1_real_LS'] = tf1_real_LS
df['tf1_imag_LS'] = tf1_imag_LS
df['tf3_real_LS'] = tf3_real_LS
df['tf3_imag_LS'] = tf3_imag_LS
df['tf5_real_LS'] = tf5_real_LS
df['tf5_imag_LS'] = tf5_imag_LS
df['tf7_real_LS'] = tf7_real_LS
df['tf7_imag_LS'] = tf7_imag_LS
df['tf9_real_LS'] = tf9_real_LS
df['tf9_imag_LS'] = tf9_imag_LS
df['tf1_var_LS'] = tf1_var_LS
df['tf3_var_LS'] = tf3_var_LS
df['tf5_var_LS'] = tf5_var_LS
df['tf7_var_LS'] = tf7_var_LS
df['tf9_var_LS'] = tf9_var_LS
df['tf1_real_RM'] = tf1_real_RM
df['tf1_imag_RM'] = tf1_imag_RM
df['tf3_real_RM'] = tf3_real_RM
df['tf3_imag_RM'] = tf3_imag_RM
df['tf5_real_RM'] = tf5_real_RM
df['tf5_imag_RM'] = tf5_imag_RM
df['tf7_real_RM'] = tf7_real_RM
df['tf7_imag_RM'] = tf7_imag_RM
df['tf9_real_RM'] = tf9_real_RM
df['tf9_imag_RM'] = tf9_imag_RM
df['tf1_var_RM'] = tf1_var_RM
df['tf3_var_RM'] = tf3_var_RM
df['tf5_var_RM'] = tf5_var_RM
df['tf7_var_RM'] = tf7_var_RM
df['tf9_var_RM'] = tf9_var_RM
df['n_iter'] = n_iters
df['n_good'] = n_goods

# TODO: write out results to avg file
df.to_csv(results_file)

#TODO: rotate elec and mag data
# Here, we need x and y and z (can't work with just one i_pair)
# Plus corresponding tripod inclinations and mag declination (one per campaign?)
# rotate.tri2geograph(tfhx,tfhy,tfhz,ix,iy,iz,dec,units='degrees'):

print('Done!')
