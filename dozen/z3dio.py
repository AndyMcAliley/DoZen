'''
Functions to read data from Zonge z3d files and Zonge calibration files
The most used function to read z3d files is z3dio.read_z3d(filename)


'''

import os
import math
import numpy as np
import struct
from statistics import mode

import pandas as pd
# compute predicted start/stop times
# compute predicted file size
# compare to actual file size

# time series record start flag
TSR_FLAG=b'\xff\xff\xff\x7f\x00\x00\x00\x80'

# default channel factor
CHFACTOR=9.536743164062e-10

def read_in_chunks(binaryfile,chunk_size=32768):
    ''' Iterable to read binary file in, one chunk at a time'''
    while True:
        data = binaryfile.read(chunk_size)
        if not data:
            break
        yield data

def find_first_record(binaryfile):
    ''' Find the position of the flag for the first data record'''
    binaryfile.seek(0)
    pos = -1
    for chunk in read_in_chunks(binaryfile):
        pos = chunk.find(TSR_FLAG)
        if pos != -1:
            break
    return pos

def find_last_record(binaryfile):
    ''' Find the position of the flag for the last data record'''
    chunk_size = 32768
    binaryfile.seek(0,os.SEEK_END)
    file_size = binaryfile.tell()
    if file_size > chunk_size:
        # ceiling is the same as upside-down floor division:
        imax = -(-file_size//chunk_size)
        for i in range(1,imax):
            binaryfile.seek(-i*chunk_size, os.SEEK_END)
            chunk = binaryfile.read(chunk_size)
            pos = chunk.rfind(TSR_FLAG)
            if pos != -1:
                pos = pos + file_size - i*chunk_size
                break
    else:
        binaryfile.seek(0)
        data = binaryfile.read()
        pos = chunk.rfind(TSR_FLAG)
    return pos

def parse_tsrecord(bytes_data,factor=2.048/(2**32)):
    '''Parse a time series record'''
    if bytes_data[:8] != TSR_FLAG:
        raise ValueError('Record does not begin with flag')
    record = {}
    record['GPS_timestamp'] = int.from_bytes(bytes_data[8:12],byteorder='little')
    # record['GPS_timestamp'] = np.frombuffer(bytes_data[8:12],dtype=np.uint32)
    # record['GPS_lat'] = np.frombuffer(bytes_data[12:20],dtype=np.float64)
    record['GPS_lat'] = struct.unpack('d',bytes_data[12:20])[0]
    record['GPS_lon'] = struct.unpack('d',bytes_data[20:28])[0]
    record['status'] = bytes_data[28:32]
    record['GPS_accuracy'] = int.from_bytes(bytes_data[32:36],byteorder='little')
    record['temperature'] = struct.unpack('f',bytes_data[36:40])[0]
    record['mystery_bytes'] = bytes_data[40:64]
    # record['ts_data'] = single(record['ts_data'] / factor)
    record['ts_data'] = np.single(np.frombuffer(bytes_data,dtype=np.int32,offset=64)*factor)
    return record

def read_next_tsrecord(binaryfile,factor=2.048/(2**32)):
    '''
    Read and parse the next record
    Assume 8 byte flag does not fall across chunk boundary
    Record may span two chunks, though
    '''
    # size of 4096 Hz record is slightly more than 16384 bytes
    # choose a number of bytes that is comfortably larger than that
    num_bytes = 65536
    found_start = False
    # uberchunk concatenates chunks, in case record spans multiple chunks
    uberchunk = []
    for chunk in read_in_chunks(binaryfile, chunk_size=num_bytes):
        flag_position = chunk.find(TSR_FLAG)
        if flag_position != -1:
            if found_start:
                # start flag found on previous iteration, so flag_position is next record's flag
                uberchunk = uberchunk + chunk[:flag_position]
                next_record = parse_tsrecord(uberchunk,factor)
                break
            else:
                found_start = True
                end_position = chunk.find(TSR_FLAG,flag_position+1)
                if end_position != -1:
                    next_record = parse_tsrecord(chunk[flag_position:end_position],factor)
                    break
                else:
                    # flag for following record not found; need to read more
                    uberchunk = chunk[flag_position:]
    return next_record

def read_header(filename):
    """Read header information from Z3D file"""
    # check file length
    filesize = os.path.getsize(filename)
    if filesize < 1024:
        raise Exception('File size is less than minimum size of 1024 bytes')
    with open(filename, 'rb') as binaryfile:
        binaryfile.seek(0)
        header = [[],[]]
        header[0] = binaryfile.read(512)
        header[1] = binaryfile.read(512)
        # headervars = [{},{},{},{}]
        headervars = {}
        for iheader in range(0,2):
            headerstring = header[iheader].decode(encoding='ascii')
            headerlines = headerstring.splitlines()
            for line in headerlines:
                name,var_string = line.partition('=')[::2]
                key = name.strip()
                if key != '' and key.isprintable():
                    if key in headervars:
                        raise Exception('duplicate key '+key+' found in header of '+filename)
                    else:
                        # store variable as float if possible
                        try:
                            var = float(var_string.strip())
                        except ValueError:
                            var = var_string.strip()
                        headervars[key] = var
        # TODO: What if there's no = in the line?
        # Are spaces in names ok?
        
        # check to see if this is tx or rx
        is_TX = False
        if 'Tx.Freq' in headervars:
            is_TX = True
        
        # read in the rest of the metadata blocks
        metadata = ''
        end_of_metadata = False
        while not(end_of_metadata):
            block = binaryfile.read(512)
            ID = ''
            try:
                blockstring = block.decode(encoding='ascii')
                ID = blockstring[:37]
            except Exception as e:
                # print(repr(e))
                print('Reached end of header before expected. Attempting to continue...')
            if ID != '\n\n\nGPS Brd339/Brd357 Metadata Record\n':
                end_of_metadata = True
                binaryfile.seek(-512,1)
            else:
                # ID string is 37 characters, records end with \x00
                blocklines = blockstring[37:-1].splitlines()
                metadata = metadata+blocklines[0]
        metadata_pairs = metadata.split('|')
        for metadata_pair in metadata_pairs:
            name,var_string = metadata_pair.partition('=')[::2]
            key = name.strip()
            var_string = var_string.strip()
            # print('-'+var_string+'-')
            if key != '' and key.isprintable():
                if key in headervars:
                    # raise Exception('Duplicate key '+key+' found in header of '+filename)
                    print('Duplicate key '+key+' found in header of '+filename)
                # else:
                try:
                    var = float(var_string)
                except ValueError:
                    var = var_string
                headervars[key] = var
        
        # read system calibrate record
        # use \x00 as ending flag
        syscalrecord = binaryfile.read(32768)
        endbyte = syscalrecord.find(b'\x00')
        if endbyte == -1:
            raise ValueError('no \x00 ending flag found in system calibrate record')
        syscalrecord = syscalrecord[:endbyte]
        # binaryfile.seek(-32768+endbyte+1,1)
        # TODO: parse syscal
        if 'syscal' in headervars:
            raise Exception('Duplicate key syscal found in header of '+filename)
        headervars['syscal'] = syscalrecord
    return headervars

def get_file_info(filename,include_record_lengths=False):
    file_size = os.path.getsize(filename)
    file_info = {}
    file_info['metadata'] = read_header(filename)
    if include_record_lengths:
        file_info['record_lengths'] = get_record_lengths(filename)
    # channel ADC factor
    # gain = file_info['metadata']['A/D Gain']
    # ch_factor = 2.048/(2**gain * 2**32)
    try:
        ch_factor = file_info['metadata']['Ch.Factor']
    except:
        ch_factor = CHFACTOR
        file_info['metadata']['Ch.Factor'] = ch_factor
    # try:
    #     ch_station = file_info['metadata']['CH.STN']
    # except:
    #     file_info['metadata']['CH.STN'] = file_info['metadata']['RX.STN']
    # try:
    #     ch_offset1 = file_info['metadata']['CH.OFFSET.XYZ1']
    # except:
    #     file_info['metadata']['CH.OFFSET.XYZ1'] = file_info['metadata']['CH.XYZ1']
    # if file_info['metadata']['CH.CMP'][0] == 'E':
    #     try:
    #         ch_offset2 = file_info['metadata']['CH.OFFSET.XYZ2']
    #     except:
    #         file_info['metadata']['CH.OFFSET.XYZ2'] = file_info['metadata']['CH.XYZ2']

    scaling_factor = ch_factor
    # TODO: apply board cal
    # TODO: apply antenna cal if applicable
    if 'TX.SENSE' in file_info['metadata']:
        ohm_sense = file_info['metadata']['TX.SENSE']
        scaling_factor = ch_factor/ohm_sense
    # I don't know why I need to add the 24 extra bytes
    record_size = 8+4+8+8+4+4+4+file_info['metadata']['A/D Rate']*4+24
    with open(filename, 'rb') as binaryfile:
        # read first record
        binaryfile.seek(0)
        first_record_pos = find_first_record(binaryfile)
        file_info['first_record_pos'] = first_record_pos
        if first_record_pos == -1:
            # raise ValueError('No record flags found')
            file_info['num_records'] = 0
        else:
            binaryfile.seek(first_record_pos)
            file_info['first_record'] = read_next_tsrecord(binaryfile,scaling_factor)
            # read second record (often it is missing 4 bytes)
            binaryfile.seek(first_record_pos+8)
            file_info['second_record'] = read_next_tsrecord(binaryfile,scaling_factor)
            # read last record
            last_record_pos = find_last_record(binaryfile)
            binaryfile.seek(last_record_pos)
            file_info['last_record'] = parse_tsrecord(binaryfile.read(),scaling_factor)
            file_info['num_records'] = (file_size - first_record_pos)/record_size
            # check for GPS week rollover
            # check fails if file holds more than 1 week of data
            if file_info['num_records'] >= 604800:
                file_info['GPS_week_rollover'] = True
                raise ValueError('More than 1 week''s worth of data; z3dio does not handle that file size yet')
            elif (file_info['last_record']['GPS_timestamp'] 
                < file_info['first_record']['GPS_timestamp']):
                file_info['GPS_week_rollover'] = True
            else :
                file_info['GPS_week_rollover'] = False
    return file_info

def get_flag_positions(file_data,asbyte=False):
    '''Find byte positions of record flags'''
    arraytype = type(file_data)
    if arraytype == np.ndarray:
        datatype = file_data.dtype
        if datatype == np.int32:
            int_flag1 = np.isin(file_data,np.frombuffer(TSR_FLAG[:4],dtype=np.int32))
            int_flag2 = np.isin(file_data,np.frombuffer(TSR_FLAG[4:8],dtype=np.int32))
            int_pos = np.where(np.logical_and(int_flag1[:-1],int_flag2[1:]))[0]
            # byte_pos = int_pos*4
        else:
            raise ValueError('Expected numpy array of dtype np.int32')
    elif arraytype == bytes:
        # not optimal: duplicates data
        # could use for loop and file_data.find(TSR_FLAG): memory-friendly but slow
        int_array = np.frombuffer(file_data,dtype=np.int32)
        int_flag1 = np.isin(int_array,np.frombuffer(TSR_FLAG[:4],dtype=np.int32))
        int_flag2 = np.isin(int_array,np.frombuffer(TSR_FLAG[4:8],dtype=np.int32))
        int_pos = np.where(np.logical_and(int_flag1[:-1],int_flag2[1:]))[0]
        # byte_pos = int_pos*4
    else:
        raise ValueError('Expected numpy array or bytes')
    if asbyte:
        return int_pos*4
    else:
        return int_pos


def _get_record_lengths(flag_positions,file_size,asbyte=False):
    '''
    Count number of data in each record in file
    If asbyte=true, returns number of bytes per record
    Otherwise, returns number of 4-byte blocks per record
    Includes metadata
    '''
    record_lengths = np.diff(np.concatenate((flag_positions,[file_size//4])))
    if asbyte:
        return record_lengths*4
    else:
        return record_lengths


def get_record_lengths(filename,asbyte=False):
    '''
    Count number of data in each record in file
    If asbyte=true, returns number of bytes per record
    Otherwise, returns number of 4-byte blocks per record
    Includes metadata
    '''
    file_size = os.path.getsize(filename)
    int_data = np.fromfile(filename,dtype=np.int32)
    flag_positions = get_flag_positions(int_data,asbyte=False)
    return _get_record_lengths(flag_positions,file_size,asbyte)


def contiguous_gps(filename):
    '''
    UNFINISHED
    Return the start and end times of the longest contiguous run of gps timestamps
    '''
    # get gps timestamps, in seconds
    gps_times = read_gps_times(filename)
    num_records = len(gps_times)
    # compute difference in seconds between consecutive timestamps
    record_durations = np.diff(gps_times)
    # find time skips
    skips = np.where(record_durations != 1.)[0]
    skips_pad = np.append(np.append([-1],skips),num_records-1)
    num_skipfree_samples = np.diff(skips_pad)-1
    max_skipfree_index = np.argmax(num_skipfree_samples)
    # max_skipfree_samples = num_skipfree_samples[max_skipfree_index]
    # max_skipfree_sec = max_skipfree_samples/sampling_rate
    # max_skipfree_tx_cycles = max_skipfree_sec*transmitter_frequency
    # max_skipfree_tx_cycless[i_pair] = max_skipfree_tx_cycles
    max_skipfree_start_gps = gps_times[skips_pad[max_skipfree_index]+1]
    max_skipfree_end_gps = gps_times[skips_pad[max_skipfree_index+1]]
    pass


def check_file(filename):
    '''
    Make sure Z3D file is the correct length and is readable
    returns a tuple: (is_valid,status)
    is_valid = True if the file is a valid z3d file
    status = number of missing bytes in the final record if file is valid
    Otherwise:
        status=-1 if the most common record length is not the expected record length
        status=-2 if the number of records based on GPS timestamps, file size, and flags are unequal
        status=-3 if there are any records of anomalous length, beside the last one
    '''
    is_valid = False
    status = 0
    fi = get_file_info(filename)
    record_lengths = get_record_lengths(filename,asbyte=True)
    # find the most common record length
    (unique_record_lengths,num_record_lengths) = np.unique(record_lengths[:-1],
                                                           return_counts=True)
    record_length_mode = unique_record_lengths[num_record_lengths.argmax()]
    # expected record length
    record_size = 8+4+8+8+4+4+4+fi['metadata']['A/D Rate']*4+24
    # is the most common record length equal to the expected record length?
    if record_length_mode == record_size:
        # compute number of records in 3 different ways
        elapsed_seconds = (fi['last_record']['GPS_timestamp']
                           -fi['first_record']['GPS_timestamp'])/1024.
        num_records_gps = int(np.round(elapsed_seconds)) + 1
        num_records_file_size = int(np.ceil(fi['num_records']))
        num_records_flags = len(record_lengths)
        # are all 3 counts of records equal?
        if num_records_gps == num_records_file_size and num_records_gps == num_records_flags:
            # are all records equal length, besides the last one?
            if len(unique_record_lengths) == 1:
                is_valid = True
                # if file is valid, status is the number of missing bytes in the last record
                status = record_size - record_lengths[-1]
            else:
                status = -3
        else:
            status = -2
    else:
        status = -1
    return (is_valid,status)


def _read_data(int_data,flag_positions,file_size,sampling_rate=-1):
    '''
    Read data from z3d file from int data
    return a single array containing only data'''
    record_header_length = 16
    # int_data = np.fromfile(filename,dtype=np.int32)
    data_starts = flag_positions + record_header_length
    # file_size = os.path.getsize(filename)
    data_ends = np.concatenate((flag_positions[1:],[file_size//4]))
    data_lengths = data_ends - data_starts
    num_records = len(data_starts)
    # guess sampling rate, if it isn't provided
    if sampling_rate == -1:
        sampling_rate = mode(data_lengths)
    # throw error if any data records are longer than the sampling rate
    if max(data_lengths)>sampling_rate:
        print('Warning: Records found with more samples than the sampling rate! Truncating.')
        # raise ValueError('Record found with more samples than the sampling rate')
        overlong_records = np.where(data_lengths>sampling_rate)
        data_lengths[overlong_records] = sampling_rate
        data_ends[overlong_records] = data_starts[overlong_records] + sampling_rate
    # # pad int_data to avoid indexing out of range during overwrite
    # # warning: this converts int_data into a float array
    # # int_data = np.concatenate([int_data,[np.nan]*sampling_rate])
    # more robust yet memory intensive: new array
    float_data = np.empty(sampling_rate*num_records)
    float_data.fill(np.nan)
    for i_record in np.arange(num_records):
        record_data = np.empty(sampling_rate)
        record_data.fill(np.nan)
        record_data[:data_lengths[i_record]] = int_data[data_starts[i_record]:data_ends[i_record]]
        stored_range_start = i_record*sampling_rate
        stored_range_end = (i_record+1)*sampling_rate
        # int_data[stored_range_start:stored_range_end] = record_data
        float_data[stored_range_start:stored_range_end] = record_data
    #TODO: deal with last record being incomplete
    # return int_data[:num_records*sampling_rate]
    return float_data

def read_data(filename,sampling_rate=-1,use_header=True):
    '''
    Read data from z3d file, return a single array containing only data
    '''
    if use_header:
        z3d_header = read_header(filename)
        sampling_rate = int(z3d_header['A/D Rate'])
        factor = z3d_header['Ch.Factor']/z3d_header['ChannelGain']
        # convert from volts to amps, if applicable
        if 'TX.SENSE' in z3d_header:
            factor = factor/z3d_header['TX.SENSE']
    else:
        factor = 1.
    int_data = np.fromfile(filename,dtype=np.int32)
    flag_positions = get_flag_positions(int_data,asbyte=False)
    file_size = os.path.getsize(filename)
    return _read_data(int_data,flag_positions,file_size,sampling_rate)*factor


def _read_gps_times(int_data,flag_positions):
    '''Read only GPS TimeStamps for each record from z3d file'''
    gps_time_offset = 2
    gps_time_pos = flag_positions + gps_time_offset
    # extract gps times from int_data
    # divide by 1024 to convert timestamp int32 to seconds
    gps_times = int_data[gps_time_pos]/1024.
    return gps_times

def read_gps_times(filename):
    '''Read only GPS TimeStamps for each record from z3d file'''
    int_data = np.fromfile(filename,dtype=np.int32)
    flag_positions = get_flag_positions(int_data,asbyte=False)
    return _read_gps_times(int_data,flag_positions)


def read_block(binaryfile):
    headerbinary = binaryfile.read(512)
    headerstring = headerbinary.decode(encoding='ascii')
    headerlines = headerstring.splitlines()
    return headerlines

def parse_block(paramlines):
    blockdict = {}
    for line in paramlines:
        # ignore lines without no '=' symbol or with invalid names
        if line.find('=')>0:
            name,var_string = line.partition('=')[::2]
            if name.strip() != '' and name.isprintable():
                # store variable as float if possible
                try:
                    var = float(var_string.strip())
                except ValueError:
                    var = var_string.strip()
                blockdict[iheader][name.strip()] = var

def read_ts(binaryfile):
    first_record_pos = find_first_record(binaryfile)
    binaryfile.seek(first_record_pos)
    file_size = x
    record_size = y
    num_records = (file_size - first_record_pos)/record_size

def read_file_full(filename):
    '''Read file record by record, storing all record metadata'''

def read_file_data(filename):
    '''Read file data only, ignoring metadata when possible'''

def read_z3d(filename):
    '''
    Read metadata and time-series data from z3d file.
    This returns a dictionary.
    To access time series data in a z3d:
        z3d = z3dio.read_z3d(filename)
        ts_data = z3d['data']
    Dictionary keys also include 'metadata' which contains header info,
    'gps_times' which is a list of gps timestamps for each record,
    and 'nan_samples', which lists missing data samples.
    'data' itself has no nans, but when a record is too short, zeros 
    fill in the missing data.
    Use 
    '''
    # read metadata
    z3d = get_file_info(filename)
    sampling_rate = int(z3d['metadata']['A/D Rate'])
    # only read file via np.fromfile once
    # get file size
    file_size = os.path.getsize(filename)
    # read intdata
    int_data = np.fromfile(filename,dtype=np.int32)
    # get flag positions, which all functions that access int_data need
    flag_positions = get_flag_positions(int_data,asbyte=False)

    # do all the things that DON'T modify int_data first
    # read gps times
    z3d['gps_times'] = _read_gps_times(int_data,flag_positions)
    # get record lengths 
    record_header_length = 16
    z3d['record_lengths'] = _get_record_lengths(flag_positions,file_size)-record_header_length
    # insert NaN at ends of incomplete record lengths
    # needed because read_data returns int array - can't hold NaNs
    nan_records = np.where(z3d['record_lengths']<sampling_rate)[0]
    num_nan_samples = sampling_rate*len(nan_records)-sum(z3d['record_lengths'][nan_records])
    nan_samples = np.empty(num_nan_samples,dtype=np.int)
    i_nan_sample = 0
    for nan_record in nan_records:
        num_nans = max(sampling_rate - z3d['record_lengths'][nan_record],0)
        nan_index_start = sampling_rate*nan_record+z3d['record_lengths'][nan_record]
        nan_indeces = np.arange(nan_index_start,nan_index_start+num_nans)
        nan_samples[i_nan_sample:i_nan_sample+num_nans]=nan_indeces
        i_nan_sample = i_nan_sample+num_nans
    z3d['nan_samples'] = nan_samples

    # finally, read data, which does modify int_data
    z3d['data'] = _read_data(int_data,flag_positions,file_size,sampling_rate)
    # insert NaN for min?
    # z3d['data'][z3d['data']==-2147483648] = numpy.nan
    # convert data from A/D counts to volts, scale by gain/attenuator
    z3d['data'] = z3d['data']*z3d['metadata']['Ch.Factor']/z3d['metadata']['ChannelGain']
    # convert from volts to amps, if applicable
    if 'TX.SENSE' in z3d['metadata']:
        z3d['data'] = z3d['data']/z3d['metadata']['TX.SENSE']
    return z3d


def read_syscal_header(filename):
    ''' Read header of Zen cal file '''
    headervars = {}
    with open(filename,'r') as f:
        # lines = f.read().splitlines()
        # for line in lines:
        for line in (x.strip() for x in f):
            # parse every line beginning with $
            if line[0] == '$':
                # if line.count(',') > 1:
                    # headervars['column_names'] = line[1:].split(',')
                # else:
                name,var_string = line[1:].partition(',')[::2]
                key = name.strip()
                var_string = var_string.strip()
                if key != '' and key.isprintable():
                    if key in headervars:
                        raise Exception('duplicate key '+key+' found in header of '+filename)
                    else:
                        # store variable as float if possible
                        try:
                            var = float(var_string)
                        except ValueError:
                            if var_string.count(',')>1:
                                var = var_string.split(',')
                            else:
                                var = var_string
                        headervars[key] = var
    return headervars


def read_syscal_raw(filename):
    '''
    Read a Zonge calibration file
    Return a dictionary with two items:
        metadata is itself a dictionary of file header
        calibration is a dataframe with calibration data
        calibrations for each board at each base cal freq and a/d freq
    Assumes 5 harmonics
    '''
    metadata = read_syscal_header(filename)
    cal_version = metadata['CAL.VER']
    valid_cal_versions = [29,31]
    if not cal_version in valid_cal_versions:
        raise Exception('Error reading '+filename+': CAL.VER = '+cal_version+', but DoZen can only read versions 29 and 31. Time to extend DoZen!')

    # specify column names
    # don't use CAL.LABEL because sometimes it doesn't have all labels
    # column_names = metadata['CAL.LABEL']
    column_names = ['CARDSNX','CALFREQ','ADFREQ','G0','ATTEN','MAG1','PHASE1','MAG3','PHASE3','MAG5','PHASE5','MAG7','PHASE7','MAG9','PHASE9']

    cal=pd.read_csv(filename,header=None,names=column_names,index_col=False,comment='$',skipinitialspace=True)

    # replace pure whitespace with nans
    cal.replace(r'^\s*$', np.nan, regex=True, inplace=True)
    # drop columns with all nans
    cal.dropna(axis='columns',how='all',inplace=True)

    # strip whitespace (maybe unnecessary)
    for key in cal:
        if cal[key].dtype == 'O':
            # assume all objects are strings
            cal[key]=cal[key].str.strip()

    # parse Card number
    if cal_version == 29:
        # OSU Zen
        card_str = cal[column_names[0]].str.split(' ',n=1,expand=True)
        cal[column_names[0]] = card_str[0]
        cal['HDWKEY'] = card_str[1]
    # elif cal_version == 31:
        # USGS Zen
        # card_no = fields[0].strip()
    # parse the rest

    zen_cal = {}
    zen_cal['metadata']=metadata
    zen_cal['calibration']=cal
    return zen_cal


def read_syscal(filename):
    '''
    Read and parse a Zonge Zen calibration file
    Return a dataframe with calibration info:
        board serial number, cal freq, A/D freq, magnitude, and phase
    Assumes 5 harmonics
    Ignore G0, ATTEN (assume both are 0)
    '''
    zcd = read_syscal_raw(filename)
    zc_raw = zcd['calibration']
    zc = zc_raw.loc[:,['CARDSNX','CALFREQ','ADFREQ','MAG1','PHASE1']]
    zc.rename(columns={'MAG1':'MAG','PHASE1':'PHASE'},inplace=True)
    # read harmonics
    for harmonic in [3,5,7,9]:
        mag_col = 'MAG' + str(harmonic)
        phase_col = 'PHASE' + str(harmonic)
        zc_harmonic = zc_raw.loc[:,['CARDSNX','CALFREQ','ADFREQ',
                                    mag_col,phase_col]]
        zc_harmonic.rename(columns={mag_col:'MAG',phase_col:'PHASE'},inplace=True)
        zc_harmonic.loc[:,'CALFREQ'] = zc_harmonic.loc[:,'CALFREQ']*harmonic
        zc = zc.append(zc_harmonic,sort=False)
    zc.sort_values(['CARDSNX','ADFREQ','CALFREQ'],inplace=True)
    zc = zc.loc[:,['CARDSNX','ADFREQ','CALFREQ','MAG','PHASE']]
    return zc


def read_antcal(filename):
    ''' 
    Read a Zonge magnetic antenna calibration file 
    Return a dataframe with calibration info:
        antenna serial number, frequency, magnitude, and phase
    '''
    block_flag = 'AMT antenna'
    comment_flag = block_flag.split(' ')[0]
    with open(filename,'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]

    # remove empty lines
    lines = [l for l in lines if l != '']
    is_header_line = [l.startswith(block_flag) for l in lines]
    if sum(is_header_line) == 0:
        raise Exception('No lines starting with "'+block_flag+'" found in '+filename)
    header_line_numbers = np.where(is_header_line)[0]
    headers = [lines[i] for i in header_line_numbers]
    cal_frequencies_unique = [float(h.lstrip(block_flag)) for h in headers]
    # replace with nearest power of 2
    true_cal_frequencies = [2**(round(math.log(freq,2))) for freq in cal_frequencies_unique]

    # read file into dataframe
    column_names=['antenna_sn','mag','phz','mag 3/2','phz 3/2']
    df = pd.read_csv(filename,names=column_names,header=None,index_col=False,sep='\s+',engine='python',comment=comment_flag,skipinitialspace=True)

    # construct column of calibration frequencies to add to data frame
    cal_frequencies = np.empty(len(df))
    block_dividers = np.append(header_line_numbers,len(lines))
    num_block_lines = np.diff(block_dividers)-1
    i_start = 0
    for i_freq,hln in enumerate(header_line_numbers):
        i_end = i_start + num_block_lines[i_freq]
        cal_frequencies[i_start:i_end] = true_cal_frequencies[i_freq]
        i_start = i_end
    # add calibration frequencies to data frame
    df = df.assign(base_frequency = cal_frequencies)

    # construct new data structure, parsing 3/2 frequency columns
    df1 = df.loc[:,['antenna_sn','base_frequency','mag','phz']]
    df1p5 = df.loc[:,['antenna_sn','base_frequency','mag 3/2','phz 3/2']]
    df1p5.rename(columns={'mag 3/2':'mag','phz 3/2':'phz'},inplace=True)
    df1p5.loc[:,'base_frequency'] = df1p5.loc[:,'base_frequency']*1.5
    df_full = df1.append(df1p5,sort=False)
    df_full.sort_values(['antenna_sn','base_frequency'],inplace=True)
    return df_full


