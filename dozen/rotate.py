import sys
import os
from statistics import mode

import numpy as np
import pandas as pd

from . import z3dio
from . import util

# global script variables
script_instructions = """Usage: python rotate.py [-d] input.inp

To rotate z3ds:
python rotate.py input_file.inp

To generate a dummy input file dummy.inp:
python rotate.py -d dummy.inp
"""


def correct_inclinations(incl_x,incl_y,incl_z,tol=0.1,units='rad'):
    # adjust inclinations so that sin^2(x)+sin^2(y)+sin^2(z)=1
    scale = util.radian_conversion_factor(units)
    incl_x_rad = incl_x*scale
    incl_y_rad = incl_y*scale
    incl_z_rad = incl_z*scale
    six = np.sin(incl_x_rad)
    siy = np.sin(incl_y_rad)
    siz = np.sin(incl_z_rad)
    i_norm = np.sqrt(six**2+siy**2+siz**2)
    if abs(i_norm-1) > tol:
        print('Warning: sine inclination angles have error above {}%'.format(tol*100))
    ax = np.arcsin(six/i_norm)/scale
    ay = np.arcsin(siy/i_norm)/scale
    az = np.arcsin(siz/i_norm)/scale
    return (ax,ay,az,i_norm)

def sine_squared_corrected_inclinations(incl_x,incl_y,incl_z,tol=0.2):
    # adjust inclinations so that sin^2(x)+sin^2(y)+sin^2(z)=1
    # return as sine squared (computationally efficient)
    ssix = np.sin(incl_x)**2
    ssiy = np.sin(incl_y)**2
    ssiz = np.sin(incl_z)**2
    i_norm_sq = ssix+ssiy+ssiz
    if abs(i_norm_sq-1) > tol:
        print('Warning: sine squared inclination angles have error above {}%'.format(tol*100))
    ssax = ssix/i_norm_sq
    ssay = ssiy/i_norm_sq
    ssaz = ssiz/i_norm_sq
    return (ssax,ssay,ssaz,i_norm_sq)

def geomag2tri_matrix(incl_x,incl_y,incl_z,units='rad'):
    # rotation matrix to transform from geomagnetic coords to tripod coords
    # convert to radians
    scale = util.radian_conversion_factor(units)
    incl_x_rad = incl_x*scale
    incl_y_rad = incl_y*scale
    incl_z_rad = incl_z*scale
    # get sine squared of each corrected angle
    ssax,ssay,ssaz,i_norm_sq = sine_squared_corrected_inclinations(incl_x_rad,incl_y_rad,incl_z_rad)
    # form square of rotation matrix
    rotsq = np.array([[ssay+ssaz,0,ssax],
                      [ssay*ssax/(1-ssax),ssaz/(1-ssax),ssay],
                      [ssaz*ssax/(1-ssax),ssay/(1-ssax),ssaz]])
    # square root, assign correct signs
    rot = np.sqrt(rotsq)
    # OSU convention
    # rot = np.multiply(rot, np.array([[ 1, 1,-1],
    #                                  [-1, 1,-1],
    #                                  [ 1, 1, 1]]))
    # USGS convention (flipped Hz)
    # rot = np.multiply(rot, np.array([[-1, 1, 1],
    #                                  [ 1,-1, 1],
    #                                  [-1,-1,-1]]))
    # USGS convention (no flipped Hz)
    rot = np.multiply(rot, np.array([[-1, 1, 1],
                                     [ 1,-1, 1],
                                     [ 1, 1, 1]]))
    return rot

def tri2geomag_matrix(incl_x,incl_y,incl_z,units='rad'):
    # rotation matrix to transform from tripod coords to geomagnetic coords
    # convert to radians
    scale = util.radian_conversion_factor(units)
    incl_x_rad = incl_x*scale
    incl_y_rad = incl_y*scale
    incl_z_rad = incl_z*scale
    # get sine squared of each corrected angle
    ssax,ssay,ssaz,i_norm_sq = sine_squared_corrected_inclinations(incl_x_rad,incl_y_rad,incl_z_rad)
    # form square of rotation matrix
    rotsq = np.array([[ssay+ssaz,ssax*ssay/(1-ssax),ssax*ssaz/(1-ssax)],
                      [0,ssaz/(1-ssax),ssay/(1-ssax)],
                      [ssax,ssay,ssaz]])
    # square root, assign correct signs
    rot = np.sqrt(rotsq)
    # OSU convention
    # rot = np.multiply(rot, np.array([[ 1, 1,-1],
    #                                  [-1, 1,-1],
    #                                  [ 1, 1, 1]]))
    # USGS convention (flipped Hz)
    # rot = np.multiply(rot, np.array([[-1, 1,-1],
    #                                  [ 1,-1,-1],
    #                                  [ 1, 1,-1]]))
    # USGS convention (no flipped Hz)
    rot = np.multiply(rot, np.array([[-1, 1, 1],
                                     [ 1,-1, 1],
                                     [ 1, 1, 1]]))
    return rot

def geomag2geograph_matrix(dec,units='rad'):
    '''
    Rotate from geomagnetic coords to geographic coords
    based on magnetic declination dec (east of north)
    '''
    # convert to radians
    scale = util.radian_conversion_factor(units)
    dec_rad = dec*scale
    cosdec = np.cos(dec_rad)
    sindec = np.sin(dec_rad)
    return np.array([[cosdec,-sindec,0.],
                    [sindec,cosdec,0.],
                    [0.,0.,1.]])

def tri2geograph_matrix(ix,iy,iz,dec,units='rad'):
    '''
    Rotate from tripod coords to geographic coords,
    based on magnetic declination dec
    NOTE: this assumes that ix,iy,iz are roughly correct, not flipped by 90-angle
    '''
    # convert to radians, depending on units
    scale = util.radian_conversion_factor(units)
    ix = scale*ix
    iy = scale*iy
    iz = scale*iz
    dec = scale*dec
    rt2g = tri2geomag_matrix(ix,iy,iz)
    rg2g = geomag2geograph_matrix(dec)
    R = np.dot(rg2g,rt2g)
    # form 3 by n matrix
    # V = np.array([vx[:],vy[:],vz[:]])
    # apply rotation
    # return np.dot(R,V)
    return R

def tri2geomag(vx,vy,vz,ix,iy,iz,units='rad'):
    '''
    Rotate vectors vx, vy, vz to geomagnetic coords
    NOTE: this assumes that ix,iy,iz are roughly correct, not flipped by 90-angle
    '''
    # convert to radians, depending on units
    scale = util.radian_conversion_factor(units)
    ix = scale*ix
    iy = scale*iy
    iz = scale*iz
    rt2g = tri2geomag_matrix(ix,iy,iz)
    # form 3 by n matrix
    V = np.array([vx[:],vy[:],vz[:]])
    # apply rotation
    return np.dot(rt2g,V)

def tri2geograph(vx,vy,vz,ix,iy,iz,dec,units='rad'):
    '''
    Rotate vectors vx, vy, vz to geographic coords,
    based on magnetic declination dec
    NOTE: this assumes that ix,iy,iz are roughly correct, not flipped by 90-angle
    '''
    R = tri2geograph_matrix(ix,iy,iz,dec,units=units)
    # form 3 by n matrix
    V = np.array([vx[:],vy[:],vz[:]])
    # apply rotation
    return np.dot(R,V)

def inc2azimuth(ix,iy,iz,dec,units='rad'):
    '''
    get azimuths of tripod coords from their inclinations and declination of x
    '''
    t2gg = tri2geograph_matrix(ix,iy,iz,dec,units=units)
    azimuths = np.arctan2(t2gg[1,:],t2gg[0,:])
    # convert back to units
    scale = util.radian_conversion_factor(units)
    azimuths /= scale
    return azimuths


class RecordReader:
    '''
    Read z3d files record by record
    '''
    def __init__(self,filename):
        self.z3d=filename
        # get record flag byte positions 
        self.flag_positions = z3dio.get_flag_positions(np.fromfile(filename,dtype=np.int32),
                                                       asbyte=True)
        self.flag_positions = np.array(self.flag_positions)
        # get record lengths
        file_size = os.path.getsize(filename)
        record_bounds = np.concatenate((self.flag_positions,[file_size]))
        self.record_lengths = np.diff(record_bounds)
        # get gps times
        self.gps_times = z3dio.read_gps_times(filename)
        self.valid_records = [True]*len(self.flag_positions)

    def iterrecord(self):
        '''
        Return a tuple at each iteration: header and data.
        For the first record, header includes all z3d metadata,
        and the first record's header.
        '''
        # get header lengths and record data lengths as byte counts from flag positions
        num_records = len(self.flag_positions)
        assert len(self.valid_records) == num_records, \
                'valid_records must be of same length as flag_positions'
        header_lengths = [self.flag_positions[0]+64]+[64]*(num_records-1)
        data_lengths = self.record_lengths-64
        # open file
        with open(self.z3d,mode='rb') as bf:
            for i_record, valid_record in enumerate(self.valid_records):
                header = bf.read(header_lengths[i_record])
                record_data_bytes = bf.read(data_lengths[i_record])
                record_data = np.single(np.frombuffer(record_data_bytes,dtype=np.int32))
                if valid_record:
                    yield (header,record_data)


def set_overlap(z1, z2, z3):
    '''
    Ensure that only overlapping records are returned by iterrecord
    z1, z2, and z3 are RecordReader objects
    This sets the valid_records property of each to only be True when 
    there is overlap in gps_times between z1, z2, and z3.
    '''

    overlaps_1 = [True]*len(z1.flag_positions)
    overlaps_2 = [True]*len(z2.flag_positions)
    overlaps_3 = [True]*len(z3.flag_positions)

    # use z1.gps_times, z2.gps_times, z3.gps_times to assign False values 
    # in these 3 arrays whenever all three gps_times do not overlap

    # assign those arrays
    z1.valid_records = overlaps_1
    z2.valid_records = overlaps_2
    z3.valid_records = overlaps_3


def rotate_df(df):
    '''
    Rotate z3ds, given inputs in dataframe df
    df has the following columns:
        file1
        file2
        file3
        inclination1
        inclination2
        inclination3
        declination1
    '''
    for row in df.itertuples(index=True,name='Pandas'):
        # count how many z3ds there are
        num_z3ds = 0
        # z3ds = []
        out_z3ds = []
        if row.file1 != '':
            num_z3ds += 1
            # z3ds.append(row.file1)
            name,ext = os.path.splitext(row.file1)
            out_z3ds.append(name+'_rot'+ext)
        if row.file2 != '':
            num_z3ds += 1
            # z3ds.append(row.file2)
            name,ext = os.path.splitext(row.file2)
            out_z3ds.append(name+'_rot'+ext)
        if row.file3 != '':
            num_z3ds += 1
            # z3ds.append(row.file3)
            name,ext = os.path.splitext(row.file3)
            out_z3ds.append(name+'_rot'+ext)
        if (num_z3ds != 3):
            print('{} z3d files listed; too few to rotate'.format(num_z3ds))
        else:

            z1 = RecordReader(row.file1)
            z2 = RecordReader(row.file2)
            z3 = RecordReader(row.file3)

            # check that gps times match among all z3ds

            # check that all record lengths are the same

            # assert all([(z1.gps_times.astype(np.int)==z2.gps_times.astype(np.int)).all(),
            #             (z1.gps_times.astype(np.int)==z3.gps_times.astype(np.int)).all()]), 'gps times mismatch'

            # synchronize data
            set_overlap(z1, z2, z3)

            # form rotation matrix
            # scale = util.radian_conversion_factor('degrees')
            # ix = scale*row.inclination1
            # iy = scale*row.inclination2
            # iz = scale*row.inclination3
            # dec = scale*row.declination1
            # rt2g = tri2geomag_matrix(ix,iy,iz)
            # rg2g = geomag2geograph_matrix(dec)
            # Rmat = np.dot(rg2g,rt2g)
            Rmat = tri2geograph_matrix(row.inclination1,
                                       row.inclination2,
                                       row.inclination3,
                                       row.declination1,
                                       units='degrees')
            # read, write out everything besides data verbatim
            with open(out_z3ds[0],'wb') as rw1, open(out_z3ds[1],'wb') as rw2, open(out_z3ds[2],'wb') as rw3:
                for record1,record2,record3 in zip(z1.iterrecord(),z2.iterrecord(),z3.iterrecord()):
                    # write headers
                    rw1.write(record1[0])
                    rw2.write(record2[0])
                    rw3.write(record3[0])

                    # apply rotation
                    nd = min([len(record1[1]),len(record2[1]),len(record3[1])])
                    V = np.array([record1[1][:nd],record2[1][:nd],record3[1][:nd]])
                    # is z sign change really hardwired in USGS equipment?
                    # V = np.array([record1[1][:nd],record2[1][:nd],-1*record3[1][:nd]])
                    rotated_data = np.dot(Rmat,V)
                    rotated_data1 = rotated_data[0].astype(np.int32)
                    rotated_data2 = rotated_data[1].astype(np.int32)
                    rotated_data3 = rotated_data[2].astype(np.int32)
                    # write rotated data
                    rotated_data1.tofile(rw1)
                    rotated_data2.tofile(rw2)
                    rotated_data3.tofile(rw3)



def main(argv):
    '''
    Entry point for script
    Parse command line arguments
    Print usage instructions
    Print empty input file
    '''

    if len(argv) == 2:
        try:
            df_input = pd.read_csv(argv[1])
        except:
            print('Error reading input file')
            exit()
    elif len(argv) == 3 and argv[1]=='-d':
        print('Writing dummy input file to '+argv[2])
        cwd = os.getcwd()
        dict_dummy={
            'file1':[cwd+'/Hx.z3d',cwd+'/Ex.z3d'],
            'file2':[cwd+'/Hy.z3d',cwd+'/Ey.z3d'],
            'file3':[cwd+'/Hz.z3d',''],
            'inclination1':[41,0],
            'inclination2':[34,0],
            'inclination3':[32,0],
            'declination1':[8.75,8.75]
        }
        df_dummy = pd.DataFrame(data=dict_dummy)
        df_dummy.to_csv(argv[2],index=False)
        exit()
    else:
        print(script_instructions)
        exit()
    rotate_df(df_input)


if __name__ == "__main__":
    main(sys.argv)

