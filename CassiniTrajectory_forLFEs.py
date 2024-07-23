"""
Created on Thu Jul 11 11:35:51 2024
@author: mrutala
"""

import spiceypy as spice
import numpy as np
import pandas as pd
import datetime
import tqdm


basefilepath = '/Users/mrutala/projects/others/CassiniLFEs/'

#   CASSINI end date: Sep. 15 2017
start_datetime = datetime.datetime(2004, 1, 1)
stop_datetime = datetime.datetime(2017, 9, 15, 11, 57) # datetime.datetime(2017, 9, 16)

#   Get an array of datetimes spanning start to stop
minutely_datetimes = np.arange(start_datetime, stop_datetime, 
                               datetime.timedelta(minutes=1),
                               dtype = datetime.datetime)

def get_CassiniEphemeris(datetimes):
    """
    

    Parameters
    ----------
    datetimes : list of datetimes
        A list of python datetimes for which the ephemeris will be returned

    Returns
    -------
    df : pandas DataFrame with columns [datetime, x_KSM, y_KSM, z_KSM, 
                                        subLat, subLon, subLST, R_KSM]
        The Cassini ephemeris for the input datetimes, including:
            - x, y, z, R in the KSM frame
            - the sub-spacecraft Latitude, Longitude, and Local Solar Time (LST)

    """
    #   Load a CASSINI SPICE metakernel 
    with spice.KernelPool('/Users/mrutala/SPICE/cassini/metakernel_cassini.txt'):
        
        #   Get SPICE code for Saturn and 3D radii
        target_id = str(spice.bodn2c('Saturn'))
        _, R_S_3 = spice.bodvrd(target_id, 'RADII', 3)
        saturn_flattening = (R_S_3[0] - R_S_3[2]) / R_S_3[0]
        
        #   Convert datetimes to ETs for SPICE
        ets = spice.datetime2et(datetimes)
        
        #   Query the spacecraft position with SPICE
        pos, lt = spice.spkpos('CASSINI', 
                               ets, 
                               'CASSINI_KSM', 
                               'None', 
                               target_id)
        
        sub_lats = []
        sub_lons = []
        sub_lsts = []
        for et in tqdm.tqdm(ets):
            
            #   Get the sub-observer (i.e. spacecraft) point on the target (planet)
            #   This does not handle array inputs, hence the loop
            subpoint, epoch, vec = spice.subpnt('INTERCEPT/ELLIPSOID',  # Method
                                                target_id,              # Saturn
                                                et,                     # When
                                                'IAU_SATURN',           # Frame
                                                'None',                 # Light travel time correction
                                                'CASSINI')              # Observer
            
            #   Convert from rectangular IAU_SATURN coords to planetographic
            lon, lat, alt = spice.recpgr(target_id, 
                                         subpoint, 
                                         R_S_3[0], 
                                         saturn_flattening)
            
            sub_lats.append(lat)
            sub_lons.append(lon)
            
            #   Convert from time to local solar time (LST), given a longitude
            lst_tuple = spice.et2lst(et, 
                                     int(target_id), 
                                     lon, 
                                     'PLANETOGRAPHIC')
            #   Convert local solar time to decimal hours
            decimal_hour_lst = lst_tuple[0] + lst_tuple[1]/60 + lst_tuple[2]/3600
            sub_lsts.append(decimal_hour_lst)
    
    #   Put everything in a dataframe
    df = pd.DataFrame({'x_KSM':pos.T[0] / R_S_3[0], 
                       'y_KSM':pos.T[1] / R_S_3[0],
                       'z_KSM':pos.T[2] / R_S_3[0], 
                       'subLat':np.array(sub_lats) * 180/np.pi,
                       'subLon':np.array(sub_lons) * 180/np.pi,
                       'subLST':np.array(sub_lsts)},
                      index = datetimes)
    
    #   Calculate the radial distances (this could be done with SPICE, too)
    df['R_KSM'] = np.sqrt(np.sum(df[['x_KSM', 'y_KSM', 'z_KSM']]**2, axis=1))
    
    return df

#   Print this to a csv (this will be pretty hefty)
eph_df = get_CassiniEphemeris(minutely_datetimes)
eph_df.to_csv(basefilepath + '/20040101000000_20170915115700_ephemeris.csv', index_label='datetime')

#   Now query times based on the original LFE list and write
LFEs_rawfilepath = basefilepath + '/2004001_2017258_start_stop_times.csv'
LFEs_raw = pd.read_csv(LFEs_rawfilepath)
LFEs_raw_datetimes = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S.%f') for t in LFEs_raw['start']]

LFEs_raw_eph_df = get_CassiniEphemeris(LFEs_raw_datetimes)
LFEs_raw_eph_df.to_csv(basefilepath + '/2004001_2017258_ephemeris.csv',
                       index_label = 'datetime')

#Subsample for the joined LFE list and write
LFEs_joinedfilepath = basefilepath + '/LFEs_joined.csv'
LFEs_joined = pd.read_csv(LFEs_joinedfilepath)
LFEs_joined_eph_datetimes = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S.%f') for t in LFEs_joined['start']]

LFEs_joined_eph_df = get_CassiniEphemeris(LFEs_joined_eph_datetimes)
LFEs_joined_eph_df.to_csv(basefilepath + '/LFEs_joined_ephemeris.csv',
                          index_label = 'datetime')



