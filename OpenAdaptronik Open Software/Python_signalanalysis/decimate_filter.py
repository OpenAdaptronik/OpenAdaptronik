import numpy as np
import scipy as sp

from scipy.signal import freqz

def decimate_filter(data,decimation_factor):
    """ decimate a time sequence by a factor of decimation_factor """
    """ apply lowpass of order filter_order to avoid aliasing """
    
    """ design filter"""
    Wpp = 0.6/decimation_factor

    [lowpass_num,lowpass_den] = sp.signal.butter(8,Wpp,'Lowpass',0,'ba')
  
    """ calculate filter FRF"""
    [w,h] = freqz(lowpass_num, lowpass_den)

    """ filter data"""
    lowpass_data_inter   = sp.signal.lfilter(lowpass_num,lowpass_den,data)
    flipped_samples      = np.arange(len(data)-1,0,-1)
    lowpass_data_reverse = sp.signal.lfilter(lowpass_num,lowpass_den,lowpass_data_inter[flipped_samples])
    flipped_samples      = np.arange(len(lowpass_data_reverse)-1,0,-1)
    lowpass_data         = lowpass_data_reverse[flipped_samples]
    """ decimate by downsampling"""
    decimation_samples   = np.arange(0,len(lowpass_data),decimation_factor)
    decimated_data       = lowpass_data[decimation_samples]
   
    return decimated_data
    