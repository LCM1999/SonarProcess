import numpy as np
import scipy
import scipy.signal

def covis_filter(data_in, filt, png):
    # Set filter defaults
    if not 'bw' in filt:
        filt['bw'] = 2
    if not 'type' in filt:
        filt['type'] = 'butterworth'
    if not 'order' in filt:
        filt['order'] = 4
    if not 'decimation' in filt:
        filt['decimation'] = 1

    # Apply Filter
    if filt['status'] == 'on':
        bw = filt['bw']    # Bandwdith is filt_bw/tau
        tau = png['hdr']['pulse_width']   # Pulse length (sec)
        fsamp = png['hdr']['sample_rate']    # Complex sampling frequency (Hz)
        order = filt['order']
        
        if filt['type'].lower() == 'butterworth':
            B, A = scipy.signal.butter(order, bw/fsamp/tau)
        else:
            print('Unknown filter type')

        # Lowpass filter
        data_out = scipy.signal.filtfilt(B, A, data_in.T)
        data_out = data_out.T

        # decimate
        R = filt['decimation']
        data_out = data_out[0:data_out.shape[0]:R,:]
    
        # update new sampling freq
        png['hdr']['sample_rate'] = fsamp / R 

    return data_out, filt, png