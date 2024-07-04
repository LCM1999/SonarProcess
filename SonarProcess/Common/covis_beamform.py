'''
注:只有fast, ftt和near暂时用不到

'''

import numpy as np
from Common import config

def covis_beamform(bfm, raw_sig):
    if not 'fc' in bfm:
        print('covis_beamform: error, fc must be defined')
        bf_sig = []
        return bfm, bf_sig
    
    # Set default type
    if not 'type' in bfm:
        bfm['type'] = 'fast'

    # Set default array length
    if not 'array_length' in bfm:
        bfm['array_length'] = 0.408    # array length (m)

    # Set default beams angles
    if not 'start_angle' in bfm:
        bfm['start_angle'] = -64
    if not 'end_angle' in bfm:
        bfm['end_angle'] = 64
    if not 'num_beams' in bfm:
        if bfm['fc'] == 200000:
            bfm['num_beams'] = 128
        if bfm['fc'] == 396000:
            bfm['num_beams'] = 256

    # Set default sound speed
    if not 'c' in bfm:
        bfm['c'] = 1500    # (m/s)

    nsamps = raw_sig.shape[0]
    nchans = raw_sig.shape[1]
    
    # set default beam angles
    # row vector containing nbeams angles at which beams are to be formed (deg)
    nbeams = bfm['num_beams']   # number of beams
    start_angle = bfm['start_angle']
    end_angle = bfm['end_angle']
    bfm['angle'] = (np.pi / 180) * np.linspace(start_angle, end_angle, nbeams)    # beam angles (rad)

    # define beamformer params
    type = bfm['type']
    angle = bfm['angle']
    L = bfm['array_length']    # array length (m)
    c = bfm['c']   # sound speed (m/s)
    fsamp = bfm['fs']    # sampling freq (hz)
    f = bfm['fc']   # center frequency (hz)
    lamb = c / f    # wave length
    k = 2 * np.pi / lamb    # wave number
    first_samp = 1

    # slant range of each sample - round trip
    dr = c / fsamp / 2
    start_range = dr * first_samp
    end_range = dr * nsamps
    range = np.arange(start_range, end_range + 0.1 * dr, dr)
    
    # set bfm structure range
    bfm['range'] = range

    # w1 = column vector containing nch shading coefficients
    w1 = np.hamming(nchans)
    
    lamb = c / f    # wave length
    k = 2 * np.pi / lamb    # wave number
    omegac = 2 * np.pi * f
    t_delay_max = L / 2 / c * np.sin(angle)
    delf = fsamp / nsamps
    freqs = delf * np.arange(0, nsamps)
    omegas = 2 * np.pi * freqs
    bf_sig = np.zeros((nsamps, nbeams))

    if config.Verbose > 2:
        print('covis_beamform: type %s, f=%.2f, fs=%.2f, c=%.2f' % type, f, fsamp,c )
        print('covis_beamform: nsamps=%d, nbeams=%d' % nsamps, nbeams)
        print('covis_beamform: start_range=%.2f, end_range=%.2f' % start_range, end_range)
    
    # phase shift beamform
    if type == 'fast':
        # phase shift to apply to each channel
        # phi is a (nchans x nbeams) matrix
        d = np.array([(L/2)*np.linspace(-1,1,nchans)]).T    # element spacing
        phi = k*d@np.array([np.sin(angle)])    # phase shift
        # phase shift and sum to form beams by matrix multiplication
        bf_sig = raw_sig@((np.array([w1]).T@np.ones((1,nbeams)))*np.exp(1j*phi))

    else:
        print('covis_beamform: UNKNOWN BEAMFORMER TYPE')

    return bfm, bf_sig

