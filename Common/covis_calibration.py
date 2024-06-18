import numpy as np
import scipy
import scipy.interpolate
from Common.franc_garr import franc_garr

def covis_calibration(data_in, bfm, ping, calibrate, T, S, pH, lat, depth):
    # The following inputs are part of the structure "ping":
    f = ping['hdr']['xmit_freq']         # xmit freq (hz)
    tau = ping['hdr']['pulse_width']     # Pulse length (sec)
    sls = ping['hdr']['power_sel']       # source level setting (dB)
    c = ping['hdr']['sound_speed']       # m/s
    fsamp = ping['hdr']['sample_rate']   # Complex sampling frequency (Hz)
    gain = ping['hdr']['gain_sel']       # gain (dB)

    # The beamformer angles and ranges are taken from "swp"
    theta_bfm = bfm['angle']   # radians
    range = bfm['range']

    cal_mode = calibrate['mode'] # 'VSS' or 'TS'

    # Attenuation (Francois-Garrison, with z = 2500, S = 35, T = 2, pH = 7)
    alpha, cc = franc_garr(f, T, S, pH, depth,lat)
    alpha = alpha/1000 # attenuation coefficient (dB/m)

    ##################################################################################

    # Fixed calibration parameters
    cali_year = 2018
    if cal_mode == 'VSS':
        if f != 396000:
            print('Error: VSS calibration only available at 396 kHz')
            return
        
        # Beam pattern for fan-beam transmitter, units are dB normalized
        # to 1 at boresight, and negative if response is less than at boresight,
        # which is the direction in which source level is measured.

        xmit = np.loadtxt('Inputs/calibration_file_%4d/TC2160_Horizontal_Beam_396kHz_%d.txt' % (cali_year,cali_year))

        # Change ram angles to angles appropriate to beam patterns
        ram_angles = xmit[:,0]
        na = np.where(ram_angles > 180)
        ram_angles[na] = ram_angles[na] - 360
        
        theta_deg = -ram_angles
        theta = np.pi/180*theta_deg

        # Interpolate beam pattern
        #xmit_patt = np.interp(theta_bfm,theta,xmit[:,1]) 不要使用np的插值，x递减时会出问题
        xmit_patt = scipy.interpolate.interp1d(theta,xmit[:,1])(theta_bfm)

        # Source level from 2018 ATF calibration
        sl_table = np.array([
            [200 , 199.8+0.26],
            [205 , 204.7+0.26],
            [210 , 209.5+0.26],
            [215 , 214.2+0.26],
            [220 , 217.0+0.26],])
        
        # interpolate source level setting to get source level
        sl = scipy.interpolate.interp1d(sl_table[:,0], sl_table[:,1])(sls)

        # Transmit 3-dB vertical full beamwidth (deg)
        theta_vr = 1.07
        # Convert to radians
        theta_vr = np.pi*theta_vr/180
        # Convert to integrated beamwidth (rad) using Gaussian approximation
        Delta_v = 0.5*np.sqrt(np.pi/np.log(2))*theta_vr
        # Coefficients of polynomial fit to receiver sensitivity
        P1 = np.loadtxt('Inputs/calibration_file_%4d/rec_sens_396_poly_coeffs_%d.txt' % (cali_year,cali_year))
        # this form of load assumes an ascii file as input and a
        #   double as output; changing to a mat-file will
        #   break the code by producing a double (kgb 9/27/10)
        # Receive sensitivity for each beam (machine units/ muPa, gain = 0 dB, no
        # TVG)
        beam_nos = np.arange(1,257)
        rs = np.polyval(P1,beam_nos)
        
        # Coefficients of polynomial fit to integrated horizontal beamwidth
        # (reverberation index = rindex)
        P2 = np.loadtxt('Inputs/calibration_file_%4d/rindex_396_poly_coeffs_%d.txt' % (cali_year,cali_year))
        # this form of load assumes an ascii file as input and a
        #   double as output; changing to a mat-file will
        #   break the code by producing a double (kgb 9/27/10)
        # Delta_h = rindex
        Delta_h = np.polyval(P2,beam_nos)

    ##################################################################################

        [nsamp,nbeams] = data_in.shape

        # Define
        rsinv = 1/rs
        Delta_hinv = 1/np.sqrt(Delta_h)
        xmitinv = np.power(10,(-xmit_patt/20))

        # The following diagonal matrix, postmultiplying the matrix data_in,
        # makes all corrections that depend on beam number
        Beam_corr = np.array([rsinv*Delta_hinv*xmitinv])

        # Transmission loss depends on range
        tl = 20*np.log10(range) + 2*alpha*range

        # The following diagonal matrix, premultiplying the matrix data_in,
        # makes the tl corrections that depends on sample number (range)
        tl_corr = np.array([np.power(10,(tl/20))])

        # The factors that don't depend upon beam
        # number are collected in
        const = 10**(-(sl+gain)/20)/np.sqrt(0.5*c*tau*Delta_v)
        data_out = const*(tl_corr.T@Beam_corr)*data_in

    ##################################################################################
        
    else:
        # TS-Fan 和 TS-Wide 暂时缺少
        print('calibration error')
        


    return data_out