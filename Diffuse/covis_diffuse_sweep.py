import numpy as np
import os
import gsw
from datetime import datetime
import scipy
import pandas as pd
import json
import scipy.interpolate
from Common import config
from Common.covis_rectgrid import covis_rectgrid
from Common.covis_read import covis_read
from Common.covis_phase_correct import covis_phase_correct
from Common.covis_filter import covis_filter
from Common.covis_beamform import covis_beamform
from Common.covis_calibration import covis_calibration
from Diffuse.covis_covar_hamming import covis_covar_hamming
from Common.l2grid import l2grid
from Common.covis_version import covis_version

def covis_diffuse_sweep(filedir, json_file = 'covis_diffuse.json'):
    head_tail = os.path.split(filedir)
    swp_path = head_tail[0]
    swp_name = head_tail[1]


    # -------------------------
    # Initialization
    # -------------------------

    # sonar's central yaw and heading
    year = swp_name[6:10]
    central_yaw = 135
    if year == '2018':
        central_head = 289
    else:
        central_head = 279

    # calibration parameters
    T = 2.41    # temperataure (degree C)
    S = 34.53    # salinity (PSU)
    pH = 7    # pH
    lat = 46    # latitude ( degree N )
    depth =1544    # depth ( m )
    p = gsw.p_from_z(-depth,lat)
    c = gsw.sound_speed_t_exact(S,T,p)    # sound speed ( m/s )

    # mask parameters
    noise_floor = 0.64    # rms noise floor (uncalibrated in machine units)
    snr_thresh = 45    # snr threshold ( dB )

    # correlation time lag (sec)
    tlag = 2

    # averaging window length (sec)
    avg_win = 4

    # bathymetry data directory
    swp_date = datetime.strptime(swp_name[6:21],'%Y%m%dT%H%M%S')

    if swp_date <= datetime(2019,7,6):
        bathy_file = 'covis_bathy_2018.mat'
    elif swp_date <= datetime(2019,11,23):
        bathy_file = 'covis_bathy_2019a.mat'
    else:
        bathy_file = 'covis_bathy_2019b.mat'

    bathy = scipy.io.loadmat('Inputs/' + bathy_file)
    

    # parse sweep.json file in data archive
    swp_file = 'sweep.json'
    if not os.path.exists(filedir + '/' + swp_file):
        print('sweep.json file does not exist in %s\n' % filedir)
        return

    with open(filedir + '/' + swp_file, 'r', encoding='utf-8') as f:  
        swp = json.load(f) 

    # save sweep path and name in swp structure
    swp['path'] = swp_path
    swp['name'] = swp_name


    # check that json input file exists
    if not os.path.exists('Inputs/' + json_file):
        print('JSON input file does not exist')
        return

    with open('Inputs/' + json_file, 'r', encoding='utf-8') as f:  
        covis = json.load(f) 


    # Define global variable Verbose and read its value from the json file
    config.Verbose = covis['user']['verbose']

    # define a 2D rectangular data grid
    for n in range(len(covis['grid'])):
        covis['grid'][n] = covis_rectgrid(covis['grid'][n])    #  corr grid
        covis['grid'][n]['name'] = swp['name']


    # set local copies of covis structs
    pos = covis['sonar']['position']
    bfm = covis['processing']['beamformer']
    cal = covis['processing']['calibrate']
    filt = covis['processing']['filter']
    cor = covis['processing']['correlation']

    # Set the type of beamforming (fast, fft, ...)
    if not 'type' in bfm:
        bfm['type'] = 'fast'
        if config.Verbose:
            print('Setting beamform type: %s' % bfm['type'])

    # Set compass declination
    if not 'declination' in pos:
        pos['declination'] = 16.0

    # calibration parameters
    if not 'mode' in cal:
        cal['mode'] = 'TS-Wide'    # 'VSS', 'TS-Wide', 'TS-Fan'
    if not 'filt' in cal:
        cal['filt'] = 1    # 1 for low-pass filtering, 0 for no filtering
        cal['filt_bw'] = 2.0    # Bandwdith is filt_bw/tau

    # sonar height off bottom
    if not 'altitude' in pos:
        pos['altitude'] = 4.2

    # correlation range window size [number of samples]
    if not 'window_size' in cor:
        cor['window_size'] = 0.001
    window_size = cor['window_size']

    # correlation window overlap [number of samples]
    if not 'window_overlap' in cor:
        cor['window_overlap'] = 0.4
    window_overlap = cor['window_overlap']


    # directory list of *.bin file
    all_files = os.listdir(filedir)
    file = [f for f in all_files if f.endswith('.bin')]
    nfiles = len(file)

    # check if there are missing json or binary files
    json_file = [f for f in all_files if f.endswith('.json')]
    if len(json_file) < nfiles:
        print('Warning: there are %d json files missing in the sweep' % (nfiles - len(json_file)))
    elif len(json_file) > nfiles:
        print('Warning: there are %d bin files missing in the sweep' % (len(json_file) - nfiles))

    # read ping meta data from json file
    with open(filedir + '/' + json_file[0], 'r', encoding='utf-8') as f:  
        jsonf = json.load(f) 
    ind = 1
    while abs(jsonf['hdr']['sample_rate'] - 34482) > 1 and ind <= len(json_file):
        with open(filedir + '/' + json_file[ind], 'r', encoding='utf-8') as f:  
            jsonf = json.load(f) 
        ind = ind + 1

    # numbers of samples in correlation window and window overlap
    fsamp = jsonf['hdr']['sample_rate']
    cwsize = np.round(fsamp * window_size)
    cwovlap = np.round(window_overlap * cwsize)

    # Read index file
    # The index file contains the sweep parameters:
    # ping,seconds,microseconds,pitch,roll,yaw,kPAngle,kRAngle,kHeading
    # in csv format
    ind_file = 'index.csv'

    # should check if file exists
    if config.Verbose:
        print('Parsing ' + filedir + '/' + ind_file)
    csv = pd.read_csv(filedir + '/' + ind_file, encoding='utf-8')
    # remove the first dummy ping that logs sonar neutral orientations
    csv = csv.drop(csv[csv['ping'] == 0].index)
    # save index data in png structure
    png = []
    nping = 0
    for row in csv.itertuples():
        if getattr(row, 'yaw') == 0:
            rot_yaw = 0
        else:
            rot_yaw = getattr(row, 'yaw') - central_yaw
        png.append({
            'num' : getattr(row, 'ping'),
            'sec' : getattr(row, 'seconds') + getattr(row, 'usecs') / 1e6,
            'rot_pitch' : getattr(row, 'pitch'),
            'sen_pitch' : getattr(row, 'kPAngle'),
            'rot_roll' : getattr(row, 'roll') / 6,
            'sen_roll' : getattr(row, 'kRAngle'),
            'rot_yaw' : rot_yaw,
            'sen_head' : getattr(row, 'kHeading'),
            'hdr' : jsonf['hdr']
        })
        nping = nping + 1

    # sonar orientation during the sweep
    pitch = (np.pi/180) * png[0]['sen_pitch']
    roll = (np.pi/180) * png[0]['sen_roll']
    yaw = (np.pi/180) * png[0]['rot_yaw']


    # -------------------------
    # Main program
    # -------------------------

    # Loop over ping files
    bad_ping = []
    bad_ping_count = 0
    count = 0

    for npg in range(nping):
        ping_num = png[npg]['num']
        bin_file = 'rec_7038_%06d.bin' % ping_num
        if not os.path.exists(filedir + '/' + bin_file):
            print('Warning: binary file missing for ping:%d' % ping_num)
            continue

        if npg == 1:
            if config.Verbose > 1:
                print(png[npg]['hdr'])    # View essential parameters

        if config.Verbose > 1:
            print('Reading %s: pitch %f, roll %f, yaw %f' % bin_file, pitch*180/np.pi, roll*180/np.pi, yaw*180/np.pi)
        
        # read raw element quadrature data
        try:
            hdr, data1 = covis_read(filedir + '/' + bin_file)
        except:
            print('Warning: error reading ping %d' % ping_num)
            bad_ping_count = bad_ping_count + 1
            bad_ping[bad_ping_count] = ping_num
            continue

        if data1.shape[0] == 0 or data1.shape[1] != 128:
            print('Warning: error reading ping %d' % ping_num)
            bad_ping_count = bad_ping_count + 1
            bad_ping[bad_ping_count] = ping_num
            continue

        count = count + 1
        if count == 1:
            data = np.full((data1.shape[0],data1.shape[1],nping),np.nan,dtype=np.complex128)
        data[:,:,npg] = data1

    # End loop on pings


    # Loop again on pings to calculate backscatter intensity and
    # scintillation parameters
    I_t1 = np.full(data.shape,np.nan,dtype=np.complex128)
    Isq_t = np.full(data.shape,np.nan,dtype=np.complex128)
    bf_sig_t = np.full(data.shape,np.nan,dtype=np.complex128)
    data0 = data[:,:,0]

    for npg in range(data.shape[2]):
        if np.all(np.isnan(data0)):
            data0 = data[:,:,npg + 1]
            continue
        data1 = data[:,:,npg]
        # Correct phase
        try:
            data1 = covis_phase_correct(png[npg], data0, data1)
        except:
            print('Warning: error in phase correction for ping: %d' % png[npg]['num'])
            continue
        
        # define beamformer parameters
        bfm['fc'] = png[npg]['hdr']['xmit_freq']
        bfm['c'] = c
        bfm['fs'] = png[npg]['hdr']['sample_rate']
        bfm['first_samp'] = hdr['first_samp'] + 1
        bfm['last_samp'] = hdr['last_samp'] + 1
        bfm['start_angle'] = -64
        bfm['end_angle'] = 64

        # Apply Filter to data
        try:
            data1, filt, png[npg] = covis_filter(data1, filt, png[npg])
        except:
            print('Warning: error in filtering ping: %d' % png[npg]['num'])
            continue

        # beamform the quadrature data
        try:
            bfm, bf_sig1 = covis_beamform(bfm, data1)
        except:
            print('Warning: error in beamforming for ping: %d' % png[npg]['num'])
            continue

        # apply calibration to beamformed data
        try:
            bf_sig1 = covis_calibration(bf_sig1, bfm, png[npg], cal,T,S,pH,lat,depth)
        except:
            print('Warning: error in calibrating ping: %d' % png[npg]['num'])
            continue
        
        bf_sig_t[:,:,npg] = bf_sig1
    
        # calculate scintillation index and log-amplitude flluctuations
        I1 = np.power(np.abs(bf_sig1),2)
        Isq1 = np.power(np.abs(bf_sig1),4)
        I_t1[:,:,npg] = I1
        Isq_t[:,:,npg] = Isq1
    
    # average over the selected time window
    ping_sec = np.array([element['sec'] for element in png[0:nping]])
    time_win1 = ping_sec[0]
    time_win2 = time_win1 + avg_win
    ii_ave = np.where((ping_sec >= time_win1) & (ping_sec <= time_win2))[0]
    nave = len(ii_ave)    # number of samples used for averaging
    I_av1 = np.nanmean(I_t1[:,:,ii_ave],2)
    X = np.log(I_t1[:,:,ii_ave]/np.tile(I_av1[:,:,np.newaxis],(1,1,len(ii_ave))))
    J = np.power(np.abs(np.nanmean(bf_sig_t[:,:,ii_ave]/np.abs(bf_sig_t[:,:,ii_ave]),2)),2)
    X_var = np.nanvar(X[:,:,ii_ave],2)
    SI = nave*np.nansum(Isq_t[:,:,ii_ave],2)/np.power(np.nansum(I_t1[:,:,ii_ave],2),2)-1
    Kp = np.nanmean(np.power(np.abs(bf_sig_t[:,:,ii_ave]),2),2)/np.power(np.nanmean(np.abs(bf_sig_t[:,:,ii_ave]),2),2)-1
    sig_phi2 = np.log(1/J)
    sig_phi2[J==0] = np.nan

    # loop over pings again to form ping-ping decorrelation
    ping_rate = png[0]['hdr']['max_ping_rate']
    bf_sig_t_sub = bf_sig_t[:,:,ii_ave]
    ping_sec_sub = ping_sec[ii_ave]
    
    for npg in range(bf_sig_t_sub.shape[2]):
        # get range and azimuthal angles
        rng = bfm['range']
        azim = bfm['angle']

        # calculate the target strength corresponding to the noise floor
        if npg == 0:
            bf_sig1 = np.squeeze(bf_sig_t[:,:,npg])
            bf_sig_noise = noise_floor*np.ones(bf_sig1.shape)
            bf_sig_noise = covis_calibration(bf_sig_noise,bfm,png[len(png)-1],cal,T,S,pH,lat,depth)
            __, E1_noise, E2_noise, __ = covis_covar_hamming(bf_sig_noise, bf_sig_noise, rng, cwsize, cwovlap)
            I_noise1 = np.power(np.abs(bf_sig_noise),2)
            I_noise2 = np.sqrt(E1_noise*E2_noise)
            I_t2 = np.full((I_noise2.shape[0],I_noise2.shape[1],data.shape[2]),np.nan)
            cov_t = np.full(I_t2.shape,np.nan)
            E1_t = np.full(I_t2.shape,np.nan)
            E2_t = np.full(I_t2.shape,np.nan)

        t1 = png[npg]['sec']
        t2 = t1 + tlag
        if (np.abs(ping_sec_sub-t2) < 1/ping_rate).any():
            npg2 = np.argmin(np.abs(ping_sec_sub-t2))
        else:
            continue
        
        # pair of pings used
        bf_sig1 = bf_sig_t_sub[:,:,npg]
        bf_sig2 = bf_sig_t_sub[:,:,npg2]

        # Correlate pings
        #  rc is the range of the center of the corr bin
        try:
            cov,E1,E2,rc = covis_covar_hamming(bf_sig1, bf_sig2, rng, cwsize, cwovlap)
        except:
            print('Warning: error in calculating decorrelation between pings %d and %d' % png[npg]['num'],png[npg2]['num'])
            continue
        I2 = np.sqrt(E1*E2)
        
        # save the quanities of interest
        cov_t[:,:,npg] = cov
        E1_t[:,:,npg] = E1
        E2_t[:,:,npg] = E2
        I_t2[:,:,npg] = I2


    # average oveer ping pairs
    cor_av = np.abs(np.nansum(cov_t,2))/np.sqrt(np.nansum(E1_t,2)*np.nansum(E2_t,2))
    I_av2 = np.nanmean(I_t2,2)
    decor_I_av = (1-cor_av)*I_av2

    # mask out data points with low SNR
    snr1 = 10*np.log10(I_av1/I_noise1)
    snr2 = 10*np.log10(I_av2/I_noise2)
    mask1 = np.full(snr1.shape,np.nan)
    mask2 = np.full(snr2.shape,np.nan)
    mask1[snr1>snr_thresh] = 1
    mask2[snr2>snr_thresh] = 1
    cor_av = cor_av*mask2
    decor_I_av = decor_I_av*mask2
    X_var = X_var*mask1
    Kp = Kp*mask1
    SI = SI*mask1
    sig_phi2 = sig_phi2*mask1


    ## Find seafloor origins of acoustic backscatter

    # load bathymetry data
    alt = covis['sonar']['position']['altitude']    # height of COVIS
    covis_bathy = bathy['covis']
    bathy_grid = covis_bathy['grid'][0][0][0][0]
    x_bathy = bathy_grid['x']
    y_bathy = bathy_grid['y']
    z_bathy = bathy_grid['v'] - alt


    # Determine the real-world coordinates of pesudo planes perpendicular along each
    # beam and the intercept of the plane on the seafloor

    ver_beam = np.linspace(-40,40,161,endpoint=True)*np.pi/180    # vertical beam width
    
    # rotation matrices
    R = np.array([[np.cos(roll), 0, -np.sin(roll)],
                  [0, 1, 0],
                  [np.sin(roll), 0, np.cos(roll)]])
    
    Y = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                  [np.sin(yaw), np.cos(yaw), 0],
                  [0, 0, 1]])
    P = np.array([[1, 0, 0],
                  [0, np.cos(pitch), np.sin(pitch)],
                  [0, -np.sin(pitch), np.cos(pitch)]])

    # COVIS central pointing direction re true north
    ang = np.pi/180*(central_head + pos['declination'] - 360)
    ROT = np.array([[np.cos(ang), np.sin(ang), 0],
                    [-np.sin(ang), np.cos(ang), 0],
                    [0, 0, 1]])
    M = ROT@Y.T@P.T@R.T
    x_out1 = np.zeros((len(rng),len(azim)))
    y_out1 = np.zeros((len(rng),len(azim)))
    x_out2 = np.zeros((len(rc),len(azim)))
    y_out2 = np.zeros((len(rc),len(azim)))


    # Python的插值函数接收一维坐标，matlab插值函数接收网格化的坐标
    # 这里要将网格化坐标转为一维的坐标
    x_bathy1d = x_bathy[0,:]
    y_bathy1d = y_bathy[:,0]
    zw2_interp = scipy.interpolate.RegularGridInterpolator((x_bathy1d, y_bathy1d), z_bathy.T, bounds_error=False)

    for j in range(len(azim)):
        # receiver coordinates of a pseudo plane perpendicular along the beam
        azim1 = azim[j]
        xr = np.array([rng]).T@np.cos([ver_beam])*np.sin(azim1)
        yr = np.array([rng]).T@np.cos([ver_beam])*np.cos(azim1)
        zr = np.array([rng]).T@np.sin([ver_beam])
        rshape = xr.shape
        rr = np.array([xr.reshape(-1), yr.reshape(-1), zr.reshape(-1)])

        # real-world coordinates of a pseudo plane
        rw = M@rr

        xw1 = rw[0,:]
        yw1 = rw[1,:]
        zw1 = rw[2,:]


        # Find the intercept of the pseudo plane on the seafloor
        zw2 = zw2_interp(np.array([xw1, yw1]).T)
        dz = np.abs(zw1-zw2)
        ii = np.where(dz < 0.1)[0]
        xw_sub = xw1[ii]
        yw_sub = yw1[ii]
        zw2_sub = zw2[ii]
        L_sub = np.sqrt(np.power(xw_sub,2)+np.power(yw_sub,2)+np.power(zw2_sub,2))

        min_diff_rc = 0.5*np.min(np.diff(rc.reshape(-1)))
        for i in range(len(rng)):
            rng1 = rng[i]
            ii = np.where(np.abs(L_sub-rng1) < min_diff_rc)[0]
            if ii.shape[0] != 0:
                x_out1[i,j] = np.mean(xw_sub[ii])
                y_out1[i,j] = np.mean(yw_sub[ii])
            else:
                x_out1[i,j] = np.nan
                y_out1[i,j] = np.nan
        
        for i in range(len(rc)):
            rc1 = rc[i]
            ii = np.where(np.abs(L_sub-rc1) < min_diff_rc)[0]
            if ii.shape[0] != 0:
                x_out2[i,j] = np.mean(xw_sub[ii])
                y_out2[i,j] = np.mean(yw_sub[ii])
            else:
                x_out2[i,j] = np.nan
                y_out2[i,j] = np.nan
    

    # interpolate onto a 2D grid
    for k in range(len(covis['grid'])):
        grd = covis['grid'][k]
        if grd['type'] == 'decorrelation intensity':
                grd['v'], grd['w'] = l2grid(x_out2, y_out2, decor_I_av, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'decorrelation':
                grd['v'], grd['w'] = l2grid(x_out2, y_out2, 1-cor_av, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'intensity':
                grd['v'], grd['w'] = l2grid(x_out1, y_out1, I_av1, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'Chi_var':
                grd['v'], grd['w'] = l2grid(x_out1, y_out1, X_var, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'SI':
                grd['v'], grd['w'] = l2grid(x_out1, y_out1, SI, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'Sig_phi2':
                grd['v'], grd['w'] = l2grid(x_out1, y_out1, sig_phi2, grd['x'], grd['y'], grd['v'], grd['w'])
        elif grd['type'] == 'Kp':
                grd['v'], grd['w'] = l2grid(x_out1, y_out1, Kp, grd['x'], grd['y'], grd['v'], grd['w'])
        covis['grid'][k] = grd

    # normalize the grid with the grid weights
    for k in range(len(covis['grid'])):
        grd = covis['grid'][k]
        grd['v'] = np.divide(grd['v'], grd['w'], where=grd['w']!=0)
        covis['grid'][k] = grd


    # save metadata into the covis structure
    covis_averaging = {
        'win' : avg_win
    }
    covis_snr = {
        'noise_floor' : noise_floor,
        'threshold' : snr_thresh
    }
    covis_vers = covis_version()
    covis['release'] = covis_vers['version_number']
    covis['sweep'] = swp
    covis['ping'] = png
    covis['sonar']['position'] = pos
    covis['processing']['beamformer'] = bfm
    covis['processing']['calibrate'] = cal
    covis['processing']['filter'] = filt
    covis['processing']['correlation']['tlag'] = tlag
    covis['processing']['averaging'] = covis_averaging
    covis['processing']['snr'] = covis_snr
    covis['bad_ping'] = bad_ping

    return covis
        





