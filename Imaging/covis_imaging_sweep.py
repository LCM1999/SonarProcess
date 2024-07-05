import numpy as np
import pandas as pd
import os
import gsw
import json
import joblib
from Common import config
from Imaging.covis_rectgrid_imaging import covis_rectgrid_imaging
from Common.covis_parse_bursts_pitch import covis_parse_bursts_pitch
from Common.covis_read import covis_read
from Common.covis_phase_correct import covis_phase_correct
from Common.covis_filter import covis_filter
from Common.covis_beamform import covis_beamform
from Common.covis_calibration import covis_calibration
from Common.covis_coords_darrell import covis_coords_darrell
from Imaging.l3grid_imaging import l3grid_imaging
from Common.covis_version import covis_version


def covis_imaging_sweep(filedir, json_file='Inputs/covis_image.json'):
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
    T = 2.41  # temperataure (degree C)
    S = 34.53  # salinity (PSU)
    pH = 7  # pH
    lat = 46  # latitude ( degree N )
    depth = 1544  # depth ( m )
    p = gsw.p_from_z(-depth, lat)
    c = gsw.sound_speed_t_exact(S, T, p)  # sound speed ( m/s )

    # snr parameters
    noise_floor = 0.1970  # rms noise floor
    snr_thresh = 50  # threshold (dB)

    # OSCFAR parameters
    clutterp = 65  # start with the 65% quantile as the clutter threshold
    scrthreshold = 20  # dB threshold for signal-to-clutter ratio

    if not os.path.exists(filedir):
        print('Error: imaging file directory not exists.')
        return

    # parse sweep.json file in data archive
    swp_file = 'sweep.json'
    full_swp_path = filedir + '/' + swp_file
    if not os.path.exists(full_swp_path):
        print('Error: could not find sweep.json in imaging file directory.')
        return
    with open(full_swp_path, 'r', encoding='utf-8') as f:
        swp = json.load(f)

        # save sweep path and name in swp structure
    swp['path'] = swp_path
    swp['name'] = swp_name

    # parsing the json input file for the user supplied parameters
    if not os.path.exists(json_file):
        print('Error: JSON input file does not exist.')
        return
    with open(json_file, 'r', encoding='utf-8') as f:
        covis = json.load(f)
    if covis['type'] != 'imaging':
        print('Error: Incorrect covis input file type. Looking for: imaging, Current type: ' + covis['type'])
        return
    print(covis)
    # Define global variable Verbose and read its value from the json file
    config.Verbose = covis['user']['verbose']
    print(covis)
    # define a 3D rectangular data grid
    covis['grid'] = covis_rectgrid_imaging(covis['grid'])
    covis['grid']['name'] = swp['name']
    grd_out = covis['grid']
    print(grd_out)

    # set local copies of covis structs
    grd = covis['grid']
    usr = covis['user']
    pos = covis['sonar']['position']
    dsp = covis['processing']
    bfm = covis['processing']['beamformer']
    cal = covis['processing']['calibrate']
    filt = covis['processing']['filter']

    # check the grid
    # Set the type of grid Ialue (intensity or complex).
    if not 'type' in grd:
        grd['type'] = 'intensity'  #default grid type
        if config.Verbose:
            print('Setting grid type: ' + grd['type'])

    # Set the type of beamforming (fast, fft, ...)
    if not 'type' in bfm:
        bfm['type'] = 'fast'
        if config.Verbose:
            print('Setting beamform type: ' + bfm['type'])

    # Set compass declination
    if not 'declination' in pos:
        pos['declination'] = 16.0

    # calibration parameters
    if not 'mode' in cal:
        cal['mode'] = 'VSS'  # 'VSS' or 'TS'

    # set position of sonar
    # should use covis.position info, but for now ...
    if not 'altitude' in pos:
        pos['altitude'] = 4.2
    alt = pos['altitude']
    origin = [0, 0, alt]  # sonar is always at (0,0,0) in world coords

    # directory list of *.bin file
    all_files = os.listdir(filedir)
    file = [f for f in all_files if f.endswith('.bin')]
    nfiles = len(file)

    # check if there are missing json files
    json_file = [f for f in all_files if f.endswith('.json')]
    if len(json_file) < nfiles:
        print('Warning: there are ' + str(nfiles - len(json_file)) + ' json files missing in the sweep')
    elif len(json_file) > nfiles:
        print('Warning: there are ' + str(len(json_file) - nfiles) + ' bin files missing in the sweep')

    # read ping meta data from json file
    with open(filedir + '/' + json_file[0], 'r', encoding='utf-8') as f:
        jsonf = json.load(f)
    ind = 1
    while abs(jsonf['hdr']['sample_rate'] - 34483) > 1 and ind <= len(json_file):
        with open(filedir + '/' + json_file[ind], 'r', encoding='utf-8') as f:
            jsonf = json.load(f)
        ind = ind + 1

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
    for row in csv.itertuples():
        if getattr(row, 'yaw') == 0:
            rot_yaw = 0
        else:
            rot_yaw = getattr(row, 'yaw') - central_yaw
        png.append({
            'num': getattr(row, 'ping'),
            'sec': getattr(row, 'seconds') + getattr(row, 'usecs') / 1e6,
            'rot_pitch': getattr(row, 'pitch'),
            'sen_pitch': getattr(row, 'kPAngle'),
            'rot_roll': getattr(row, 'roll') / 6,
            'sen_roll': getattr(row, 'kRAngle'),
            'rot_yaw': rot_yaw,
            'sen_head': getattr(row, 'kHeading'),
            'hdr': jsonf['hdr']
        })

    # group the pings transmitted at the same elevation angles into bursts
    burst = covis_parse_bursts_pitch(png)

    # range of elev angles to process
    elev_start = covis['processing']['bounds']['pitch']['start']
    elev_stop = covis['processing']['bounds']['pitch']['stop']

    # -------------------------
    # Main program
    # -------------------------

    # loop over bursts
    bad_ping = []
    bad_ping_count = 0
    nbursts = len(burst)
    burst_count = 0

    for nb in range(nbursts):
        # check elevation
        if (burst[nb]['pitch'] < elev_start) or (burst[nb]['pitch'] > elev_stop):
            continue
        if config.Verbose:
            print('Burst %d: pitch %f' % (nb, burst[nb]['pitch']))

        npings = burst[nb]['npings']

        # check that there's enough pings in burst
        if npings < 2:
            print('Warning: not enough pings in burst')
            continue

        # loop over pings in a burst
        ping_count = 0
        for npg in range(npings):
            ping_num = burst[nb]['ping'][npg]  # ping number
            # read the corresponding binary file
            bin_file = 'rec_7038_%06d.bin' % ping_num
            if not os.path.exists(filedir + '/' + bin_file):
                print('Warning: binary file missing for ping:%d' % ping_num)
                continue
            for i in range(len(png)):
                if png[i]['num'] == ping_num:
                    ip = i

            # define sonar orientation based on TCM readings
            pitch = (np.pi / 180) * png[ip]['sen_pitch']
            roll = (np.pi / 180) * png[ip]['sen_roll']
            yaw = (np.pi / 180) * png[ip]['rot_yaw']

            if config.Verbose > 1:
                print('Reading: ' + filedir + '/' + bin_file)
            if config.Verbose > 2:
                print('Pitch %f, roll %f, yaw %f\n' % pitch * 180 / np.pi, roll * 180 / np.pi, yaw * 180 / np.pi)

            # read raw element quadrature data
            try:
                hdr, data = covis_read(filedir + '/' + bin_file)
            except:
                print('Warning: error reading ping %d at pitch %f' % ping_num, burst[nb]['pitch'])
                bad_ping_count = bad_ping_count + 1
                bad_ping.append(ping_num)
                continue

            if data.size == 0 or data.shape[1] != 256:
                print('Warning: error reading ping %d at pitch %f' % ping_num, burst[nb]['pitch'])
                bad_ping_count = bad_ping_count + 1
                bad_ping.append(ping_num)
                continue

            ping_count = ping_count + 1
            if ping_count == 1:
                monitor = data
                bf_sig_out = np.full((data.shape[0], data.shape[1], npings), np.nan, dtype=np.complex128)
                I_out = np.full(bf_sig_out.shape, np.nan, dtype=np.complex128)
                Isq_out = np.full(bf_sig_out.shape, np.nan, dtype=np.complex128)

            # Correct phase using first ping as reference
            try:
                data = covis_phase_correct(png[ip], monitor, data)
            except:
                print('Warning: error in phase correction for ping %d at pitch %f' % ping_num, burst[nb]['pitch'])
                bad_ping_count = bad_ping_count + 1
                bad_ping.append(ping_num)
                continue

            # Apply Filter to data
            try:
                data, filt, png[len(png) - 1] = covis_filter(data, filt, png[len(png) - 1])
            except:
                print('Warning: error in filtering ping %d at pitch %f' % ping_num, burst[nb]['pitch'])
                bad_ping_count = bad_ping_count + 1
                bad_ping.append(ping_num)
                continue

            # define beamformer parameters
            bfm['fc'] = png[ip]['hdr']['xmit_freq']
            bfm['c'] = c
            bfm['fs'] = png[ip]['hdr']['sample_rate']
            bfm['first_samp'] = hdr['first_samp'] + 1
            bfm['last_samp'] = hdr['last_samp']
            bfm['start_angle'] = -54
            bfm['end_angle'] = 54

            # conduct beamforming
            try:
                bfm, bf_sig = covis_beamform(bfm, data)
            except:
                print('Warning: error in beamforming for ping %d at pitch %f' % ping_num, burst[nb]['pitch'])
                bad_ping_count = bad_ping_count + 1
                bad_ping.append(ping_num)
                continue

            bf_sig_out[:, :, npg] = bf_sig
            I_out[:, :, npg] = np.power(np.abs(bf_sig), 2)
            Isq_out[:, :, npg] = np.power(np.abs(bf_sig), 4)

        if np.all(np.isnan(bf_sig_out)):
            print('Warning: no valid pings at pitch %f' % burst[nb].pitch)
            continue

        burst_count = burst_count + 1
        if burst_count == 1:
            xv_out = np.full((bf_sig.shape[0], bf_sig.shape[1], nbursts), np.nan)
            yv_out = np.full(xv_out.shape, np.nan)
            zv_out = np.full(xv_out.shape, np.nan)
            Ia_out = np.full(xv_out.shape, np.nan)
            Id_out = np.full(xv_out.shape, np.nan)
            Ia_filt_out = np.full(xv_out.shape, np.nan)
            Id_filt_out = np.full(xv_out.shape, np.nan)
            Kp_out = np.full(xv_out.shape, np.nan)

        # calculate SI
        Kp = np.nanmean(np.power(np.abs(bf_sig_out), 2), 2) / np.power(np.nanmean(np.abs(bf_sig_out), 2), 2) - 1

        # ping average
        average = np.nanmean(bf_sig_out, 2)
        bf_sig_a = np.abs(average)
        # remove the average to enhance plume signals
        bf_sig_d = np.sqrt(np.nanmean(
            np.power(np.abs(bf_sig_out - np.tile(average[:, :, np.newaxis], (1, 1, bf_sig_out.shape[2]))), 2), 2))

        # calculate signal to noise ratio
        snr_a = 20 * np.log10(np.abs(bf_sig_a) / noise_floor)
        snr_d = 20 * np.log10(np.abs(bf_sig_d) / noise_floor)

        # OSCFAR detection section
        # magnitude squared of the beamformed output
        mag2 = bf_sig_d * np.conj(bf_sig_d)
        # determine the specified quantile at each range (time) step
        clutter = np.percentile(mag2, clutterp, 1)
        dBclutter = 10 * np.log10(np.array([clutter]).T) @ np.ones((1, len(bfm['angle'])))
        dBmag = 10 * np.log10(mag2)
        dBSCR = dBmag - dBclutter  # signal to clutter ratio (SCR)
        indexD_d = dBSCR > scrthreshold  # detected samples

        # magnitude squared of the beamformed output
        mag2 = bf_sig_a * np.conj(bf_sig_a)
        # determine the specified quantile at each range (time) step
        clutter = np.percentile(mag2, clutterp, 1)
        dBclutter = 10 * np.log10(np.array([clutter]).T) @ np.ones((1, len(bfm['angle'])))
        dBmag = 10 * np.log10(mag2)
        dBSCR = dBmag - dBclutter  # signal to clutter ratio (SCR)
        indexD_a = dBSCR > scrthreshold  # detected samples

        # calibration
        try:
            bf_sig_d_cal = covis_calibration(bf_sig_d, bfm, png[len(png) - 1], cal, T, S, pH, lat, depth)
            bf_sig_a_cal = covis_calibration(bf_sig_a, bfm, png[len(png) - 1], cal, T, S, pH, lat, depth)
        except:
            print('Warning: error in calibration at pitch %f', burst[nb]['pitch'])
            continue

        Id = np.power(np.abs(bf_sig_d_cal), 2)
        Ia = np.power(np.abs(bf_sig_a_cal), 2)
        Id_filt = Id.copy()
        Ia_filt = Ia.copy()

        # mask out data with low snr
        Id_filt[snr_d < snr_thresh] = 10 ** -9
        Id_filt[indexD_d == False] = 10 ** -9
        Id[snr_d < snr_thresh] = 10 ** -9
        Ia_filt[snr_a < snr_thresh] = 10 ** -9
        Ia_filt[indexD_a == False] = 10 ** -9
        Ia[snr_a < snr_thresh] = 10 ** -9
        Kp[snr_d < snr_thresh] = 0

        # transform sonar coords into world coords
        rng = bfm['range']
        azim = bfm['angle']

        xv, yv, zv = covis_coords_darrell(origin, rng, azim, yaw, roll, pitch, central_head, pos['declination'])

        xv_out[:, :, nb] = xv
        yv_out[:, :, nb] = yv
        zv_out[:, :, nb] = zv
        Ia_out[:, :, nb] = Ia
        Id_out[:, :, nb] = Id
        Ia_filt_out[:, :, nb] = Ia_filt
        Id_filt_out[:, :, nb] = Id_filt
        Kp_out[:, :, nb] = Kp

    grd_in = {
        'x': xv_out,  # .reshape(-1).astype('float32'),
        'y': yv_out,  # .reshape(-1).astype('float32'),
        'z': zv_out,  # .reshape(-1).astype('float32'),
        'Ia': Ia_out,
        'Id': Id_out,  # .reshape(-1).astype('float32'),
        'Ia_filt': Ia_filt_out,
        'Id_filt': Id_filt_out,
        'Kp': Kp_out
    }

    '''
    with open('data/grd_in.bin', "wb") as file:
        joblib.dump(grd_in,file)
    with open('data/grd_out.bin', "wb") as file:
        joblib.dump(grd_out,file)
    exit(0)
    '''

    grd_out = l3grid_imaging(grd_in, grd_out)

    # normalize the grid with the grid weights
    grd_out['Ia'] = np.divide(grd_out['Ia'], grd_out['w'], where=grd_out['w'] != 0)
    grd_out['Id'] = np.divide(grd_out['Id'], grd_out['w'], where=grd_out['w'] != 0)
    grd_out['Ia_filt'] = np.divide(grd_out['Ia_filt'], grd_out['w'], where=grd_out['w'] != 0)
    grd_out['Id_filt'] = np.divide(grd_out['Id_filt'], grd_out['w'], where=grd_out['w'] != 0)
    grd_out['Kp'] = np.divide(grd_out['Kp'], grd_out['w'], where=grd_out['w'] != 0)

    # save local copies of covis structs
    covis_vers = covis_version()
    covis_snr = {
        'noise_floor': noise_floor,
        'threshold': snr_thresh
    }
    covis_oscfar = {
        'clutterp': clutterp,
        'scrthreshold': scrthreshold
    }
    covis['release'] = covis_vers['version_number']
    covis['sweep'] = swp
    covis['grid'] = grd_out
    covis['ping'] = png
    covis['sonar']['position'] = pos
    covis['processing']['beamformer'] = bfm
    covis['processing']['calibrate'] = cal
    covis['processing']['filter'] = filt
    covis['processing']['snr'] = covis_snr
    covis['processing']['oscfar'] = covis_oscfar
    covis['burst'] = ['burst']
    covis['bad_ping'] = ['bad_ping']

    return covis
