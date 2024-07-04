import numpy as np
from Common import config

def covis_parse_bursts_pitch(png):
    pitch_rot = [d['rot_pitch'] for d in png]
    pitch_tcm = [d['sen_pitch'] for d in png]
    ping_no = [d['num'] for d in png]
    a, ia = np.unique(pitch_rot, return_index=True)
    ia.sort()
    burst = []
    for i in range(len(ia)):
        if i != len(ia) - 1:
            bping = ping_no[ia[i] : ia[i+1]]
            bnpings = ia[i+1] - ia[i]
        else:
            bping = ping_no[ia[i] : len(ping_no)]
            bnpings = len(ping_no) - ia[i]
        burst.append({
            'pitch' : pitch_tcm[ia[i]],
            'ping' : bping,
            'npings' : bnpings
        })
        if config.Verbose > 1:
            print('Burst %d: pitch=%f, npings=%d, start_ping=%d' % (i, burst[i]['pitch'], burst[i]['npings'], burst[i]['ping'][0]))
    return burst