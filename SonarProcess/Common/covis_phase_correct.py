import numpy as np

def covis_phase_correct(ping, data1, data2):
    chn = 10    # of the reference channel

    sample_rate = ping['hdr']['sample_rate']
    pulse_width = ping['hdr']['pulse_width']

    # Number of samples to keep from monitor channel
    nkeep = round(5 * sample_rate * pulse_width)

    # The monitor channel is the last channel
    monitor1 = data1[0:nkeep,chn-1]
    monitor2 = data2[0:nkeep,chn-1]

    # Cross-correlate monitor signals to determine phase jitter
    phase_jitter = np.angle(np.sum(monitor2 * monitor1.conj()))

    # Correct phase of data2 to agree with data1
    data2_corr = np.exp(-1j * phase_jitter) * data2

    return data2_corr