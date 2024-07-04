import numpy as np
import struct

def covis_read(file):
    hdr = {}
    data = []
    element = []
    with open(file, "rb") as binary:
        context = binary.read(8)
        hdr['serial_number'] = struct.unpack('Q', context)[0]
        context = binary.read(4)
        hdr['number'] = struct.unpack('I', context)[0]
        reserved = binary.read(2)    # reserved
        context = binary.read(2)
        hdr['total_nelems'] = struct.unpack('H', context)[0]
        context = binary.read(4)
        hdr['nsamps'] = struct.unpack('I', context)[0]
        context = binary.read(2)
        hdr['nelems'] = struct.unpack('H', context)[0]
        context = binary.read(4)
        hdr['first_samp'] = struct.unpack('I', context)[0]
        context = binary.read(4)
        hdr['last_samp'] = struct.unpack('I', context)[0]
        context = binary.read(2)
        hdr['samp_type'] = struct.unpack('H', context)[0]
        reserved = binary.read(28)    # reserved

        nsamps = hdr['nsamps']
        nelems = hdr['nelems']

        # read element number list
        context = binary.read(2 * nelems)
        for elem in struct.unpack(str(nelems)+'H', context):
            element.append(elem)
        hdr['element'] = element

        # read data (I&Q pairs) for all elements 
        context = binary.read(2 * 2 * nelems * nsamps)
        for d in struct.unpack(str(2 * nelems * nsamps)+'h', context):
            data.append(d)
        
    # 注意：matlab的矩阵和Python的矩阵行列是相反的
    #data = np.array(data).reshape(2 * nelems, nsamps)
    data = np.array(data).reshape(nsamps, 2 * nelems)

    # If the sample type field is 16, the data value is stored in
    # the upper 12 bits rather than the lower 12 bits.
    if hdr['samp_type'] == 16:
        data = data / 16

    if(data.shape[0] < nsamps):
        print('Error reading %d samples' % nsamps)

    # make a complex number from real and imag components
    end = data.shape[1]
    data = data[:,0:end-1:2] + 1j * data[:,1:end:2]

    return hdr, data