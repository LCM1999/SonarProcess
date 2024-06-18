import os
import os.path
from io import open
import struct
import numpy as np
import vtk

serial_number = 0
number = 0
# reserved = 0
total_nelems = 0
nsamps = 0
nelems = 0
first_samp = 0
last_samp = 0
samp_type = 0

element = []
data = []


with open("F:/Scivis2024/COVIS-20180729T160020-diffuse/COVIS-20180729T160020-diffuse/rec_7038_000001.bin", "rb") as binary:
    context = binary.read(8)
    serial_number = struct.unpack('Q', context)[0]
    context = binary.read(4)
    number = struct.unpack('I', context)[0]
    context = binary.read(2)    #reserved
    # reserved = struct.unpack('H', context)[0]
    context = binary.read(2)
    total_nelems = struct.unpack('H', context)[0]
    context = binary.read(4)
    nsamps = struct.unpack('I', context)[0]
    context = binary.read(2)
    nelems = struct.unpack('H', context)[0]
    context = binary.read(4)
    first_samp = struct.unpack('I', context)[0]
    context = binary.read(4)
    last_samp = struct.unpack('I', context)[0]
    context = binary.read(2)
    samp_type = struct.unpack('H', context)[0]
    reserved = binary.read(28)  #reserved

    context = binary.read(2 * nelems)
    for elem in struct.unpack(str(nelems)+'H', context):
        element.append(elem)

    context = binary.read(2 * 2 * nelems * nsamps)
    for d in struct.unpack(str(2 * nelems * nsamps)+'h', context):
        data.append(d)

print(nelems)
print(total_nelems)
print(nsamps)
print(samp_type)
print(element)

# print(data)

data = np.array(data).reshape(2 * nelems, nsamps)
if samp_type == 16:
    data = data / 16

print(data)


