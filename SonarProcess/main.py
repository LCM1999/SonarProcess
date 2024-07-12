from covis_raw_sweep import covis_raw_sweep
from Common import config
from Imaging.covis_imaging_sweep import covis_imaging_sweep
from Common.covis_read import covis_read
from Imaging.l3grid_imaging import l3grid_imaging
import numpy as np
import joblib


if __name__ == "__main__":
    covis = covis_raw_sweep('F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20111001T210757.973Z-IMAGING.tar.gz')
    #with open('data/covis.bin', "rb") as file:
    #    covis = joblib.load(file)

    print(covis)








# test.txt
#covis = covis_imaging_sweep('raw/COVIS-20190719T040001-imaging1')
#with open('data/covis.bin', "wb") as file:
#    joblib.dump(covis,file)
#print(covis)
#hdr,data=covis_read('raw/COVIS-20190719T040001-imaging1/rec_7038_000016.bin')


#with open('data/grd_in.bin', "rb") as file:
#    grd_in = joblib.load(file)
#with open('data/grd_out.bin', "rb") as file:
#    grd_out = joblib.load(file)

#grd_out2 = l3grid_imaging(grd_in, grd_out)
#with open('data/grd_out2.bin', "wb") as file:
#    joblib.dump(grd_out2,file)
#print(grd_out2['Ia'])

#with open('data/covis.bin', "rb") as file:
#    covis = joblib.load(file)

#print(covis['grid']['Id_filt'][177][105][40])
