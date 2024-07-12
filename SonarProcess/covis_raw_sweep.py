import numpy as np
import os, shutil
import tarfile
from Imaging.covis_imaging_sweep import covis_imaging_sweep
from Diffuse.covis_diffuse_sweep import covis_diffuse_sweep


def covis_raw_sweep(filepath):
    # 解压文件
    root1 = os.path.splitext(filepath)
    root2 = os.path.splitext(root1[0])
    head_tail = os.path.split(root2[0])
    if root1[1] == '.gz':
        tar = tarfile.open(filepath, "r:gz")
        tar.extractall(head_tail[0])
        tar.close()
    else:
        print('Error: unrecognized data format!')
        return
    
    swp_name = head_tail[1]
    swp_type = swp_name[22:25]

    if swp_type == 'ima':
        print('imaging: ' + root2[0])
        # ---------------------
        # 此处进入imaging数据处理:
        covis = covis_imaging_sweep(root2[0])
        # ---------------------
    elif swp_type == 'dif':
        print('diffuse: ' + root2[0])
        # ---------------------
        # 此处进入diffuse数据处理:
        covis = covis_diffuse_sweep(root2[0])
        # ---------------------
    elif swp_type == 'bat':
        print('bathymetry' + root2[0])
        # ---------------------
        # 此处进入bathymetry数据处理:
        # covis = covis_bathy_sweep(root2[0])
        # ---------------------

    # 删除解压缩文件
    shutil.rmtree(root2[0])

    return covis

if __name__ == "__main__":
    covis = covis_imaging_sweep("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20111001T210757.973Z-IMAGING")

  
