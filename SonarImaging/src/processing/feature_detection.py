import numpy as np
import vtk, vtkmodules
from vtkmodules.util import numpy_support
import scipy
import queue
import config


marks = []
regions = []
th = 3e-5



def region_growing(grp, seed):
    global marks
    size = config.imaging['data_cut'].shape
    que = queue.Queue()
    que.put(seed)
    marks[seed] = -1
    region = {
        'id': grp,
        'count': 0,
        'max_x' : 0,
        'max_y' : 0,
        'max_z' : 0,
        'min_x' : size[0],
        'min_y' : size[1],
        'min_z' : size[2]
    }

    while not que.empty():
        pnt = que.get()
        if config.imaging['data_cut'][pnt] > th:
            marks[pnt] = grp
            region['count'] = region['count'] + 1
            region['max_x'] = max(region['max_x'],pnt[0])
            region['min_x'] = min(region['min_x'],pnt[0])
            region['max_y'] = max(region['max_y'],pnt[1])
            region['min_y'] = min(region['min_y'],pnt[1])
            region['max_z'] = max(region['max_z'],pnt[2])
            region['min_z'] = min(region['min_z'],pnt[2])
        else:
            continue
        dx = [1,-1,0,0,0,0]
        dy = [0,0,1,-1,0,0]
        dz = [0,0,0,0,1,-1]
        max_x = size[0] - 1
        max_y = size[1] - 1
        max_z = size[2] - 1
        for ii in range(6):
            x = pnt[0] + dx[ii]
            y = pnt[1] + dy[ii]
            z = pnt[2] + dz[ii]
            if 0<=x<=max_x and 0<=y<=max_y and 0<=z<=max_z and marks[x][y][z]==0:
                que.put((x,y,z))
                marks[x][y][z] = -1
    return region


def caculate_region(interval=4):
    global marks, regions
    regions = []
    size = config.imaging['data_cut'].shape
    marks = np.zeros(size,dtype=np.int32)
    for i in range(0,size[0],interval):
        for j in range(0,size[1],interval):
            for k in range(0,size[2],interval):
                if(marks[i][j][k] == 0):
                    region = region_growing(len(regions)+1,(i,j,k))
                    if region['count'] > 0:
                        regions.append(region)


def get_center_points(region_id, interval=1):
    global marks, regions
    points = []
    region = regions[region_id-1]
    min_x = region['min_x']
    max_x = region['max_x']
    min_y = region['min_y']
    max_y = region['max_y']
    min_z = region['min_z']
    max_z = region['max_z']
    bounding_data = config.imaging['data_cut'][min_x:max_x+1, min_y:max_y+1, min_z:max_z+1]
    bounding_marks = marks[min_x:max_x+1, min_y:max_y+1, min_z:max_z+1]
    for z in range(0,bounding_data.shape[2],interval):
        plane_data = bounding_data[:,:,z]
        plane_marks = bounding_marks[:,:,z]
        plane_data[plane_marks!=region_id] = 0
        max_v = np.max(plane_data)
        max_points = np.where(plane_data==max_v)
        points.append((max_points[0][0]+min_x,max_points[1][0]+min_y,z+min_z))
    return points

                        
