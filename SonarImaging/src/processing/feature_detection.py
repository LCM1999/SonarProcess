import math

import numpy
import numpy as np
import scipy.optimize
from sklearn.preprocessing import normalize
import vtk, vtkmodules
from vtkmodules.util import numpy_support
import scipy
import queue
# import config

from SonarImaging.src.processing.read_from_mat import *

# marks = []
# regions = []
th = 3e-5
lamb = 1.06
lamb2 = lamb * lamb
alpha = 0.1
alpha_T = 1.32e-4
ro_ref = 1.04e3
Cp = 3.92e3
g = 9.8


# def region_growing(grp, seed):
#     global marks
#     size = config.imaging['data_cut'].shape
#     que = queue.Queue()
#     que.put(seed)
#     marks[seed] = -1
#     region = {
#         'id': grp,
#         'count': 0,
#         'max_x': 0,
#         'max_y': 0,
#         'max_z': 0,
#         'min_x': size[0],
#         'min_y': size[1],
#         'min_z': size[2]
#     }
#
#     while not que.empty():
#         pnt = que.get()
#         if config.imaging['data_cut'][pnt] > th:
#             marks[pnt] = grp
#             region['count'] = region['count'] + 1
#             region['max_x'] = max(region['max_x'], pnt[0])
#             region['min_x'] = min(region['min_x'], pnt[0])
#             region['max_y'] = max(region['max_y'], pnt[1])
#             region['min_y'] = min(region['min_y'], pnt[1])
#             region['max_z'] = max(region['max_z'], pnt[2])
#             region['min_z'] = min(region['min_z'], pnt[2])
#         else:
#             continue
#         dx = [1, -1, 0, 0, 0, 0]
#         dy = [0, 0, 1, -1, 0, 0]
#         dz = [0, 0, 0, 0, 1, -1]
#         max_x = size[0] - 1
#         max_y = size[1] - 1
#         max_z = size[2] - 1
#         for ii in range(6):
#             x = pnt[0] + dx[ii]
#             y = pnt[1] + dy[ii]
#             z = pnt[2] + dz[ii]
#             if 0 <= x <= max_x and 0 <= y <= max_y and 0 <= z <= max_z and marks[x][y][z] == 0:
#                 que.put((x, y, z))
#                 marks[x][y][z] = -1
#     return region
#
#
# def caculate_region(interval=4):
#     global marks, regions
#     regions = []
#     size = config.imaging['data_cut'].shape
#     marks = np.zeros(size, dtype=np.int32)
#     for i in range(0, size[0], interval):
#         for j in range(0, size[1], interval):
#             for k in range(0, size[2], interval):
#                 if (marks[i][j][k] == 0):
#                     region = region_growing(len(regions) + 1, (i, j, k))
#                     if region['count'] > 0:
#                         regions.append(region)
#
#
# def get_center_points(region_id, interval=1):
#     global marks, regions
#     points = []
#     region = regions[region_id - 1]
#     min_x = region['min_x']
#     max_x = region['max_x']
#     min_y = region['min_y']
#     max_y = region['max_y']
#     min_z = region['min_z']
#     max_z = region['max_z']
#     bounding_data = config.imaging['data_cut'][min_x:max_x + 1, min_y:max_y + 1, min_z:max_z + 1]
#     bounding_marks = marks[min_x:max_x + 1, min_y:max_y + 1, min_z:max_z + 1]
#     for z in range(0, bounding_data.shape[2], interval):
#         plane_data = bounding_data[:, :, z]
#         plane_marks = bounding_marks[:, :, z]
#         plane_data[plane_marks != region_id] = 0
#         max_v = np.max(plane_data)
#         max_points = np.where(plane_data == max_v)
#         points.append((max_points[0][0] + min_x, max_points[1][0] + min_y, z + min_z))
#     return points


def curve_fitting(points):
    npoints = np.array(points)
    x = npoints[:, 0]
    y = npoints[:, 1]
    z = npoints[:, 2]

    def curve_func(t, a0, a1, a2):
        return a0 + a1 * t + a2 * t ** 2

    p0 = [1, 1, 1]
    params_x, pcov_x = scipy.optimize.curve_fit(curve_func, z, x, p0)
    params_y, pcov_y = scipy.optimize.curve_fit(curve_func, z, y, p0)

    xx = curve_func(z, *params_x)
    yy = curve_func(z, *params_y)

    new_points = []
    for i in range(len(points)):
        new_points.append((xx[i], yy[i], z[i]))

    return new_points


def temp_region_growing(grid_in, scalar_name, seed, marks, region_id):
    shape = grid_in['shape']
    que = queue.Queue()
    que.put(seed)
    marks[seed] = -1

    region = {
        'id': region_id,
        'count': 0,
        'bounds': [shape[0], 0, shape[1], 0, shape[2], 0]
    }

    while not que.empty():
        pnt = que.get()
        if grid_in[scalar_name][pnt] > th:
            marks[pnt] = region_id
            region['count'] = region['count'] + 1
            region['bounds'][0] = min(region['bounds'][0], pnt[0])
            region['bounds'][1] = max(region['bounds'][1], pnt[0])
            region['bounds'][2] = min(region['bounds'][2], pnt[1])
            region['bounds'][3] = max(region['bounds'][3], pnt[1])
            region['bounds'][4] = min(region['bounds'][4], pnt[2])
            region['bounds'][5] = max(region['bounds'][5], pnt[2])
        else:
            continue
        dx = [1, -1, 0, 0, 0, 0]
        dy = [0, 0, 1, -1, 0, 0]
        dz = [0, 0, 0, 0, 1, -1]
        max_x = shape[0] - 1
        max_y = shape[1] - 1
        max_z = shape[2] - 1
        for ii in range(6):
            x = pnt[0] + dx[ii]
            y = pnt[1] + dy[ii]
            z = pnt[2] + dz[ii]
            if 0 <= x <= max_x and 0 <= y <= max_y and 0 <= z <= max_z and marks[x][y][z] == 0:
                que.put((x, y, z))
                marks[x][y][z] = -1
    return region


def get_plume_region(grid_in, array_name, threshold=100, interval=4):
    shape = grid_in['shape']
    regions = []
    plume_region_id = 0
    max_count = 0
    marks = np.zeros(shape, dtype=np.int32)
    for i in range(0, shape[0], interval):
        for j in range(0, shape[1], interval):
            for k in range(0, shape[2], interval):
                if marks[i][j][k] == 0:
                    temp_region = temp_region_growing(grid_in, array_name, (i, j, k), marks, len(regions) + 1)
                    if temp_region['count'] > 0:
                        regions.append(temp_region)
                        if temp_region['count'] > max_count:
                            max_count = temp_region['count']
                            plume_region_id = temp_region['id']
    return regions[plume_region_id - 1], marks


def get_plume_center_line(grid_in, array_name, marks, region: dict, interval=1):
    points = []
    grid_bounds = grid_in['bounds']
    region_bounds = region['bounds']
    # cut_box = vtk.vtkBox()
    # cut_box.SetXMin(bounds[0], bounds[2], bounds[3])
    # cut_box.SetXMax(bounds[1], bounds[3], bounds[5])
    # cutter = vtk.vtkClipDataSet()
    # cutter.SetClipFunction(cut_box)
    # cutter.SetInputData(vtk_dataset)
    # cutter.Update()
    # cut_dataset = cutter.GetOutput()
    # cut_shape = cutter.GetDimensions()
    bounding_data = grid_in[array_name][region_bounds[0]: region_bounds[1] + 1, region_bounds[2]: region_bounds[3] + 1,
                    region_bounds[4]: region_bounds[5] + 1]
    bounding_marks = marks[region_bounds[0]:region_bounds[1] + 1, region_bounds[2]:region_bounds[3] + 1,
                     region_bounds[4]:region_bounds[5] + 1]
    for z in range(0, bounding_data.shape[2], interval):
        plane_data = bounding_data[:, :, z]
        plane_marks = bounding_marks[:, :, z]
        plane_data[plane_marks != region['id']] = 0
        max_v = np.max(plane_data)
        max_points = np.where(plane_data == max_v)
        points.append((grid_bounds[0] + (region_bounds[0] + max_points[0][0]) * grid_in['spacing'][0],
                       grid_bounds[2] + (region_bounds[2] + max_points[1][0]) * grid_in['spacing'][1],
                       grid_bounds[4] + (region_bounds[4] + z) * grid_in['spacing'][2]))

    npoints = np.array(points)
    x = npoints[:, 0]
    y = npoints[:, 1]
    z = npoints[:, 2]

    def curve_func(t, a0, a1, a2):
        return a0 + a1 * t + a2 * t ** 2

    p0 = [1, 1, 1]
    params_x, pcov_x = scipy.optimize.curve_fit(curve_func, z, x, p0)
    params_y, pcov_y = scipy.optimize.curve_fit(curve_func, z, y, p0)

    xx = curve_func(z, *params_x)
    yy = curve_func(z, *params_y)

    new_points = []
    for i in range(len(points)):
        new_points.append(np.array([min(max(xx[i], grid_bounds[0]), grid_bounds[1]),
                                    min(max(yy[i], grid_bounds[2]), grid_bounds[3]),
                                    min(max(z[i], grid_bounds[4]), grid_bounds[5])]))

    return np.array(new_points)


def interpulate_sight_velocity(doppler, center_line: dict, grid_bounds, grid_spacing):
    npoints = np.array(center_line['center_points'])

    i = np.array([np.floor((npoints[:, 0] - grid_bounds[0]) / grid_spacing[0]),
                  np.ceil((npoints[:, 0] - grid_bounds[0]) / grid_spacing[0])]).astype(np.int32)
    j = np.array([np.floor((npoints[:, 1] - grid_bounds[2]) / grid_spacing[1]),
                  np.ceil((npoints[:, 1] - grid_bounds[2]) / grid_spacing[1])]).astype(np.int32)
    k = np.array([np.floor((npoints[:, 2] - grid_bounds[4]) / grid_spacing[2]),
                  np.ceil((npoints[:, 2] - grid_bounds[4]) / grid_spacing[2])]).astype(np.int32)
    wx = np.zeros(shape=i.shape)
    wx[1, :] = np.abs((i[0, :] * grid_spacing[0] + grid_bounds[0]) - npoints[:, 0]) / grid_spacing[0]
    wx[0, :] = 1 - wx[1, :]
    wy = np.zeros(shape=i.shape)
    wy[1, :] = np.abs((j[0, :] * grid_spacing[1] + grid_bounds[2]) - npoints[:, 1]) / grid_spacing[1]
    wy[0, :] = 1 - wy[1, :]
    wz = np.zeros(shape=i.shape)
    wz[1, :] = np.abs((k[0, :] * grid_spacing[2] + grid_bounds[4]) - npoints[:, 2]) / grid_spacing[2]
    wz[0, :] = 1 - wz[1, :]
    points_v = np.zeros(len(points))
    for n in range(2):
        for m in range(2):
            for l in range(2):
                w = wx[l, :] * wy[m, :] * wz[n, :]
                points_v = points_v + doppler[i[l, :], j[m, :], k[n, :]] * w

    center_line['sight_velocity'] = points_v
    return
    # return points_v


def calculate_center_velocity(center_line: dict, origin=(0.0, 0.0, 0.0)):
    size = center_line['size']
    points = center_line['center_points']
    center_velocity = np.zeros(size)

    tan_vec = np.zeros((size, 3))
    tan_vec[0] = points[1] - points[0]
    tan_vec[-1] = points[size - 1] - points[size - 2]
    tan_vec[1: size - 1] = ((points[1: size - 1] - points[0: size - 2]) + (points[2:] - points[1:size - 1])) / 2

    sight_vec = np.array(points - origin)

    unit_tan_vec = normalize(tan_vec, axis=1, norm='l2')
    unit_sight_vec = normalize(sight_vec, axis=1, norm='l2')

    center_velocity = center_line['sight_velocity'] / np.sum(unit_tan_vec * unit_sight_vec, axis=1)
    center_line['center_velocity'] = center_velocity
    center_line['tan_vec'] = unit_tan_vec
    center_line['sight_vec'] = unit_sight_vec

    return


def calculate_velocity_feild(doppler_grid, marks, plume_region, center_line: dict):
    shape = doppler_grid['shape']
    velocity_field = np.zeros((shape + [3]))
    grid_bounds = doppler_grid['bounds']
    region_bounds = plume_region['bounds']
    spacing = doppler_grid['spacing']

    size = center_line['size']
    points = center_line['center_points']
    center_v = center_line['center_velocity']
    center_tan_vec = center_line['tan_vec']

    for k in range(size):
        z = points[k][2]
        zi = k + region_bounds[4]
        v_vec = center_v[k] * center_tan_vec[k]
        vh = np.array([v_vec[0], v_vec[1], 0.0])
        for i in range(region_bounds[0], region_bounds[1]):
            for j in range(region_bounds[2], region_bounds[3]):
                if marks[i][j][zi] != plume_region['id']:
                    continue
                vp = doppler_grid['v'][i][j][zi]
                p = np.array([grid_bounds[0] + i * spacing[0],
                              grid_bounds[1] + j * spacing[1],
                              z])
                unit_p = p / np.linalg.norm(p)
                wp = (vp - np.dot(vh, unit_p)) / unit_p[2]
                velocity_field[i][j][zi] = [v_vec[0], v_vec[1], wp]

    return velocity_field


def Q_M_estimate(w_field, plume_region, center_line, spacing = 0.5):
    region_bounds = plume_region['bounds']

    size = center_line['size']
    points = center_line['center_points']
    center_v = center_line['center_velocity']
    center_tan_vec = center_line['tan_vec']

    Q = np.zeros(size)
    M = np.zeros(size)

    center_w = np.einsum('i, ij->ij', center_v, center_tan_vec)[:, 2]

    for k in range(size):
        plane_w_field = w_field[:, :, region_bounds[4] + k]
        # plane_marks = marks[:, :, region_bounds[4] + size]
        plane_w_field[plane_w_field < (center_w[k] * 0.1)] = 0
        Q[k] = np.sum(plane_w_field * np.power(spacing, 2)) / 0.9
        M[k] = np.sum(np.power(plane_w_field, 2) * np.power(spacing, 2)) / 0.99

    return Q, M


def be_estimate(Q, M):
    return Q / np.sqrt(2 * np.pi * M)


def Zi_estimate(be):
    return (5 * be) / (6 * alpha)


def B0_estimate(Q, Zi):
    return np.power(Q, 3) / ((3 * np.pi * (1 + lamb2)) / (2 * np.power(5 / (6 * alpha), 4)) * np.power(Zi, 5))


def H0_estimate(B0):
    return (Cp * ro_ref) / (g * alpha_T) * B0


def get_H0(v_field, plume_region, center_line, spacing=0.5):
    w_field = v_field[:, :, :, 2]
    Q, M = Q_M_estimate(w_field, plume_region, center_line, spacing)
    be = be_estimate(Q, M)
    Zi = Zi_estimate(be)
    B0 = B0_estimate(Q, Zi)
    H0 = H0_estimate(B0)
    return H0


def calculate_H_field(v_field, spacing=0.5):
    w_field = v_field[:, :, :, 2]
    dS = np.power(spacing, 2)
    Q_field = w_field * dS
    # M_field = np.power(w_field, 2) * dS
    be = 1 / np.sqrt(2 * np.pi)
    Zi = (5 * be) / (6 * alpha)
    B_field = B0_estimate(Q_field, Zi)
    H_field = H0_estimate(B_field)
    return H_field



if __name__ == "__main__":
    doppler = get_doppler_data("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20141006T210836.397Z-DOPPLER.mat")
    imaging = get_doppler_image("F:/Covis/pythonProject/Inputs/APLUWCOVISMBSONAR001_20141006T210836.397Z-DOPPLER.mat")
    th = 1e-7
    plume_region, marks = get_plume_region(imaging, 'v')
    points = get_plume_center_line(imaging, 'v', marks, plume_region)
    center_line = {
        'center_points': points,
        'size': len(points)
    }
    interpulate_sight_velocity(doppler['v'], center_line, imaging['bounds'], imaging['spacing'])
    calculate_center_velocity(center_line)

    v_field = calculate_velocity_feild(doppler, marks, plume_region, center_line)

    H = get_H0(v_field, plume_region, center_line)
    H_field = calculate_H_field(v_field)

    # v_array = vtkmodules.util.numpy_support.numpy_to_vtk(v_field.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    # v_array.SetName('velocity')
    # v_array.SetNumberOfComponents(3)

    velocity_image = build_image_grid(
        array_dict={'v': v_field, 'H': H_field}, bounds=doppler['bounds'], spacing=doppler['spacing'])
    save_image(velocity_image, "F:/data/", 'velocity.vti')

    sight_v_array = vtkmodules.util.numpy_support.numpy_to_vtk(center_line['sight_velocity'].ravel(), deep=True,
                                                               array_type=vtk.VTK_DOUBLE)
    sight_v_array.SetName('sight_v')
    sight_v_array.SetNumberOfComponents(1)

    center_v_array = vtkmodules.util.numpy_support.numpy_to_vtk(center_line['center_velocity'].ravel(), deep=True,
                                                                array_type=vtk.VTK_DOUBLE)
    center_v_array.SetName('center_v')
    center_v_array.SetNumberOfComponents(1)

    H0_array = vtkmodules.util.numpy_support.numpy_to_vtk(H.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    H0_array.SetName('H0')
    H0_array.SetNumberOfComponents(1)

    vtkPoints = vtk.vtkPoints()
    for point in points:
        vtkPoints.InsertNextPoint(point)
    polyLine = vtk.vtkPolyLine()
    polyLine.GetPointIds().SetNumberOfIds(vtkPoints.GetNumberOfPoints())
    for i in range(vtkPoints.GetNumberOfPoints()):
        polyLine.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyLine)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(vtkPoints)
    polyData.SetLines(cells)
    polyData.GetPointData().AddArray(sight_v_array)
    polyData.GetPointData().AddArray(center_v_array)
    polyData.GetPointData().AddArray(H0_array)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(polyData)
    writer.SetFileName("F:/data/center.vtp")
    writer.Write()
