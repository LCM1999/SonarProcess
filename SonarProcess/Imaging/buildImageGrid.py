import math

import numpy as np
import scipy
import vtk
import vtkmodules


def build_image_grid(array_dict, bounds, spacing, need_smooth=True):
    # points = [(x[i], y[i], z[i]) for i in range(len(x))]
    image = vtk.vtkImageData()
    image.SetDimensions(math.ceil((bounds[1] - bounds[0]) / spacing[0]) + 1,
                        math.ceil((bounds[3] - bounds[2]) / spacing[1]) + 1,
                        math.ceil((bounds[5] - bounds[4]) / spacing[2]) + 1)

    image.SetOrigin(bounds[0], bounds[2], bounds[4])
    image.SetSpacing(spacing[0], spacing[1], spacing[2])

    dims = image.GetDimensions()
    pointsNum = dims[0] * dims[1] * dims[2]
    size = pointsNum

    vtk_arrays = []
    input_arrays = []

    array_num = 0
    for array_name, array_value in array_dict.items():
        if need_smooth and len(array_value.shape) <= 3:
            array_value = scipy.ndimage.gaussian_filter(array_value, sigma=1.1, mode='nearest')
        if len(array_value.shape) < 3:
            temp_array = vtkmodules.util.numpy_support.numpy_to_vtk(array_value, deep=True, array_type=vtk.VTK_FLOAT)
            temp_array.SetName(array_name)
            temp_array.SetNumberOfComponents(array_value.shape[-1])
            vtk_arrays.append(temp_array)
        else:
            vtk_arrays.append(vtk.vtkDoubleArray())
            vtk_arrays[array_num].SetName(array_name)
            if len(array_value.shape) == 4:
                vtk_arrays[array_num].SetNumberOfComponents(array_value.shape[3])
                vtk_arrays[array_num].SetNumberOfTuples(size)
            else:
                vtk_arrays[array_num].SetNumberOfComponents(1)
                vtk_arrays[array_num].SetNumberOfTuples(size)
            index = 0
            for z in range(dims[2]):
                for y in range(dims[1]):
                    for x in range(dims[0]):
                        if len(array_value.shape) == 4:
                            vtk_arrays[-1].InsertTuple(index, array_value[x][y][z])
                        else:
                            vtk_arrays[-1].InsertTuple1(index, array_value[x][y][z])
                        index += 1
        array_num += 1

    #             # hash_id = neighborSearcher.get_cell_id([x, y, z])
    #             # assert hash_id >= 0
    #             # pList = neighborSearcher.get_in_cell_list(hash_id)
    #             # if len(pList) == 0:
    #             #     Ia_array.InsertNextTuple1(0.0)
    #             #     Id_array.InsertNextTuple1(0.0)
    #             #     Ia_filt_array.InsertNextTuple1(0.0)
    #             #     Id_filt_array.InsertNextTuple1(0.0)
    #             #     Kp_array.InsertNextTuple1(0.0)
    #             #     continue

    for array in vtk_arrays:
        image.GetPointData().AddArray(array)
    return image


def save_image(image, path_dir, file_name):
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInputData(image)
    writer.SetFileName(path_dir + file_name)
    # writer.SetFileTypeToBinary()
    writer.Write()
