import scipy

import vtk


def build_structured_grid(xv, yv, zv, array_dict, shape, need_smooth=True):
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(shape)

    size = shape[0] * shape[1] * shape[2]
    points = vtk.vtkPoints()
    points.Allocate(size)

    vtk_arrays = []
    input_arrays = []

    array_num = 0
    for array_name, array_value in array_dict.items():
        temp_array = array_value
        if need_smooth and len(temp_array.shape) <= 3:
            temp_array = scipy.ndimage.gaussian_filter(array_value, sigma=1.1, mode='nearest')
        vtk_arrays.append(vtk.vtkDoubleArray())
        input_arrays.append(temp_array)
        vtk_arrays[array_num].SetName(array_name)
        if len(temp_array.shape) == 4:
            vtk_arrays[array_num].SetNumberOfComponents(temp_array.shape[3])
            vtk_arrays[array_num].SetNumberOfTuples(size * temp_array.shape[3])
        else:
            vtk_arrays[array_num].SetNumberOfComponents(1)
            vtk_arrays[array_num].SetNumberOfTuples(size)

    index = 0
    for z in range(shape[2]):
        for y in range(shape[1]):
            for x in range(shape[0]):
                points.InsertPoint(index, (xv[x][y][z], yv[x][y][z], zv[x][y][z]))
                for i in range(len(vtk_arrays)):
                    vtk_arrays[i].InsertTuple1(index, input_arrays[i][x][y][z])
                index += 1

    sgrid.SetPoints(points)
    for array in vtk_arrays:
        sgrid.GetPointData().AddArray(array)
    return sgrid


def save_structured_grid(sgrid, path_dir, file_name, need_binary=True):
    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(path_dir + file_name)
    writer.SetInputData(sgrid)
    if need_binary:
        writer.SetFileTypeToBinary()
    writer.Write()
