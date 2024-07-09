import math
import scipy
import vtk


def build_image_grid(grid_dict, need_smooth=True):
    # points = [(x[i], y[i], z[i]) for i in range(len(x))]
    image = vtk.vtkImageData()
    bounding = [grid_dict['bounds']['xmin'], grid_dict['bounds']['xmax'],
                grid_dict['bounds']['ymin'], grid_dict['bounds']['ymax'],
                grid_dict['bounds']['zmin'], grid_dict['bounds']['zmax']]
    spacing = [grid_dict['spacing']['dx'], grid_dict['spacing']['dy'], grid_dict['spacing']['dz']]
    image.SetDimensions(math.ceil((bounding[1] - bounding[0]) / spacing[0]) + 1,
                        math.ceil((bounding[3] - bounding[2]) / spacing[1]) + 1,
                        math.ceil((bounding[5] - bounding[4]) / spacing[2]) + 1)

    image.SetOrigin(bounding[0], bounding[2], bounding[4])
    image.SetSpacing(spacing[0], spacing[1], spacing[2])

    dims = image.GetDimensions()
    pointsNum = dims[0] * dims[1] * dims[2]

    Ia = grid_dict['Ia']
    Id = grid_dict['Id']
    Ia_filt = grid_dict['Ia_filt']
    Id_filt = grid_dict['Id_filt']
    Kp = grid_dict['Kp']
    if need_smooth:
        Ia = scipy.ndimage.gaussian_filter(Ia, sigma=1.1, mode='nearest')
        Id = scipy.ndimage.gaussian_filter(Id, sigma=1.1, mode='nearest')
        Ia_filt = scipy.ndimage.gaussian_filter(Ia_filt, sigma=1.1, mode='nearest')
        Id_filt = scipy.ndimage.gaussian_filter(Id_filt, sigma=1.1, mode='nearest')
        Kp = scipy.ndimage.gaussian_filter(Kp, sigma=1.1, mode='nearest')

    Ia_array = vtk.vtkFloatArray()
    Ia_array.SetName("Ia")
    Ia_array.SetNumberOfComponents(1)
    Ia_array.SetNumberOfTuples(pointsNum)

    Id_array = vtk.vtkFloatArray()
    Id_array.SetName("Id")
    Id_array.SetNumberOfComponents(1)
    Id_array.SetNumberOfTuples(pointsNum)

    Ia_filt_array = vtk.vtkFloatArray()
    Ia_filt_array.SetName("Ia_filt")
    Ia_filt_array.SetNumberOfComponents(1)
    Ia_filt_array.SetNumberOfTuples(pointsNum)

    Id_filt_array = vtk.vtkFloatArray()
    Id_filt_array.SetName("Id_filt")
    Id_filt_array.SetNumberOfComponents(1)
    Id_filt_array.SetNumberOfTuples(pointsNum)

    Kp_array = vtk.vtkFloatArray()
    Kp_array.SetName("Kp")
    Kp_array.SetNumberOfComponents(1)
    Kp_array.SetNumberOfTuples(pointsNum)

    index = 0
    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                # hash_id = neighborSearcher.get_cell_id([x, y, z])
                # assert hash_id >= 0
                # pList = neighborSearcher.get_in_cell_list(hash_id)
                # if len(pList) == 0:
                #     Ia_array.InsertNextTuple1(0.0)
                #     Id_array.InsertNextTuple1(0.0)
                #     Ia_filt_array.InsertNextTuple1(0.0)
                #     Id_filt_array.InsertNextTuple1(0.0)
                #     Kp_array.InsertNextTuple1(0.0)
                #     continue
                Ia_array.InsertTuple1(index, Ia[x][y][z])
                Id_array.InsertTuple1(index, Id[x][y][z])
                Ia_filt_array.InsertTuple1(index, Ia_filt[x][y][z])
                Id_filt_array.InsertTuple1(index, Id_filt[x][y][z])
                Kp_array.InsertTuple1(index, Kp[x][y][z])
                index += 1

    image.GetPointData().AddArray(Ia_array)
    image.GetPointData().AddArray(Id_array)
    image.GetPointData().AddArray(Ia_filt_array)
    image.GetPointData().AddArray(Id_filt_array)
    image.GetPointData().AddArray(Kp_array)

    return image


def save_image(image, path_dir, file_name):
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInputData(image)
    writer.SetFileName(path_dir + file_name)
    # writer.SetFileTypeToBinary()
    writer.Write()
