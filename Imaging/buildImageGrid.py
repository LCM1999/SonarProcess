import math
import scipy
import vtk

def build_image_grid(x, y, z, grd_out, need_smooth=False):
    points = [(x[i], y[i], z[i]) for i in range(len(x))]
    image = vtk.vtkImageData()
    bounding = [grd_out['bounds']['xmin'], grd_out['bounds']['xmax'],
                grd_out['bounds']['ymin'], grd_out['bounds']['ymax'],
                grd_out['bounds']['zmin'], grd_out['bounds']['zmax']]
    spacing = [grd_out['spacing']['dx'], grd_out['spacing']['dy'], grd_out['spacing']['dz']]
    image.SetDimensions(math.ceil((bounding[1] - bounding[0]) / spacing[0]) + 1,
                        math.ceil((bounding[3] - bounding[2]) / spacing[1]) + 1,
                        math.ceil((bounding[5] - bounding[4]) / spacing[2]) + 1)

    image.SetOrigin(bounding[0], bounding[2], bounding[4])
    image.SetSpacing(spacing[0], spacing[1], spacing[2])

    dims = image.GetDimensions()
    pointsNum = dims[0] * dims[1] * dims[2]

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
                Ia_array.InsertTuple1(index, grd_out['Ia'][index])
                Id_array.InsertTuple1(index, grd_out['Id'][index])
                Ia_filt_array.InsertTuple1(index, grd_out['Ia_filt'][index])
                Id_filt_array.InsertTuple1(index, grd_out['Id_filt'][index])
                Kp_array.InsertTuple1(index, grd_out['Kp'][index])
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
