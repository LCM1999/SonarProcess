import scipy

import vtk


def build_structured_grid(xv, yv, zv, Ia=None, Id=None, Ia_filt=None, Id_filt=None, Kp=None, need_smooth=False):
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(xv.shape)

    points = vtk.vtkPoints()
    points.Allocate(xv.size)

    Ia_array = None
    Id_array = None
    Ia_filt_array = None
    Id_filt_array = None
    Kp_array = None

    if Ia is not None:
        if need_smooth:
            Ia = scipy.ndimage.gaussian_filter(Ia, sigma=1.1, mode='nearest')
        Ia_array = vtk.vtkFloatArray()
        Ia_array.SetName('Ia')
        Ia_array.SetNumberOfComponents(1)
        Ia_array.SetNumberOfTuples(Ia.size)

    if Id is not None:
        if need_smooth:
            Id = scipy.ndimage.gaussian_filter(Id, sigma=1.1, mode='nearest')
        Id_array = vtk.vtkFloatArray()
        Id_array.SetName('Id')
        Id_array.SetNumberOfComponents(1)
        Id_array.SetNumberOfTuples(Id.size)

    if Ia_filt is not None:
        if need_smooth:
            Ia_filt = scipy.ndimage.gaussian_filter(Ia_filt, sigma=1.1, mode='nearest')
        Ia_filt_array = vtk.vtkFloatArray()
        Ia_filt_array.SetName('Ia_filt')
        Ia_filt_array.SetNumberOfComponents(1)
        Ia_filt_array.SetNumberOfTuples(Ia_filt.size)

    if Id_filt is not None:
        if need_smooth:
            Id_filt = scipy.ndimage.gaussian_filter(Id_filt, sigma=1.1, mode='nearest')
        Id_filt_array = vtk.vtkFloatArray()
        Id_filt_array.SetName('Id_filt')
        Id_filt_array.SetNumberOfComponents(1)
        Id_filt_array.SetNumberOfTuples(Id_filt.size)

    if Kp is not None:
        if need_smooth:
            Kp = scipy.ndimage.gaussian_filter(Kp, sigma=1.1, mode='nearest')
        Kp_array = vtk.vtkFloatArray()
        Kp_array.SetName('Kp')
        Kp_array.SetNumberOfComponents(1)
        Kp_array.SetNumberOfTuples(Kp.size)

    index = 0
    for z in range(xv.shape[2]):
        for y in range(xv.shape[1]):
            for x in range(xv.shape[0]):
                points.InsertPoint(index, (xv[x][y][z], yv[x][y][z], zv[x][y][z]))
                if Ia_array is not None:
                    Ia_array.InsertTuple1(index, Ia[x][y][z])
                if Id_array is not None:
                    Id_array.InsertTuple1(index, Id[x][y][z])
                if Ia_filt_array is not None:
                    Ia_filt_array.InsertTuple1(index, Ia_filt[x][y][z])
                if Id_filt_array is not None:
                    Id_filt_array.InsertTuple1(index, Id_filt[x][y][z])
                if Kp_array is not None:
                    Kp_array.InsertTuple1(index, Kp[x][y][z])
                index += 1

    sgrid.SetPoints(points)
    sgrid.GetPointData().SetVectors(Id_array)
    sgrid.GetPointData().AddArray(Ia_array)
    sgrid.GetPointData().AddArray(Id_array)
    sgrid.GetPointData().AddArray(Ia_filt_array)
    sgrid.GetPointData().AddArray(Id_filt_array)
    sgrid.GetPointData().AddArray(Kp_array)

    return sgrid


def save_structured_grid(sgrid, path_dir, file_name, need_binary=True):
    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(path_dir + file_name)
    writer.SetInputData(sgrid)
    if need_binary:
        writer.SetFileTypeToBinary()
    writer.Write()
