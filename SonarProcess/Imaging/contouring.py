import vtk


# return a vtkPolyData
def get_isosurface(grid, iso_value, attribute_name):
    grid.GetPointData().SetActiveAttribute(attribute_name, 0)
    surface = vtk.vtkContourFilter()
    surface.SetInputData(grid)
    surface.ComputeScalarsOff()
    surface.ComputeNormalsOff()
    surface.ComputeGradientsOff()
    surface.SetValue(0, iso_value)
    surface.Update()
    return surface.GetOutput()


def save_isosurface(surface, path_dir, file_name):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path_dir + file_name)
    writer.SetInputData(surface)
    writer.Write()
