import vtk


def get_isosurface(grid, iso_value, attribute_name):
    grid.GetOutputData().SetActiveAttribute(attribute_name, 0)
    surface = vtk.vtkContourFilter()
    surface.SetInputData(grid)
    surface.ComputeScalarsOff()
    surface.ComputeNormalsOff()
    surface.ComputeGradientsOff()
    surface.SetValue(0, iso_value)
    surface.Update()
    return surface.GetOutput()

