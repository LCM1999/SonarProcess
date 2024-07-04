import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# render
render = vtk.vtkRenderer()

# interactor
interactor = []

# data
imaging = {
    'data' : [],
    'bounds' : [],
    'spacing' : [],
    'data_cut' : [],
    'cut_bounds': [],
    'vtk_data' : []
}

bathy = {
    'data' : [],
    'bounds' : [],
    'spacing' : []
}


# actors
volume_actor = vtk.vtkVolume()
axis_actor = vtk.vtkActor()
contour_actor = vtk.vtkActor()
centerline_actors = []
bathy_actor = vtk.vtkActor()