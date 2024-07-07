import vtk, vtkmodules
import numpy as np
import joblib
import cv2
from vtkmodules.util import numpy_support, vtkImageImportFromArray
import scipy
from processing import read_from_mat, feature_detection
import time
import config

def axis(L=100):
    axis = vtk.vtkAxesActor()
    axis.SetTotalLength(L, L, L)
    axis.AxisLabelsOff()
    axis.SetAxisLabels(1)
    axis.SetShaftTypeToCylinder()
    axis.SetCylinderRadius(0.001)
    axis.GetXAxisShaftProperty().SetColor(0, 1, 0)
    axis.GetYAxisShaftProperty().SetColor(0, 1, 0)
    axis.GetZAxisShaftProperty().SetColor(0, 1, 0)
    axis.GetXAxisTipProperty().SetColor(0, 1, 0)
    axis.GetYAxisTipProperty().SetColor(0, 1, 0)
    axis.GetZAxisTipProperty().SetColor(0, 1, 0)
    return axis


def volume(vtk_data):
    color_func = vtk.vtkColorTransferFunction()
    color_func.AddRGBPoint(0.0, 0.0, 0.0, 0.0)    
    color_func.AddRGBPoint(0.00001, 1.0, 1.0, 1.0)

    opacity_func = vtk.vtkPiecewiseFunction()
    opacity_func.AddPoint(0, 0)
    opacity_func.AddPoint(0.00005, 0.8)

    volume_mapper = vtk.vtkGPUVolumeRayCastMapper()   
    volume_mapper.SetInputData(vtk_data)

    volume_property = vtk.vtkVolumeProperty()
    volume_property.SetColor(color_func)
    volume_property.SetScalarOpacity(opacity_func)
    volume_property.SetInterpolationTypeToLinear()
    volume_property.ShadeOn()

    vol = vtk.vtkVolume()
    vol.SetMapper(volume_mapper)
    vol.SetProperty(volume_property)   
    return vol 



def contour(vtk_data, values, color=(0.0,1.0,0.0)):
    contour = vtk.vtkContourFilter()
    contour.SetInputData(vtk_data)
    for i in range(len(values)):
        contour.SetValue(i, values[i])

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(contour.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    return actor



def lines(points, width=2, color=(1.0,1.0,1.0)):
    line_source = vtk.vtkLineSource()
    vtk_points = vtk.vtkPoints()
    for p in points:
        vtk_points.InsertNextPoint(p)
    line_source.SetPoints(vtk_points)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line_source.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetLineWidth(width)
    actor.GetProperty().SetColor(color)
    return actor


def centerline(th=100, color=(1.0,0.0,0.0)):
    spacing = config.imaging['spacing']
    bounds = config.imaging['cut_bounds']
    regions = feature_detection.regions
    actors = []
    for r in regions:
        if r['count'] > th:
            points = feature_detection.get_center_points(r['id'])
            points = feature_detection.curve_fitting(points)
            for i in range(len(points)):
                points[i] = (points[i][0]*spacing[0]+bounds[0],points[i][1]*spacing[1]+bounds[2],points[i][2]*spacing[2]+bounds[4])
            act = lines(points, color=color)
            actors.append(act)
    return actors


def bathy(color=(1.0,1.0,0.0)):
    bathy = config.bathy['data']
    spacing = config.bathy['spacing']
    bounds = config.bathy['bounds']
    points = vtk.vtkPoints()
    points.Allocate(bathy.ravel().shape[0])
    for i in range(bathy.shape[0]):
        for j in range(bathy.shape[1]):
            points.InsertNextPoint(i*spacing[0]+bounds[0], j*spacing[0]+bounds[2], bathy[i][j])

    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points)
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(poly_data)
    delaunay.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(delaunay.GetOutputPort())
    config.bathy_actor.SetMapper(mapper)
    config.bathy_actor.GetProperty().SetColor(color)

    



