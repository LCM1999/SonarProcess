import sys
from PyQt6.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout
from PyQt6.QtWidgets import QPushButton
import vtk, vtkmodules
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import config
from processing import data_process, feature_detection, read_from_mat, rendering_actors





class VTKWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)

        config.interactor = QVTKRenderWindowInteractor()
        config.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        layout.addWidget(config.interactor)
        config.interactor.GetRenderWindow().AddRenderer(config.render)
        #------------------------------------------------------------------------------------------------------------
 
        # imaging
        read_from_mat.get_imaging('./data/mat/imaging/COVIS-20220508T120002-imaging1.mat', 'Id_filt')
        config.imaging['data'] = config.imaging['data']
        data_process.imaging_cutting([-20,0,-20,0,0,20])
        config.imaging['vtk_data'] = data_process.to_vtk_data(config.imaging['data_cut'])
        

        config.volume_actor = rendering_actors.volume(config.imaging['vtk_data'])
        config.render.AddActor(config.volume_actor)


        # feature region
        feature_detection.caculate_region()


        # bathy(seabed)
        read_from_mat.get_bathy_from_diffuse('./data/mat/diffuse/COVIS-20220508T120002-diffuse_fr1.mat')

        rendering_actors.bathy()
        config.render.AddActor(config.bathy_actor)
 
        #------------------------------------------------------------------------------------------------------------
        config.render.ResetCamera()
        config.interactor.Initialize()
