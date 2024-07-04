import sys
from PyQt6.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout
from PyQt6.QtWidgets import QPushButton, QLabel
import vtk, vtkmodules
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import config
from processing import feature_detection, read_from_mat, rendering_actors


class VolumePanel(QWidget):
    is_draw = True

    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        label = QLabel('Volume:')
        main_layout.addWidget(label)

        button_layout = QHBoxLayout()
        main_layout.addLayout(button_layout)
        draw_button = QPushButton('draw')
        draw_button.clicked.connect(self.draw_volume)
        button_layout.addWidget(draw_button)
        delete_button = QPushButton('delete')
        delete_button.clicked.connect(self.delete_volume)
        button_layout.addWidget(delete_button)

    def draw_volume(self):
        if not self.is_draw:
            config.volume_actor = rendering_actors.volume(config.imaging['vtk_data'])
            config.render.AddActor(config.volume_actor)
            config.interactor.GetRenderWindow().Render()
            self.is_draw = True
    def delete_volume(self):
        config.render.RemoveActor(config.volume_actor)
        config.interactor.GetRenderWindow().Render()
        self.is_draw = False

class AxisPanel(QWidget):
    is_draw = False

    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        label = QLabel('Axis:')
        main_layout.addWidget(label)

        button_layout = QHBoxLayout()
        main_layout.addLayout(button_layout)
        draw_button = QPushButton('draw')
        draw_button.clicked.connect(self.draw_axis)
        button_layout.addWidget(draw_button)
        delete_button = QPushButton('delete')
        delete_button.clicked.connect(self.delete_axis)
        button_layout.addWidget(delete_button)

    def draw_axis(self):
        if not self.is_draw:
            config.axis_actor = rendering_actors.axis()
            config.render.AddActor(config.axis_actor)
            config.interactor.GetRenderWindow().Render()
            self.is_draw = True
    def delete_axis(self):
        config.render.RemoveActor(config.axis_actor)
        config.interactor.GetRenderWindow().Render()
        self.is_draw = False


class ContourPanel(QWidget):
    is_draw = False

    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        label = QLabel('Contour:')
        main_layout.addWidget(label)

        button_layout = QHBoxLayout()
        main_layout.addLayout(button_layout)
        draw_button = QPushButton('draw')
        draw_button.clicked.connect(self.draw_contour)
        button_layout.addWidget(draw_button)
        delete_button = QPushButton('delete')
        delete_button.clicked.connect(self.delete_contour)
        button_layout.addWidget(delete_button)

    def draw_contour(self):
        if not self.is_draw:
            config.contour_actor = rendering_actors.contour(config.imaging['vtk_data'], [0.001])
            config.render.AddActor(config.contour_actor)
            config.interactor.GetRenderWindow().Render()
            self.is_draw = True
    def delete_contour(self):
        config.render.RemoveActor(config.contour_actor)
        config.interactor.GetRenderWindow().Render()
        self.is_draw = False
       

class CenterlinePanel(QWidget):
    is_draw = False
    pth = 100

    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        label = QLabel('Center Line:')
        main_layout.addWidget(label)

        button_layout = QHBoxLayout()
        main_layout.addLayout(button_layout)
        draw_button = QPushButton('draw')
        draw_button.clicked.connect(self.draw_centerline)
        button_layout.addWidget(draw_button)
        delete_button = QPushButton('delete')
        delete_button.clicked.connect(self.delete_centerline)
        button_layout.addWidget(delete_button)

    def draw_centerline(self):
        if not self.is_draw:
            config.centerline_actors = rendering_actors.centerline()
            for act in config.centerline_actors:
                config.render.AddActor(act)
            config.interactor.GetRenderWindow().Render()
            self.is_draw = True
    def delete_centerline(self):
        for act in config.centerline_actors:
            config.render.RemoveActor(act)
        config.interactor.GetRenderWindow().Render()
        self.is_draw = False



class RightPanel(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)

        volume_pannel = VolumePanel()
        layout.addWidget(volume_pannel)

        axis_pannel = AxisPanel()
        layout.addWidget(axis_pannel)

        contour_panel = ContourPanel()
        layout.addWidget(contour_panel)

        centerline_panel = CenterlinePanel()
        layout.addWidget(centerline_panel)

        

