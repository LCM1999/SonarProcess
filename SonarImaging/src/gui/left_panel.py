import sys
from PyQt6.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout
from PyQt6.QtWidgets import QPushButton, QSlider
from PyQt6.QtCore import Qt
import vtk, vtkmodules
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import config
from processing import feature_detection, read_from_mat, rendering_actors

def change_value(value):
    print(value)




class SliderWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)

        sldxmin = QSlider(Qt.Orientation.Horizontal, self)
        sldxmin.setRange(0,500)
        sldxmin.valueChanged.connect(change_value)
        layout.addWidget(sldxmin)

        sldxmax = QSlider(Qt.Orientation.Horizontal, self)
        sldxmax.setRange(0,500)
        sldxmax.valueChanged.connect(change_value)
        layout.addWidget(sldxmax)

        sldymin = QSlider(Qt.Orientation.Horizontal, self)
        sldymin.setRange(0,500)
        sldymin.valueChanged.connect(change_value)
        layout.addWidget(sldymin)

        sldymax = QSlider(Qt.Orientation.Horizontal, self)
        sldymax.setRange(0,500)
        sldymax.valueChanged.connect(change_value)
        layout.addWidget(sldymax)




class LeftPanel(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        slider_widget = SliderWidget()
        layout.addWidget(slider_widget)
        self.setGeometry(0,0,50,50)

       