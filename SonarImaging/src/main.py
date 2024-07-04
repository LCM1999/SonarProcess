import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout
from PyQt6.QtWidgets import QPushButton
import vtk, vtkmodules
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import config
from processing import feature_detection, read_from_mat, rendering_actors
from gui.left_panel import LeftPanel
from gui.right_panel import RightPanel
from gui.vtk_widget import VTKWidget



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(550, 300, 1000, 600)
        self.setWindowTitle("main")
        central_widget = QWidget()
        central_widget.setStyleSheet('background-color:rgb(230,230,230);')
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout()
        central_widget.setLayout(layout)

        menu_bar = self.menuBar()
        file = menu_bar.addMenu('File')
        file.addAction('Load Imaging')
        file.addAction('Load Diffuse')

        #left_panel = LeftPanel()
        vtk_widget = VTKWidget()
        right_panel = RightPanel()
        #layout.addWidget(left_panel)
        layout.addWidget(vtk_widget)
        layout.addWidget(right_panel)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())
