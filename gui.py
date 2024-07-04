import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QVBoxLayout, QWidget, QLabel, QComboBox, QLineEdit, QPushButton
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
from matplotlib.widgets import RectangleSelector
from trajectory import TrajectoryProcessor
from analysis import CDFAnalysis
from config import Config
from plotting import Plotting

class CDFGui(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.config = Config()
        self.traj_processor = None
        self.analysis = None
        self.plotting = None
        self.cdf_data = None

    def initUI(self):
        self.setWindowTitle('CDF Visualizer')
        self.setGeometry(100, 100, 1200, 800)

        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        self.layout = QVBoxLayout(self.main_widget)

        self.file_label = QLabel('Select Topology and Trajectory Files:')
        self.layout.addWidget(self.file_label)

        self.topology_button = QPushButton('Select Topology File')
        self.topology_button.clicked.connect(self.select_topology_file)
        self.layout.addWidget(self.topology_button)

        self.trajectory_button = QPushButton('Select Trajectory File')
        self.trajectory_button.clicked.connect(self.select_trajectory_file)
        self.layout.addWidget(self.trajectory_button)

        self.parameter_label = QLabel('Set Parameters:')
        self.layout.addWidget(self.parameter_label)

        self.mode_dropdown = QComboBox(self)
        self.mode_dropdown.addItems(['default', 'hybrid'])
        self.layout.addWidget(self.mode_dropdown)

        self.sel_dropdown = QComboBox(self)
        self.sel_dropdown.addItems(['com', 'name', 'element'])
        self.layout.addWidget(self.sel_dropdown)

        self.selname_entry = QLineEdit(self)
        self.selname_entry.setPlaceholderText('Selection Name')
        self.layout.addWidget(self.selname_entry)

        self.selelement_entry = QLineEdit(self)
        self.selelement_entry.setPlaceholderText('Selection Element')
        self.layout.addWidget(self.selelement_entry)

        self.length_entry = QLineEdit(self)
        self.length_entry.setPlaceholderText('Cutoff Length (default: 40)')
        self.layout.addWidget(self.length_entry)

        self.radius_entry = QLineEdit(self)
        self.radius_entry.setPlaceholderText('Cutoff Radius (default: 15)')
        self.layout.addWidget(self.radius_entry)

        self.res_entry = QLineEdit(self)
        self.res_entry.setPlaceholderText('Resolution (default: 4)')
        self.layout.addWidget(self.res_entry)

        self.plot_button = QPushButton('Plot CDF')
        self.plot_button.clicked.connect(self.plot_cdf)
        self.layout.addWidget(self.plot_button)

        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

    def select_topology_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Select Topology File", "", "All Files (*);;Topology Files (*.topo)", options=options)
        if file_name:
            self.config.args.topology = file_name
            self.file_label.setText(f'Selected Topology File: {file_name}')

    def select_trajectory_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Select Trajectory File", "", "All Files (*);;Trajectory Files (*.traj)", options=options)
        if file_name:
            self.config.args.trajectory = file_name
            self.file_label.setText(f'Selected Trajectory File: {file_name}')

    def plot_cdf(self):
        self.config.args.selection_mode = self.mode_dropdown.currentText()
        self.config.args.selection = self.sel_dropdown.currentText()
        self.config.args.selection_name = self.selname_entry.text()
        self.config.args.selection_element = self.selelement_entry.text()
        self.config.args.cutoff_length = int(self.length_entry.text() or 40)
        self.config.args.cutoff_radius = int(self.radius_entry.text() or 15)
        self.config.args.res = int(self.res_entry.text() or 4)

        self.traj_processor = TrajectoryProcessor(self.config)
        self.analysis = CDFAnalysis(self.config, self.traj_processor)
        self.analysis.calculate_distances()
        self.analysis.cylindrical_pcf()
        
        self.plotting = Plotting(self.config, self.analysis)
        self.cdf_data = self.plotting.generate_cdf()
    
        self.ax.clear()
        self.ax.imshow(self.cdf_data, cmap='viridis')
        self.canvas.draw()

        self.selector = RectangleSelector(self.ax, self.onselect, drawtype='box', useblit=True, button=[1], minspanx=5, minspany=5, spancoords='pixels', interactive=True)

    def onselect(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, y2 = erelease.ydata
        print(f"Selected region: ({x1}, {y1}) to ({x2}, {y2})")
        # Implement logic to find molecular pairs within this region
        self.find_molecular_pairs(x1, y1, x2, y2)

    def find_molecular_pairs(self, x1, y1, x2, y2):
        # Implement the logic from lookup.py to find and display molecular pairs
        # For demonstration purposes, we'll just print the range
        print(f"Finding pairs in range: L:({x1}, {x2}), R:({y1}, {y2})")
        # Replace with actual logic from lookup.py

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = CDFGui()
    ex.show()
    sys.exit(app.exec_())
