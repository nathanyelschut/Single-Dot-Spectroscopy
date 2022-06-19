from PyQt5.QtWidgets import QProgressBar, QFileDialog, QGroupBox, QFormLayout, QLineEdit, QMainWindow, QApplication, QPushButton, QSpinBox, QDoubleSpinBox, QWidget, QTabWidget, QGridLayout, QVBoxLayout, QHBoxLayout, QComboBox, QLabel
import sys 
import numpy as np 
from math import ceil, floor
import pyqtgraph as pg 
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import QDot_spectroscopy_Model

# Graph settings
pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "k")

class Worker(QObject):
    """Class for assigning a longer task to a thread. 
    
    In this case the data analysis is assigned to the worker
    thread.
    """
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def run(self, CCD_width, CCD_height, CCD_width_corr, CCD_height_corr, lower_wavelength, upper_wavelength, filenames, meas_corr_curve_file_path):
        """Method which actually runs the data analysis.
        
        Arguments:
            CCD_width (int): number of pixels along the vertical axis of the CCD.
            CCD_height (int): number of pixels along the vertical axis of the CCD.
            CCD_width_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
            CCD_height_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
            lower_wavelength (float): the lowest wavelength of the spectrum in nanometers.
            upper_wavelength (float): the highest wavelength of the spectrum in nanometers.
            filenames (list of strings or paths): paths of the datafiles to be analyzed.
            meas_corr_curve_file_path (string or path): path of the correction file for correcting the measurement data.

        Returns:
            The data for making the plots and histograms.
        """

        # Run file analysis and show plots
        all_files_data, FWHM_hist_data, central_energy_hist_data = QDot_spectroscopy_Model.multiple_files_analysis(lower_wavelength, upper_wavelength, CCD_height, CCD_width, CCD_height_corr, CCD_width_corr, filenames, meas_corr_curve_file_path, self.progress)

        self.finished.emit()

        return all_files_data, FWHM_hist_data, central_energy_hist_data

class Window(QMainWindow): 
    """Class for making a GUI for a data anlysis program for single quantum dot spectroscopy.
    
    The Class creates a gui with 3 tabs to view several 2d maps of the CCD image, the spectrum of
    quantum dots and to see histograms of the FWHM and central energy of all the quantum dots in the
    PL-image.
    """
    def __init__(self): 
        """Initializes global variables and creates the frame of the application."""
        
        super().__init__() 

        self.corr_file_loaded = False
        self.files_loaded = False
        self.data_is_available = False
        self.counter = 0

        self.all_files_data = []
        self.FWHM_hist_data = []
        self.central_energy_hist_data = []

        self.titles = ["Raw data", "No background noise, non-corrected data", "Raw corrected data", "No background noise, corrected data"]

        # Setting geometry 
        self.setGeometry(100, 50, 1700, 900) 

        # Setting title 
        self.window = QMainWindow(self)
        self.setCentralWidget(self.window)
        self.setWindowTitle("Single Quantum Dot Spectroscopy") 
  
        # Calling method 
        self.UiComponents() 
  
        # Showing all the widgets 
        self.show() 
  
    def UiComponents(self): 
        """Method for creating UI components such as tabs, color maps and sidepanel items.
        
        This method specifically includes a color map for 2d maps, 3 tabbed windows (one for 2d maps, one
        for graphs and one for histograms) and a sidepanel containing options to select files and to fill in measurement
        settings.
        """
  
        # Creating a widget object 
        widget = QWidget() 
  
        # Setting configuration options 
        pg.setConfigOptions(antialias=True) 

        # --- Color map (for 2d maps) ---
        self.colors = [ 
            (68.08602, 1.24287, 84.000825),
            (64.75342500000001, 67.63977, 135.145665),
            (41.724374999999995, 120.13891500000001, 142.32774),
            (34.34646, 167.95218, 132.000495),
            (121.76352, 209.46821999999997, 81.139725),
            (253.27824, 231.070035, 36.703680000000006)
        ]
        self.cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=self.colors) 
    	# --------------------------------


        # Creating a grid layout and adding plots to grid
        layout = QHBoxLayout()
        file_box = QGroupBox("Load Files")
        settings_box = QGroupBox("Measurement Settings")

        # --- Tabbed windows ---
        # Initialize tabs
        tabwidget = QTabWidget()
        tabwidget.setMinimumWidth(1200)
        tabwidget.setMaximumWidth(1201)
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()

        # Add tabs to tab widget
        tabwidget.addTab(self.tab1, "2d maps")
        tabwidget.addTab(self.tab2, "Data fitting")
        tabwidget.addTab(self.tab3, "Histogram")

        # Create first tab
        self.tab1.layout = QGridLayout()
        self.tab2.layout = QGridLayout()
        self.tab3.layout = QVBoxLayout()
        self.tab3.top_layout = QHBoxLayout()
        self.tab3.bottom_layout = QHBoxLayout()
        self.tab3.layout.addLayout(self.tab3.top_layout)
        self.tab3.layout.addLayout(self.tab3.bottom_layout)

        # Add layouts to tabs
        self.tab1.setLayout(self.tab1.layout)
        self.tab2.setLayout(self.tab2.layout)
        self.tab3.setLayout(self.tab3.layout)
        # ----------------------
        

        # --- User interface Sidepanel ---
        # Create sidepanel items
        self.files_combobox = QComboBox() 
        self.files_combobox.currentIndexChanged.connect(self.update_plots)

        # Buttons / status updates
        self.button_groupbox = QGroupBox()
        self.file_label = QLabel("Select file for analysis")
        self.open_files_button = QPushButton("Load Files", default=False, autoDefault=False)
        self.open_files_button.clicked.connect(self.load_files)
        self.open_corr_file_button = QPushButton("Load Correction File", default=False, autoDefault=False)
        self.open_corr_file_button.clicked.connect(self.load_corr_file)
        self.corr_file_status = QLineEdit()
        self.corr_file_status.setReadOnly(True)
        self.corr_file_status.setText("No correction file is loaded")
        self.files_status = QLineEdit()
        self.files_status.setReadOnly(True)
        self.files_status.setText("No files loaded")
        self.run_analysis_button = QPushButton("Run Analysis")
        self.run_analysis_button.clicked.connect(self.start_analysis)
        self.pbar = QProgressBar()

        # Measurement settings user interface
        self.lower_wavelength_label = QLabel("Lower wavelength (in nm)")
        self.lower_wavelength_dspinbox = QDoubleSpinBox()
        self.lower_wavelength_dspinbox.setMaximum(1000000)
        self.lower_wavelength_dspinbox.setDecimals(8)
        self.lower_wavelength_dspinbox.setValue(610.00929128)

        self.upper_wavelength_label = QLabel("Upper wavelength (in nm)")
        self.upper_wavelength_dspinbox = QDoubleSpinBox()
        self.upper_wavelength_dspinbox.setMaximum(1000000)
        self.upper_wavelength_dspinbox.setDecimals(8)
        self.upper_wavelength_dspinbox.setValue(886.93139111)

        self.CCD_resolution_label = QLabel("CCD resolution of data files (x, y)")
        self.coordinate_label1 = QLabel("(")
        self.coordinate_label1.setMaximumWidth(5)
        self.coordinate_label2 = QLabel(",")
        self.coordinate_label2.setMaximumWidth(5)
        self.coordinate_label3 = QLabel(")")
        self.coordinate_label3.setMaximumWidth(5)
        self.CCD_resolution_x_spinbox = QSpinBox()
        self.CCD_resolution_x_spinbox.setMaximumWidth(60)
        self.CCD_resolution_x_spinbox.setMaximum(100000)
        self.CCD_resolution_x_spinbox.setValue(167)
        self.CCD_resolution_y_spinbox = QSpinBox()
        self.CCD_resolution_y_spinbox.setMaximumWidth(60)
        self.CCD_resolution_y_spinbox.setMaximum(100000)
        self.CCD_resolution_y_spinbox.setValue(200)

        self.CCD_resolution_corr_label = QLabel("CCD resolution of correction files (x, y)")
        self.coordinate_label4 = QLabel("(")
        self.coordinate_label4.setMaximumWidth(5)
        self.coordinate_label5 = QLabel(",")
        self.coordinate_label5.setMaximumWidth(5)
        self.coordinate_label6 = QLabel(")")
        self.coordinate_label6.setMaximumWidth(5)
        self.CCD_resolution_corr_x_spinbox = QSpinBox()
        self.CCD_resolution_corr_x_spinbox.setMaximumWidth(60)
        self.CCD_resolution_corr_x_spinbox.setMaximum(100000)
        self.CCD_resolution_corr_x_spinbox.setValue(1340)
        self.CCD_resolution_corr_y_spinbox = QSpinBox()
        self.CCD_resolution_corr_y_spinbox.setMaximumWidth(60)
        self.CCD_resolution_corr_y_spinbox.setMaximum(100000)
        self.CCD_resolution_corr_y_spinbox.setValue(400)

        # Create sidepanel
        self.sidepanel = QWidget()
        self.sidepanel.layout = QVBoxLayout()                # general layout
        self.sidepanel.file_layout = QVBoxLayout()           # file selection layout (for groupbox)
        self.sidepanel.select_file_layout = QFormLayout()    # sublayout for drop down menu
        self.sidepanel.settings_layout = QVBoxLayout()       # settings layout
        self.sidepanel.settings_layout1 = QFormLayout()      # sublayout for wavelength settings
        self.sidepanel.settings_layout2 = QGridLayout()      # sublayout for CCD settings
        self.sidepanel.run_analysis_layout = QVBoxLayout()   # sublayout for analysis button and pbar

        self.sidepanel.file_layout.addWidget(self.open_corr_file_button)
        self.sidepanel.file_layout.addWidget(self.corr_file_status)
        self.sidepanel.file_layout.addWidget(self.open_files_button)
        self.sidepanel.file_layout.addWidget(self.files_status)
        self.sidepanel.file_layout.addLayout(self.sidepanel.select_file_layout)
        file_box.setLayout(self.sidepanel.file_layout)
        self.sidepanel.layout.addWidget(file_box)
        self.sidepanel.layout.addSpacing(20)
        self.sidepanel.layout.addWidget(settings_box)
        self.sidepanel.layout.addSpacing(20)
        self.sidepanel.run_analysis_layout.addWidget(self.run_analysis_button)
        self.sidepanel.layout.addLayout(self.sidepanel.run_analysis_layout)

        self.sidepanel.select_file_layout.addRow(self.file_label, self.files_combobox)
        self.sidepanel.settings_layout1.addRow(self.lower_wavelength_label, self.lower_wavelength_dspinbox)
        self.sidepanel.settings_layout1.addRow(self.upper_wavelength_label, self.upper_wavelength_dspinbox)
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_label, 0, 0)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label1, 0, 1)
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_x_spinbox, 0, 2)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label2, 0, 3)
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_y_spinbox, 0, 4)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label3, 0, 5)
        
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_corr_label, 1, 0)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label4, 1, 1)
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_corr_x_spinbox, 1, 2)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label5, 1, 3)
        self.sidepanel.settings_layout2.addWidget(self.CCD_resolution_corr_y_spinbox, 1, 4)
        self.sidepanel.settings_layout2.addWidget(self.coordinate_label6, 1, 5)

        self.sidepanel.settings_layout.addLayout(self.sidepanel.settings_layout1)
        self.sidepanel.settings_layout.addLayout(self.sidepanel.settings_layout2)
        settings_box.setLayout(self.sidepanel.settings_layout)

        self.sidepanel.layout.insertSpacing(0, 5)
        self.sidepanel.layout.addStretch()
        self.sidepanel.setLayout(self.sidepanel.layout)
        # ---------------------------------
        

        # Add tab widget and sidepanel to central widget
        layout.addWidget(tabwidget)
        layout.addWidget(self.sidepanel)
        
        # Stting this layout to the widget 
        widget.setLayout(layout)
        
        # Setting this widget as central widget of the main widow 
        self.setCentralWidget(widget)
    
    def load_corr_file(self):
        """Method for loading in the correction file."""

        # Dialog window
        self.meas_corr_curve_file_path, _ = QFileDialog.getOpenFileName(filter="CSV files(*.csv)")
        
        if self.meas_corr_curve_file_path != '':
            # Status box update
            shortfilename = get_file_name_from_path(self.meas_corr_curve_file_path)
            self.corr_file_status.setText(f"Succesfully loaded {shortfilename}")
            self.files_status.setText("No files loaded")
            self.corr_file_loaded = True

        else:
            self.corr_file_status.setText("No correction file is loaded")
            self.corr_file_loaded = False

    def load_files(self):
        """Method for loading in data files."""

        # Dialog window
        self.filenames, _ = QFileDialog.getOpenFileNames(filter="CSV files(*.csv)")
        
        # Dropdown menu
        self.shortfilenames = []
        for filename in self.filenames:
            shortfilename = get_file_name_from_path(filename)
            self.shortfilenames.append(shortfilename)
        
        self.files_combobox.clear()
        self.files_combobox.addItems(self.shortfilenames)

        # Update status
        if len(self.filenames) > 0:
            self.files_status.setText(f"Succesfully loaded files")
            self.files_loaded = True

        else:
            self.files_status.setText(f"No files loaded")
            self.files_loaded = False
            
        self.clear_plots(self.tab1.layout)
        self.clear_plots(self.tab2.layout)
        
        if self.tab3.top_layout.count() > 1:
            self.clear_plots(self.tab3.sidepanel_top)
            self.tab3.top_layout.removeItem(self.tab3.sidepanel_top)
            self.clear_plots(self.tab3.top_layout)
            self.clear_plots(self.tab3.sidepanel_bottom)
            self.tab3.bottom_layout.removeItem(self.tab3.sidepanel_bottom)
            self.clear_plots(self.tab3.bottom_layout)

    def get_spinbox_values(self):
        """Method to retreive the current values of the measurement spinboxes.
        
        Returns:
            All the values in the spinboxes, contained in a tuple.
        """
        self.CCD_width = self.CCD_resolution_x_spinbox.value()
        self.CCD_height = self.CCD_resolution_y_spinbox.value()

        self.CCD_width_corr = self.CCD_resolution_corr_x_spinbox.value()
        self.CCD_height_corr = self.CCD_resolution_corr_y_spinbox.value()

        self.lower_wavelength = self.lower_wavelength_dspinbox.value()
        self.upper_wavelength = self.upper_wavelength_dspinbox.value()
        
        return self.CCD_width, self.CCD_height, self.CCD_width_corr, self.CCD_height_corr, self.lower_wavelength, self.upper_wavelength

    def start_analysis(self):
        """Method for running the analysis on the files.
        
        This method starts a worker thread for the data analysis on the data files, and adds
        a progress bar to the GUI when the data is being analyzed.
        """

        if self.corr_file_loaded == True and self.files_loaded == True:
            self.CCD_width, self.CCD_height, self.CCD_width_corr, self.CCD_height_corr, self.lower_wavelength, self.upper_wavelength = self.get_spinbox_values()
            self.sidepanel.run_analysis_layout.addWidget(self.pbar)

            self.thread = QThread()
            self.worker = Worker()
            self.worker.moveToThread(self.thread)

            def start_long_func():
                self.all_files_data, self.FWHM_hist_data, self.central_energy_hist_data, = self.worker.run(self.CCD_width, self.CCD_height, self.CCD_width_corr, self.CCD_height_corr, self.lower_wavelength, self.upper_wavelength, self.filenames, self.meas_corr_curve_file_path)
            
            self.thread.started.connect(start_long_func)
            self.worker.finished.connect(self.thread.quit)
            self.worker.finished.connect(self.worker.deleteLater)
            self.thread.finished.connect(self.thread.deleteLater)
            self.worker.progress.connect(self.report_progress)

            self.thread.start()
            self.run_analysis_button.setEnabled(False)
            self.thread.finished.connect(
                lambda: self.run_analysis_button.setEnabled(True)
            )
            
            def set_boolean():
                self.data_is_available = True

            def remove_pbar():
                self.pbar.setParent(None)

            self.thread.finished.connect(remove_pbar)
            self.thread.finished.connect(set_boolean)
            self.thread.finished.connect(self.update_plots)

    def update_plots(self):
        """Method for updating plot widgets.
        
        This method creates the plot items for the 2d maps and the QDot graphs.
        """

        if self.data_is_available is True:
            file_index = self.files_combobox.currentIndex()
            file_data = self.all_files_data[file_index]
            
            self.mapdata = [file_data[i].values for i in range(4)]
            self.Q_Dot_plot_data = file_data[4]

            for i in reversed(range(self.tab1.layout.count())): 
                self.tab1.layout.itemAt(i).widget().setParent(None)

            # Create ticks for x axis
            tick_strings_lower = floor(file_data[0].reset_index()['index'].loc[0] * 10) / 10
            tick_strings_upper = ceil(file_data[0].reset_index()['index'].loc[self.CCD_width - 1] * 10) / 10
            tick_strings = np.arange(tick_strings_lower, tick_strings_upper, -0.1)
            
            energies = file_data[0].reset_index()['index'].to_numpy()
            
            tick_values = []
            for E in tick_strings:
                tick_value = (E - energies[0]) * (self.CCD_width / (energies[-1] - energies[0]))
                tick_values.append(tick_value)

            tick_tuples = []
            for i in range(len(tick_values)):
                tick_tuples.append((float(tick_values[i]), "{:.1f}".format(tick_strings[i])))

            for i in reversed(range(self.tab2.layout.count())): 
                self.tab2.layout.itemAt(i).widget().setParent(None)

            # Create 2d maps image views
            rows = [1, 1, 4, 4]
            columns = [1, 2, 1, 2]
            plots_2dmap = {}
            axis_2dmap = {}
            imv_2dmap = {}
            for i in range(4):
                axis_2dmap[i] = pg.AxisItem('bottom')
                axis_2dmap[i].setTicks([tick_tuples])

                plots_2dmap[i] = pg.PlotItem(title=self.titles[i], axisItems={"bottom": axis_2dmap[i]})
                plots_2dmap[i].setLabel(axis='left', text='Slit position')
                plots_2dmap[i].setLabel(axis='bottom', text='Energy (eV)')
                imv_2dmap[i] = pg.ImageView(view=plots_2dmap[i]) 
                imv_2dmap[i].setImage(self.mapdata[i]) 
                imv_2dmap[i].setColorMap(self.cmap) 
                # imv_2dmap[i].ui.histogram.hide()
                imv_2dmap[i].ui.roiBtn.hide()
                imv_2dmap[i].ui.menuBtn.hide()

                # Add plots to layout of tab1
                self.tab1.layout.addWidget(imv_2dmap[i], rows[i], columns[i], 3, 1)

            # Create plot items for each QDot
            self.Q_Dot_plot_data_arr = self.Q_Dot_plot_data.swapaxes('index', 'columns').to_numpy()
            column_headers = list(self.Q_Dot_plot_data.columns)
            counter = 1
            plots_Qdot = {}
            axis_Qdot = {}
            for QdotIndex in range(0, len(self.Q_Dot_plot_data_arr), 2):
                row = floor((counter - 1) / 3)
                column = counter - ((row * 3)) - 1

                axis_Qdot[QdotIndex / 2] = pg.AxisItem('bottom')
                axis_Qdot[QdotIndex / 2].setTicks([tick_tuples])

                plots_Qdot[QdotIndex / 2] = pg.PlotWidget(title=f"Slit position: {column_headers[QdotIndex]}", axisItems={"bottom": axis_Qdot[QdotIndex / 2]})
                plots_Qdot[QdotIndex / 2].addLegend()
                plots_Qdot[QdotIndex / 2].plot(self.Q_Dot_plot_data_arr[QdotIndex], symbol="o", symbolSize=6, pen=None, name="PL-spectrum")
                plots_Qdot[QdotIndex / 2].plot(self.Q_Dot_plot_data_arr[QdotIndex + 1], pen=pg.mkPen('r'), name="Optimal gaussian fit")
                plots_Qdot[QdotIndex / 2].setLabel(axis='left', text='Intensity')
                plots_Qdot[QdotIndex / 2].setLabel(axis='bottom', text='Energy (eV)')
                

                # Add plots to layout of tab2
                self.tab2.layout.addWidget(plots_Qdot[QdotIndex / 2], row, column, 1, 1)
                counter += 1

            if self.tab3.top_layout.count() == 0:
                self.make_histogram()
    
    def clear_plots(self, layout):
        for i in reversed(range(layout.count())): 
            layout.itemAt(i).widget().setParent(None)

    def make_histogram(self):
        """Method for creating the FWHM and Central energy histograms."""
       
        # Create sidepanel for top histogram
        self.tab3.sidepanel_top = QFormLayout(self)
        
        self.FWHM_bin_spinbox = QSpinBox()
        self.FWHM_bin_spinbox.setRange(1, 100)
        self.FWHM_bin_spinbox.setValue(30)
        self.FWHM_bin_spinbox.valueChanged.connect(self.update_histograms)
        
        self.FWHM_min_spinbox = QDoubleSpinBox()
        self.FWHM_min_spinbox.setDecimals(2)
        self.FWHM_min_spinbox.setRange(0.00, 999.99)
        self.FWHM_min_spinbox.setSingleStep(0.1)
        self.FWHM_min_spinbox.setValue(0)
        self.FWHM_min_spinbox.valueChanged.connect(self.update_histograms)

        self.FWHM_max_spinbox = QDoubleSpinBox()
        self.FWHM_max_spinbox.setDecimals(2)
        self.FWHM_max_spinbox.setRange(0.01, 1000.00)
        self.FWHM_max_spinbox.setSingleStep(0.1)
        self.FWHM_max_spinbox.setValue(1000.00)
        self.FWHM_max_spinbox.valueChanged.connect(self.update_histograms)
        
        self.tab3.sidepanel_top.addRow(QLabel("Bins:"), self.FWHM_bin_spinbox)
        self.tab3.sidepanel_top.addRow(QLabel("Minimum FWHM:"), self.FWHM_min_spinbox)
        self.tab3.sidepanel_top.addRow(QLabel("Maximum FWHM:"), self.FWHM_max_spinbox)

        # Create histogram for FWHM
        self.FWHM_hist = pg.PlotWidget(title="FWHM histogram of all QDots")
        FWHM_bins = self.FWHM_bin_spinbox.value()
        FWHM_min = self.FWHM_min_spinbox.value()
        FWHM_max = self.FWHM_max_spinbox.value()
        x, y, x_new, best_fit = QDot_spectroscopy_Model.histogram(self.FWHM_hist_data, FWHM_bins, FWHM_min, FWHM_max)
        
        FWHM_curve = pg.PlotCurveItem(x, y, stepMode=True, fillLevel=0, brush=(0, 0, 255, 80))
        FWHM_fit = pg.PlotDataItem(x_new, best_fit, pen=pg.mkPen('r'))
        self.FWHM_hist.addItem(FWHM_curve)
        self.FWHM_hist.addItem(FWHM_fit)
        self.FWHM_hist.setLabel(axis='left', text='Frequency')
        self.FWHM_hist.setLabel(axis='bottom', text='FWHM (eV)')
        self.tab3.top_layout.addWidget(self.FWHM_hist)
        self.tab3.top_layout.addLayout(self.tab3.sidepanel_top)

        # Create sidepanel for bottom histogram
        self.tab3.sidepanel_bottom = QFormLayout(self)
        
        self.c_energy_bin_spinbox = QSpinBox()
        self.c_energy_bin_spinbox.setRange(1, 100)
        self.c_energy_bin_spinbox.setValue(30)
        self.c_energy_bin_spinbox.valueChanged.connect(self.update_histograms)
        
        self.c_energy_min_spinbox = QDoubleSpinBox()
        self.c_energy_min_spinbox.setDecimals(2)
        self.c_energy_min_spinbox.setRange(-1000.00, 999.99)
        self.c_energy_min_spinbox.setSingleStep(0.1)
        self.c_energy_min_spinbox.setValue(-100)
        self.c_energy_min_spinbox.valueChanged.connect(self.update_histograms)

        self.c_energy_max_spinbox = QDoubleSpinBox()
        self.c_energy_max_spinbox.setDecimals(2)
        self.c_energy_max_spinbox.setRange(0.01, 1000.00)
        self.c_energy_max_spinbox.setSingleStep(0.1)
        self.c_energy_max_spinbox.setValue(1000.00)
        self.c_energy_max_spinbox.valueChanged.connect(self.update_histograms)
        
        self.tab3.sidepanel_bottom.addRow(QLabel("Bins:"), self.c_energy_bin_spinbox)
        self.tab3.sidepanel_bottom.addRow(QLabel("Minimum central energy:"), self.c_energy_min_spinbox)
        self.tab3.sidepanel_bottom.addRow(QLabel("Maximum central energy:"), self.c_energy_max_spinbox)

        # Create histogram for central energies
        self.c_energy_hist = pg.PlotWidget(title="Central energy histogram of all QDots")
        c_energy_bins = self.c_energy_bin_spinbox.value()
        c_energy_min = self.c_energy_min_spinbox.value()
        c_energy_max = self.c_energy_max_spinbox.value()
        x, y, x_new, best_fit = QDot_spectroscopy_Model.histogram(self.central_energy_hist_data, c_energy_bins, c_energy_min, c_energy_max)
        
        c_energy_curve = pg.PlotCurveItem(x, y, stepMode=True, fillLevel=0, brush=(0, 0, 255, 80))
        c_energy_fit = pg.PlotDataItem(x_new, best_fit, pen=pg.mkPen('r'))
        self.c_energy_hist.addItem(c_energy_curve)
        self.c_energy_hist.addItem(c_energy_fit)
        self.c_energy_hist.setLabel(axis='left', text='Frequency')
        self.c_energy_hist.setLabel(axis='bottom', text='Central energy (eV)')

        self.tab3.bottom_layout.addWidget(self.c_energy_hist)
        self.tab3.bottom_layout.addLayout(self.tab3.sidepanel_bottom)

    def update_histograms(self):
        """Method for updating histograms based on sidepanel settings."""

        self.FWHM_hist.clear()

        FWHM_bins = self.FWHM_bin_spinbox.value()
        FWHM_min = self.FWHM_min_spinbox.value()
        FWHM_max = self.FWHM_max_spinbox.value()

        x, y, x_new, best_fit = QDot_spectroscopy_Model.histogram(self.FWHM_hist_data, FWHM_bins, FWHM_min, FWHM_max)

        FWHM_curve = pg.PlotCurveItem(x, y, stepMode=True, fillLevel=0, brush=(0, 0, 255, 80))
        FWHM_fit = pg.PlotDataItem(x_new, best_fit, pen=pg.mkPen('r'))
        self.FWHM_hist.addItem(FWHM_curve)
        self.FWHM_hist.addItem(FWHM_fit)

        self.c_energy_hist.clear()

        c_energy_bins = self.c_energy_bin_spinbox.value()
        c_energy_min = self.c_energy_min_spinbox.value()
        c_energy_max = self.c_energy_max_spinbox.value()

        x, y, x_new, best_fit = QDot_spectroscopy_Model.histogram(self.central_energy_hist_data, c_energy_bins, c_energy_min, c_energy_max)

        c_energy_curve = pg.PlotCurveItem(x, y, stepMode=True, fillLevel=0, brush=(0, 0, 255, 80))
        c_energy_fit = pg.PlotDataItem(x_new, best_fit, pen=pg.mkPen('r'))
        self.c_energy_hist.addItem(c_energy_curve)
        self.c_energy_hist.addItem(c_energy_fit)

    def report_progress(self, n):
        """Method for updating the progress bar in the GUI."""
        self.pbar.setValue(n)

def get_file_name_from_path(path):
    """Returns file name from the full path
    
    Arguments:
        path (str): full path name or similar to parent_directory/directory/file_name.

    Returns:
        A string containing only the file name
    """
    pos_slash = [pos for pos, char in enumerate(path) if char == '/']
    file_name = path[pos_slash[-1]+1:]

    return file_name

# create pyqt5 app 
App = QApplication(sys.argv) 
  
# create the instance of our Window 
window = Window() 
  
# start the app 
sys.exit(App.exec()) 