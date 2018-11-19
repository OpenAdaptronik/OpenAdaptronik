""" Try to embed the RD estimation and plotting int PyQT5 """
import sys

import numpy as np
import pyqtgraph as pg
import scipy as sp
from PyQt5.QtCore import QRect, Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QAction, QApplication, QCheckBox, QFileDialog,
                             QGridLayout, QHBoxLayout, QLabel, QLineEdit,
                             QMainWindow, QPushButton, QRadioButton, QSlider,
                             QSpinBox, QVBoxLayout, QWidget, qApp)

import decimate_filter
import rd_estim

""" Set some global parameters"""
frame_size = 1024*30
pre_time = 0.2

cutoff_highpass = 5.0*2*np.pi

RD_order = 256*2
trigger_level = 0.01
mirror = 1

global counter
counter = 1.0

raw_data = []
RD_sequence = np.zeros(RD_order)

global win
global p0
global c0

""" Open PyQtGraphWindow"""
win=pg.GraphicsWindow(title="Data Window")
pg.setConfigOptions(antialias=True)
""" Plot: Acquired data with zoom window"""
p0 = win.addPlot(title="Time Data")
c0 = p0.plot([0,1],[0,0])
limits = pg.LinearRegionItem()
p0.addItem(limits)
win.nextRow()
""" Plot: Zoomed data for analysis"""
p1 = win.addPlot(title="Analysed Time Frame")
win.nextRow()
""" Plot: RD signature """
p2 = win.addPlot(title="Random Decrement Estimation")
win.nextRow()
""" Autopower spectrum"""
p3 = win.addPlot(title="Autopower Spectrum")



class AppForm(QMainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('RD Autocorrelation Estimator')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.RDTextbox.setText('512')
        self.spinBox.setValue(1)

        self.trigger_condition = 3

        self.time_vec = []
        self.x_vec = []
        self.y_vec = []
        self.z_vec = []
        self.activeaxis = 'z'
        self.time_zoom = [0, 0]

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"

        path = (QFileDialog.getSaveFileName(self,
                        'Save file', '',
                        file_choices))
        if path:
            exporter = pg.exporters.ImageExporter(p1)
            exporter.export(path)

    def import_csv_OA(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.windowtitle = filename

        csv_header = np.genfromtxt(filename, delimiter = ',',max_rows = 10)

        amplitude_factor = csv_header[8,1]
        sampling_factor = csv_header[6,1]

        csv_data = np.genfromtxt(filename, delimiter = ',',skip_header = 10)

        self.time_vec = csv_data[1:,0]*1E-6
        self.x_vec = csv_data[1:,1]*amplitude_factor/2.0**15
        self.y_vec = csv_data[1:,2]*amplitude_factor/2.0**15
        self.z_vec = csv_data[1:,3]*amplitude_factor/2.0**15

        """ Set the sampling rate"""
        self.rate_in = sampling_factor

        self.RecLength_textbox.setText('{:.2f}'.format(self.time_vec[-1]))
        self.Sampling_textbox.setText('{:.2f}'.format(self.rate_in))

        c0.clear()

    def import_csv_PP(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.windowtitle = filename

        csv_data = np.genfromtxt(filename, delimiter = '\t')
        self.time_vec = csv_data[1:,0]
        self.x_vec = csv_data[1:,1]
        self.y_vec = csv_data[1:,2]
        self.z_vec = csv_data[1:,3]

        """ Estimate the sampling rate"""
        t_jitter = -self.time_vec[0:len(self.time_vec)-1]+self.time_vec[1:len(self.time_vec)]
        dT = np.mean(t_jitter)
        max_jitter = max(t_jitter-dT)

        print ("Time sampling jitter: ",str(max_jitter))

        self.rate_in = 1.0/dT

        self.RecLength_textbox.setText('{:.2f}'.format(self.time_vec[-1]))
        self.Sampling_textbox.setText('{:.2f}'.format(self.rate_in))

        c0.clear()

    def on_exit(self):
        win.close()
        """app.closeAllWindows()"""
        qApp.quit
        sys.exit()

    def level_crossing(self, ):
        self.trigger_condition = 1
        self.statusBar().showMessage('Level Crossing Trigger')

    def local_extremum(self):
        self.trigger_condition = 2
        self.statusBar().showMessage('Local Extremum Trigger')

    def positive_point(self):
        self.trigger_condition = 3
        self.statusBar().showMessage('Positive Point Trigger')

    def on_draw(self):

        """ Re-Analyses the data and redraws the figure"""

        global RD_sequence
        global counter
        global c0
        global p0
        global win

        """get zoom on time series region"""
        self.time_zoom = limits.getRegion()
        lower_limit = (np.abs(self.time_vec-self.time_zoom[0])).argmin()
        upper_limit = (np.abs(self.time_vec-self.time_zoom[1])).argmin()

        """get RD order from textbox"""
        str = self.RDTextbox.text()
        RD_order = int(str)
        """get Decimation factor from textbox"""
        str = self.spinBox.text()

        dec_factor = self.spinBox.value()
        rate_down = self.rate_in/dec_factor

        """ design a highpass filter"""

        num_highpass_s = [1.0, 0]
        den_highpass_s = [1.0, cutoff_highpass]

        [num_highpass, den_highpass] = \
            sp.signal.bilinear(num_highpass_s, den_highpass_s, self.rate_in)

        RD_sequence = np.zeros(RD_order)
        RD_time = np.arange(0, RD_order)
        RD_freq = np.linspace(0, int(rate_down/2), int(RD_order/2))

        self.statusBar().showMessage('Calculating')

        if self.activeaxis == 'x':
             raw_data = self.x_vec
        elif self.activeaxis == 'y':
             raw_data = self.y_vec
        elif self.activeaxis == 'z':
             raw_data = self.z_vec

        trigger_level = self.slider.value() / 100.0 * np.max(raw_data)

        raw_data_plot = raw_data
        raw_data = raw_data[lower_limit:upper_limit]

        if self.checkbox.isChecked():
            """ kill an offset"""
            compensated_data = sp.signal.lfilter(num_highpass, den_highpass,raw_data)
        else:
            compensated_data = raw_data

        if dec_factor != 1:
            """ Downsampling"""
            decimated_data = decimate_filter.decimate_filter(
                compensated_data, dec_factor)
        else:
            decimated_data = compensated_data

        data_time = np.arange(0, float(len(decimated_data))/rate_down, 1.0/rate_down)

        self.statusBar().showMessage('RD Estimation')

        """ Estimate the RD signature"""

        [RD_sequence, counter_out, trigger_event_time, trigger_event_amp] = \
             rd_estim.RD_frame(decimated_data, trigger_level, mirror, self.trigger_condition,\
                           RD_order, counter, RD_sequence)

        """ Calculate spectrum"""
        #Spectrum of the Random Decrement Signature
        mean_sequence = np.var(decimated_data)
        RD_spectrum = sp.fft(RD_sequence) * \
            (2*np.pi*mean_sequence)/trigger_level
        RD_spectrum = RD_spectrum[0:int(RD_order/2)]

        """ pyqtgraph variante"""

        win.setWindowTitle(self.windowtitle)
        c0.clear()
        c0 = p0.plot(self.time_vec, raw_data_plot)
        p1.clear()
        p1.plot(data_time, decimated_data, pen='g', name="Decimated Data ")
        p1.plot(trigger_event_time/rate_down, trigger_event_amp,
                pen=None, symbol='t', symbolPen=None, symbolBrush='r')
        p2.clear()
        p2.plot(RD_time, RD_sequence, pen='g')
        p3.clear()
        p3.plot(RD_freq, abs(RD_spectrum), pen='g')

        self.statusBar().showMessage('Ready')

    def on_radio_button_toggled(self):
            radiobutton = self.sender()

            if radiobutton.isChecked():
                self.activeaxis = radiobutton.axis

    def create_main_frame(self):
        self.main_frame = QWidget()

        # Create the navigation toolbar, tied to the canvas
        # Other GUI controls
        
        """RD order"""
        RDTextbox_label = QLabel('RD order:')
        self.RDTextbox = QLineEdit()
        self.RDTextbox.setMinimumWidth(0.5)

        RecLength_label = QLabel('Record length [s]:')
        self.RecLength_textbox = QLineEdit()
        self.RecLength_textbox.setMinimumWidth(0.5)
        self.RecLength_textbox.setReadOnly(True)

        spinBox_label = QLabel('Decimation')
        self.spinBox = QSpinBox()
        self.spinBox.setGeometry(QRect(440, 260, 71, 22))
        self.spinBox.setMinimum(1)
        self.spinBox.setMaximum(10)
        self.spinBox.setSingleStep(1)
        self.spinBox.setObjectName("spinBox")

        Sampling_label = QLabel('Sampling Frequency [Hz]')
        self.Sampling_textbox = QLineEdit()
        self.Sampling_textbox.setMinimumWidth(0.5)
        self.Sampling_textbox.setReadOnly = (True)

        self.draw_button = QPushButton("&GO!")
        self.draw_button.clicked.connect(self.on_draw)

        slider_label = QLabel('Trigger Level (%):')
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(1, 100)
        self.slider.setValue(25)
        self.slider.setTracking(True)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.slider.setMinimumWidth(0.8)

        checkbox_label = QLabel('Highpass')
        self.checkbox = QCheckBox()
        self.checkbox.toggle()

        #
        # Layout with box sizers
        #

        control_box = QHBoxLayout()

        grid_box = QGridLayout()

        radiobutton = QRadioButton("x")
        radiobutton.axis = "x"
        radiobutton.toggled.connect(self.on_radio_button_toggled)
        grid_box.addWidget(radiobutton, 0, 0)

        radiobutton = QRadioButton("y")
        radiobutton.axis = "y"
        radiobutton.toggled.connect(self.on_radio_button_toggled)
        grid_box.addWidget(radiobutton, 1, 0)

        radiobutton = QRadioButton("z")
        radiobutton.setChecked(True)
        radiobutton.axis = "z"
        radiobutton.toggled.connect(self.on_radio_button_toggled)
        grid_box.addWidget(radiobutton, 2, 0)

        """ Infoboxes Input Data"""
        grid_box.addWidget(Sampling_label, 0,1)
        grid_box.addWidget(self.Sampling_textbox, 1,1)

        grid_box.addWidget(RecLength_label, 0,2)
        grid_box.addWidget(self.RecLength_textbox, 1,2)

        control_box.addLayout(grid_box)

        """ Signal Processing"""
        grid_box.addWidget(checkbox_label, 0,3)
        grid_box.addWidget(self.checkbox, 1,3)

        grid_box.addWidget(spinBox_label, 0,4)
        grid_box.addWidget(self.spinBox, 1,4)

        """ RD Estimation"""
        grid_box.addWidget(slider_label, 0,5)
        grid_box.addWidget(self.slider, 1,5)

        grid_box.addWidget(RDTextbox_label, 0,6)
        grid_box.addWidget(self.RDTextbox, 1,6)

        control_box.addWidget(self.draw_button)

        vbox = QVBoxLayout()
        vbox.addLayout(control_box)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):

        self.statusBar().showMessage('RD Autocorrelation Estimator')

    def create_menu(self):
        """ File Menu"""
        """Exit"""
        exitAction = QAction(QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.on_exit)
        """Import a CSV file from Phyphox"""
        importAction_PP = QAction(QIcon('PP_icon.png'), '&Import from PP', self)
        importAction_PP.setStatusTip('Import CSV file from Phyphox')
        importAction_PP.triggered.connect(self.import_csv_PP)

        """estimate RD signature and draw"""
        importAction_PP.triggered.connect(self.on_draw)

        """Import a CSV file from OA Datalogger"""
        importAction_OA = QAction(QIcon('openadaptronik_icon.png'), '&Import from OA', self)
        importAction_OA.setStatusTip('Import CSV file from OA Datalogger')
        importAction_OA.triggered.connect(self.import_csv_OA)

        """estimate RD signature and draw"""
        importAction_OA.triggered.connect(self.on_draw)

        """Save a plot"""
        SaveFileAction = QAction(QIcon('filesave.png'), '&Save plot', self)
        SaveFileAction.setShortcut('Ctrl+S')
        SaveFileAction.setStatusTip('Save the plot')
        SaveFileAction.triggered.connect(self.save_plot)

        """ Trigger Menu"""
        TriggerLevelCrossingAction = QAction('&Level Crossing', self)
        TriggerLevelCrossingAction.setStatusTip(
            'Set Trigger to Level Crossing')
        TriggerLevelCrossingAction.triggered.connect(self.level_crossing)

        TriggerLocalExtremumAction = QAction('&Local Extremum', self)
        TriggerLocalExtremumAction.setStatusTip(
            'Set Trigger to Local Extremum')
        TriggerLocalExtremumAction.triggered.connect(self.local_extremum)

        TriggerPositivePointAction = QAction('&PositivePoint', self)
        TriggerPositivePointAction.setStatusTip(
            'Set Trigger to Positive Point')
        TriggerPositivePointAction.triggered.connect(self.positive_point)

        menubar = self.menuBar()

        self.fileMenu = menubar.addMenu('&File')
        self.fileMenu.addAction(exitAction)
        self.fileMenu.addAction(importAction_OA)
        self.fileMenu.addAction(importAction_PP)
        self.fileMenu.addAction(SaveFileAction)

        self.triggerMenu = menubar.addMenu("&Trigger")
        self.triggerMenu.addAction(TriggerLevelCrossingAction)
        self.triggerMenu.addAction(TriggerLocalExtremumAction)
        self.triggerMenu.addAction(TriggerPositivePointAction)


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()




if __name__ == "__main__":
    main()
