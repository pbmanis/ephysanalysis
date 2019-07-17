#!/usr/bin/env python3

"""
Brige balance tool
Version 0.1

Graphical interface
Part of Ephysanalysis package

Usage:
notch
"""

import os
import sys
import argparse
from pathlib import Path
import pathlib
import numpy as np
import ephysanalysis as EP
import pylibrary.fileselector as FS
import pyqtgraph as pg
from pyqtgraph.parametertree import Parameter, ParameterTree
import mapanalysistools.digital_filters as FILT

# datadir = '/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107Ai32Het'
# dbfile = 'NF107Ai32Het_bcorr2.pkl'

class TraceAnalyzer(pg.QtGui.QWidget):
    def __init__(self):
        self.datadir = '/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107Ai32Het'
        self.AR = EP.acq4read.Acq4Read()  # make our own private cersion of the analysis and reader
        self.SP = EP.SpikeAnalysis.SpikeAnalysis()
        self.RM = EP.RmTauAnalysis.RmTauAnalysis()
        self.LPF = 5000.
        self.HPF = 0.
        self.notch_freqs = [60., 120., 180., 240.]
        self.notch_Q = 30.
        self.curves = []
        self.AR.traces = None
        self.maxT = 0.6
        
        self.n_adjusted = 0
    
    def getProtocolDir(self):
        sel = FS.FileSelector(dialogtype='dir', startingdir=self.datadir)
        print(sel.fileName)
        self.clampfiles = []
        if sel.fileName is not None:
            self.pdirs = Path(sel.fileName).glob('**/MultiClamp1.ma')
            for p in self.pdirs:
                self.clampfiles.append(p)
                # print(p)
        self.w1.slider.setValue(1)
        self.w1.slider.setMinimum(1)
        self.w1.slider.setMaximum(len(self.clampfiles)+1)
        self.protocolPath = sel.fileName
        print('protocolpath: ', sel.fileName)
        self.updateTraces()

    def setProtocol(self, date, sliceno, cellno, protocolName):
        # create an IV protocol path:
        self.newbr = 0.
        self.protocolBridge = 0.
        self.date = date
        self.slice = sliceno
        self.cell = cellno
        if not '_' in date:
            self.date = date+'_000'
        if isinstance(sliceno, int):
            self.slice = 'slice_{0:03d}'.format(sliceno)
        if isinstance(cellno, int):
            self.cell = 'cell_{0:03d}'.format(cellno)
        self.protocolName = protocolName
        self.protocolPath = Path(self.datadir, self.date, self.slice, self.cell, self.protocolName)
        self.protocolKey = Path(self.date, self.slice, self.cell, self.protocolName)
        if not self.protocolPath.is_dir():
            print('dir not found: ', str(self.protocolPath))
            exit(1)

    def updateTraces(self):
        self.AR.setProtocol(self.protocolPath)  # define the protocol path where the data is
        if self.AR.getData():  # get that data.
            self.update_traces()

    def quit(self):
        exit(0)

    def getProtocols(self):
        thisdata = self.df.index[(self.df['date'] == self.date) &
                                (self.df['slice_slice'] == self.slice) &
                                (self.df['cell_cell'] == self.cell)].tolist()
        if len(thisdata) > 1:
            raise ValueError('Search for data resulted in more than one entry!')
        ivprots = self.df.iloc[thisdata]['IV'].values[0]  # all the protocols in the dict
        return thisdata, ivprots
        
    def getProtocol(self, protocolName):
        thisdata, ivprots = self.getIVProtocols()
        if protocolName not in ivprots.keys():
            return None
        else:
            return ivprots[protocolName]

    def set_window(self, parent=None):
        super(TraceAnalyzer, self).__init__(parent=parent)
        self.win = pg.GraphicsWindow(title="TraceAnalyzer")
        layout = pg.QtGui.QGridLayout()
        layout.setSpacing(10)
        self.win.setLayout(layout)
        self.win.resize(1280, 800)
        
        self.buttons = pg.QtGui.QGridLayout()
        layout.addLayout(self.buttons, 0, 0, 7, 1)
        self.b3 = pg.QtGui.QPushButton("Get Protocol Directory")
        self.buttons.addWidget(self.b3)
        self.b3.clicked.connect(self.getProtocolDir)
        
        self.notch_button = pg.QtGui.QDoubleSpinBox()
        self.notch_button.setValue(60.)
        self.notch_button.setMinimum(10.)
        self.notch_button.setMaximum(5000.)
        self.buttons.addWidget(self.notch_button)
        self.notch_button.valueChanged.connect(self.update_traces)
        self.notch_check = pg.QtGui.QCheckBox("Notch")
        self.notch_check.stateChanged.connect(self.update_traces)
        self.buttons.addWidget(self.notch_check)
        
        
        self.lpf_button = pg.QtGui.QDoubleSpinBox()
        self.lpf_button.setValue(5000.)
        self.lpf_button.setMinimum(200.)
        self.lpf_button.setMaximum(10000.)
        self.buttons.addWidget(self.lpf_button)
        self.lpf_button.valueChanged.connect(self.update_traces)
        self.buttons.addWidget(self.lpf_button)
        self.lpf_check = pg.QtGui.QCheckBox("LPF")
        self.lpf_check.stateChanged.connect(self.update_traces)
        self.buttons.addWidget(self.lpf_check)
        
        spacerItem1 = pg.QtGui.QSpacerItem(0, 400, pg.QtGui.QSizePolicy.Expanding, pg.QtGui.QSizePolicy.Minimum)
        self.buttons.addItem(spacerItem1)
        
        self.b2 = pg.QtGui.QPushButton("Quit")
        self.buttons.addWidget(self.b2)
        self.b2.clicked.connect(self.quit)

        # spacerItem = pg.QtGui.QSpacerItem(0, 10, pg.QtGui.QSizePolicy.Expanding, pg.QtGui.QSizePolicy.Minimum)
        # self.buttons.addItem(spacerItem)

        self.dataplot = pg.PlotWidget()
        layout.addWidget(self.dataplot, 0, 1, 6, 7)
        self.w1 = Slider(0, 250, scalar=1.)
        self.w1.slider.valueChanged.connect(self.update_traces)
        layout.addWidget(self.w1, 7, 1, 1, 7)
        layout.setColumnStretch(0, 1)  # reduce width of LHS column of buttons
        layout.setColumnStretch(1, 7)  # and stretch out the data dispaly
        
    def update_traces(self):
        if self.AR.traces is None:
            return
        self.current_trace = int(self.w1.x)
        self.dataplot.setTitle(f'Trace: {self.current_trace:d}')
        for c in self.curves:
            c.clear()
        self.curves = []
        notchfr = self.notch_button.value()
        i = self.current_trace
        if (i > self.AR.traces.shape[0]):
            self.dataplot.setTitle(f'Trace > Max traces: {self.AR.traces.shape[0]:d}')
            return
        imax = int(self.maxT*self.AR.sample_rate[0])
        print('imax: ', imax)
        mod_data = self.AR.traces[i,:].copy()[:imax]
        if self.notch_check.checkState():
            mod_data =  FILT.NotchFilterZP(mod_data, notchf=[notchfr], Q=self.notch_Q,
                               QScale=False, samplefreq=self.AR.sample_rate[0])
        if self.lpf_check.checkState():
            mod_data = FILT.SignalFilter_LPFBessel(mod_data, self.lpf_button.value(), samplefreq=self.AR.sample_rate[0], NPole=8)
        self.curves.append(self.dataplot.plot(self.AR.time_base[:imax]*1e3,
                            # self.AR.traces[i,:],
                            mod_data, 
                           # FILT.NotchFilterZP(self.AR.traces[i,:], notchf=[notchfr], Q=self.notch_Q,
 #                               QScale=False, samplefreq=self.AR.sample_rate[0]),
                           pen=pg.intColor(1)))

class FloatSlider(pg.QtGui.QSlider):
    def __init__(self, parent, decimals=3, *args, **kargs):
        super(FloatSlider, self).__init__(parent, *args, **kargs)
        self._multi = 10 ** decimals
        self.setMinimum(self.minimum())
        self.setMaximum(self.maximum())
 
    def value(self):
        return float(super(FloatSlider, self).value()) / self._multi
 
    def setMinimum(self, value):
        self.min_val = value
        return super(FloatSlider, self).setMinimum(value * self._multi)
 
    def setMaximum(self, value):
        self.max_val = value
        return super(FloatSlider, self).setMaximum(value * self._multi)
 
    def setValue(self, value):
        super(FloatSlider, self).setValue(int((value-self.min_val) * self._multi))


class Slider(pg.QtGui.QWidget):
    def __init__(self, minimum, maximum, scalar=1.0, parent=None):
        super(Slider, self).__init__(parent=parent)
        self.verticalLayout = pg.QtGui.QVBoxLayout(self)
        self.label = pg.QtGui.QLabel(self)
        self.verticalLayout.addWidget(self.label, alignment=pg.QtCore.Qt.AlignHCenter)
        self.horizontalLayout = pg.QtGui.QHBoxLayout()
        spacerItem = pg.QtGui.QSpacerItem(0, 20, pg.QtGui.QSizePolicy.Expanding, pg.QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.slider = FloatSlider(self, decimals=2)
        self.slider.setOrientation(pg.QtCore.Qt.Horizontal)
        self.horizontalLayout.addWidget(self.slider)
        spacerItem1 = pg.QtGui.QSpacerItem(0, 20, pg.QtGui.QSizePolicy.Expanding, pg.QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.resize(self.sizeHint())

        self.minimum = minimum*scalar
        self.maximum = maximum*scalar
        self.scalar = scalar
        self.slider.setMinimum(self.minimum)
        self.slider.setMaximum(self.maximum)
        #self.slider.setRange(self.minimum, self.maximum)
        self.slider.valueChanged.connect(self.setLabelValue)
        self.setLabelValue(self.slider.value())

    def setLabelValue(self, value):
        self.x = int((self.minimum + (float(value) / (self.slider.maximum() - self.slider.minimum())) * (
        self.maximum - self.minimum)) /self.scalar)
        self.label.setText(f"{self.x:4d}")
    
    def getPosValue(self, x):
        return int((x-self.minimum)*(self.slider.maximum() - self.slider.minimum())/(self.maximum - self.minimum))

def main():

    TA = TraceAnalyzer()

    app = pg.mkQApp()
    app.aboutToQuit.connect(TA.quit)  # prevent python exception when closing window with system control
    TA.set_window()
    # TA.show()

    if sys.flags.interactive == 0:
     app.exec_()


if __name__ == '__main__':
    main()

