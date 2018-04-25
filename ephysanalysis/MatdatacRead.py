# -*- coding: utf-8 -*-
from collections import OrderedDict
import os, re, pprint
import sys
import user

import itertools
import pickle
import pyqtgraph as pg

from scipy.io import matlab
import numpy as np

import matplotlib.pyplot as mpl
from pyqtgraph import metaarray
# from acq4.util import DataManager
# import acq4.analysis.tools.Utility as Utility  # pbm's utilities...
# import acq4.analysis.tools.Fitting as Fitting  # pbm's fitting stuff...
import Utility
import Fitting
# from acq4.analysis.dataModels import PatchEPhys
# DM = DataManager
import pylibrary.PlotHelpers as PH
try:
    from PyQt4 import QtCore, QtGui
except:
    try:
        from PyQt5 import QtCore, QtGui
    except:
        raise Importerror('Cannot find Pyqt5 or pyqt4')



import SpikeAnalysis
import RmTauAnalysis

class DatacBrowser(QtGui.QSplitter):
    def __init__(self, path):
        self.path = path
        
        QtGui.QSplitter.__init__(self, QtCore.Qt.Horizontal)
        self.split1 = QtGui.QSplitter(QtCore.Qt.Vertical, self)
        self.addWidget(self.split1)
        self.split2 = QtGui.QSplitter(QtCore.Qt.Vertical, self)
        self.addWidget(self.split2)
        
        self.tree = DirTreeWidget(parent=self, baseDir=path)
        self.split1.addWidget(self.tree)
        self.tree.sigSelectionChanged.connect(self.loadFile)
        
        self.blockList = QtGui.QListWidget()
        self.split1.addWidget(self.blockList)
        self.blockList.currentItemChanged.connect(self.blockSelected)
        
        self.desc = QtGui.QPlainTextEdit()
        self.split1.addWidget(self.desc)
        
        self.plt1 = pg.PlotWidget(labels={'bottom': ('Time', 's')})
        self.plt2 = pg.PlotWidget()
        self.plt3 = pg.PlotWidget()
#        self.plt2.setXLink(self.plt1)
#        self.plt3.setXLink(self.plt1)
        for p in [self.plt1, self.plt2, self.plt3]:
#            p.disableAutoRange()
            self.split2.addWidget(p)
#            p.setClipToView(True)
            #p.setDownsampling(auto=True, mode='peak')
            
        self.split2.setStretchFactor(0, 50)
        self.split2.setStretchFactor(1, 10)
        self.split2.setStretchFactor(2, 10)
        self.plotItems = []
        self.clear()

    def loadFile(self, selfile=None):
        """
        load a file. If selfile is None, then we try to get the file from the browser tree
        else, we use sel as the name of the file to be loaded
        """
        if selfile is None:
            sel = self.tree.selectedFile()
            if sel is None or not sel.endswith('.mat'):
                return
            self.clear()
        else:
            if isinstance(selfile, DirTreeWidget):  # need to parse from the tree
                try:
                    selfile = os.path.join(self.path, selfile.selectedItems()[0].text(0))
                except:
                    print 'no file selected'
                    return 
            if not os.path.isfile(selfile) or not selfile.endswith('.mat'):
                print 'File not found: %s' % selfile
                return
            self.clear()
            sel = selfile
#        print 'SEL: ', sel
            
        with pg.BusyCursor():
#            print 'reading file'
            self.loaded = DatacFile(sel)
            self.desc.setPlainText(self.loaded.summary())
        
        for i in self.loaded.items:
#            print 'items: %d', i
            item = QtGui.QListWidgetItem(str(i.type))
            item.datac_item = i
            self.blockList.addItem(item)
#        print 'blockList: ', self.blockList

    def clear(self):
        self.loaded = None
        self.blockList.clear()
        self.clearBlock()

    def clearBlock(self):
        for i in self.plotItems:
            i.scene().removeItem(i)
        self.plotItems = []
        self.desc.setPlainText('')

    def getDataMode(self):
        datamode = self.currentBlock.dfile['Data_Mode']['v']
        #print "current block data mode : ", datamode
        if datamode in ['cc', 'ic', 'CC', 'IC']:
            datamode = 'cc'
            chmap = [1, 0, 2]
        elif datamode in ['vc', 'VC']:
            datamode = 'vc'
            chmap = [0, 1, 2]
        else:
            raise ValueError('Data mode <%s> is not known' % datamode)
        return datamode, chmap

    def blockSelected(self, item):
        self.clearBlock()
        if item is None:
            return
        self.currentBlock = item.datac_item
        datamode = self.getDataMode()
        
        datamode = self.currentBlock.dfile['Data_Mode']['v']
        #print "current block data mode : ", datamode
        if datamode in ['cc', 'ic', 'CC', 'IC']:
            datamode = 'cc'
            chmap = [1, 0, 2]
        elif datamode in ['vc', 'VC']:
            datamode = 'vc'
            chmap = [0, 1, 2]
        else:
            raise ValueError('Data mode <%s> is not known' % datamode)
        self.Clamps = GetClamps(self, self.path)
        self.Clamps.getClampData(chmap, item)
        
        self.desc.setPlainText(self.currentBlock.summary())
        if isinstance(self.currentBlock, Block):
            with pg.BusyCursor():
                rate, recs = self.currentBlock.data()
                nchannels = len(self.currentBlock.dfile['Channels']['v'])
                dt = nchannels / rate

                for i, x in enumerate(recs):
                    pen = (i, len(recs)*1.5)
                    p1 = self.plt1.plot(np.arange(len(x[chmap[0]]))*dt, x[chmap[0]]*1e-3, pen=pen)
                    p2 = self.plt2.plot(np.arange(len(x[chmap[1]]))*dt, x[chmap[1]]*1e-12, pen=pen)
                    cmd = np.interp(np.arange(len(x[chmap[0]]))*dt, np.arange(0, len(x[chmap[2]]))/1000., x[chmap[2]])
                    # print 'len(x[chmap[2]): ', len(x[chmap[2]])
                    # print 'len(cmd): ', len(cmd)
                    # print 'cmd min max: ', np.min(cmd), np.max(cmd)
                    p3 = self.plt3.plot(np.arange(len(x[chmap[0]]))*dt, cmd*1e-12, pen=pen)
                    self.plotItems.extend([p1, p2, p3])

        print np.max(len(x[chmap[0]]))*dt

        self.plt1.setXRange(0, np.max(len(x[chmap[0]]))*dt)
        self.plt2.setXRange(0, np.max(len(x[chmap[1]]))*dt)
        self.plt3.setXRange(0, np.max(len(x[chmap[0]]))*dt)
        self.plt1.setYRange(-0.120, 0.060)
        self.plt2.setYRange(-1e-9, 1e-9)
        self.plt3.setYRange(-1e-9, 1e-9)
        
    
class GetClamps(DatacBrowser):
    
    def __init__(self, datac, path):
        super(GetClamps, self).__init__(path)
        self.datac = datac


    def getClampData(self, chmap=None, item=None):
        """
        create a Clamp structure for use in SpikeAnalysis and RMTauAnalysis from acq4.
        Fills in the fields that are returned by PatchEPhys getClamps:
        clampInfo['dirs]
        clampInfo['missingData']
        self.time_base
        self.values
        self.traceStartTimes = np.zeros(0)
        self.sequence
        self.clampValues (sequence)
        self.nclamp = len(self.clmapVlaues
        self.repc
        self.nrepc
        self.data_mode
        self.model_mode = False
        self.command_scale_factor
        self.command_units
        self.devicesUsed
        self.clampDevices
        self.holding
        self.clampState
        self.sample_interval
        self.RSeriesUncomp
        self.amplifeirSettings['WCCompValid', 'WCEmabled', 'CompEnabled', 'WCSeriesResistance']
        self.cmd_wave
        self.commandLevels (np.array(self.values))
        self.traces = MetaArray(traces, info=info)
        self.tstart
        self.tdur
        self.tend
        self.spikecount = np.zeros(len...) if in vcmode.
        
        Info from an example data file:
        [{'name': 'Channel', 'cols': [{'units': 'A', 'name': 'Command'}, {'units': 'V', 'name': 'primary'}, {'units': 'A', 'name': 'secondary'}]},
        {'units': 's', 'values': array([ 0.00000000e+00, 2.50000000e-05, 5.00000000e-05, ..., 6.99925000e-01, 6.99950000e-01, 6.99975000e-01]),
        'name': 'Time'}, {'ClampState': {'primaryGain': 10.0, 'ClampParams': {'OutputZeroEnable': 0, 'PipetteOffset': 0.05197399854660034,
        'Holding': -1.525747063413352e-11, 'PrimarySignalHPF': 0.0, 'BridgeBalResist': 20757020.0, 'PrimarySignalLPF': 20000.0, 'RsCompBandwidth':
        8.413395979806202e-42, 'WholeCellCompResist': 8.413395979806202e-42, 'WholeCellCompEnable': 6004, 'LeakSubResist': 8.413395979806202e-42,
        'HoldingEnable': 1, 'FastCompTau': 8.413395979806202e-42, 'SlowCompCap': 8.413395979806202e-42, 'WholeCellCompCap': 8.413395979806202e-42,
        'LeakSubEnable': 6004, 'NeutralizationCap': 1.9578947837994853e-12, 'BridgeBalEnable': 1, 'RsCompCorrection': 8.413395979806202e-42,
        'NeutralizationEnable': 1, 'RsCompEnable': 6004, 'OutputZeroAmplitude': -0.0009990156395360827, 'FastCompCap': 8.413395979806202e-42,
        'SlowCompTau': 8.413395979806202e-42}, 'secondarySignal': 'Command Current', 'secondaryGain': 1.0, 'secondaryScaleFactor': 2e-09,
        'primarySignal': 'Membrane Potential', 'extCmdScale': 4e-10, 'mode': 'IC', 'holding': 0.0, 'primaryUnits': 'V', 'LPFCutoff': 20000.0,
        'secondaryUnits': 'A', 'primaryScaleFactor': 0.1, 'membraneCapacitance': 0.0}, 'Protocol': {'recordState': True, 'secondary': None,
        'primary': None, 'mode': 'IC'}, 'DAQ': {'command': {'numPts': 28000, 'rate': 40000.0, 'type': 'ao', 'startTime': 1296241556.7347913},
        'primary': {'numPts': 28000, 'rate': 40000.0, 'type': 'ai', 'startTime': 1296241556.7347913}, 'secondary': {'numPts': 28000, 'rate':
        40000.0, 'type': 'ai', 'startTime': 1296241556.7347913}}, 'startTime': 1296241556.7347913}]

        )
        """
        if item is not None:
            protocol = item.text()
        else:
            protocol = 'unknown'
        
        rate, recs = self.datac.currentBlock.data()
        self.dfile = self.datac.currentBlock.dfile
        (ddir, fname) = os.path.split(self.datac.currentBlock.datac.fname)
        points = self.dfile['Points']['v']
        nchannels = len(self.dfile['Channels']['v'])
        dt = nchannels / rate
        data_mode = self.datac.getDataMode()
        # print len(recs)
        # print recs[0]
        self.data_mode = data_mode[0]
        self.time_base = np.arange(points)*dt
        datac_traces = np.zeros((len(recs), nchannels+1, points))
        for i in range(len(recs)):
            for j in range(nchannels+1):
                if j == 2:
                    # print len(recs[i][j])
                    # print len(np.arange(0, len(recs[i][j])))
                    cmd = np.interp(self.time_base, np.arange(0, len(recs[i][j]))/1000., recs[i][j])
                    datac_traces[i, j, :] = cmd
                else:
                    # print 'i,j', i, j
                    datac_traces[i, j, :] = recs[i][j]
                    if j == 1:
                        datac_traces[i, j, :] *= 1e-3 # convert to V from mV
        #np.array([[x[chmap[i]] for x in recs] for i in range(len(chmap))])  # remap channels
        #self.values = np.array([np.mean(datac_traces[1,i]) for i in range(len(recs))])
        # print self.values
        self.repc = 1
        self.nrepc = 1
        self.model_mode = False
        self.command_scale_factor = 1
        self.command_units = 'A'
        self.devicesUsed = None
        self.clampDevices = None
        self.holding = 0.
        self.amplfierSettings = {'WCCompValid': False, 'WCEnabled': False, 
                'CompEnabled': False, 'WCSeriesResistance': 0.}
        self.clampState = None
        self.sample_interval = dt
        self.RSeriesUncomp = 0.
        self.cmd_wave = np.squeeze(datac_traces[:, 2, :])
        self.protoTimes = {'drugtestiv': [0.21, 0.51], 'ap-iv2': [0.01, 0.5]}
        self.tstart = 0.01
        self.tdur = 0.500
        if protocol in self.protoTimes:
            self.tstart = self.protoTimes[protocol][0]
            self.tdur = self.protoTimes[protocol][1]
            
        self.tend = self.tstart + self.tdur
        self.traces = np.array(datac_traces)
        cmds = np.squeeze(self.traces[:,0,:])
        t0 = int(self.tstart/dt)
        t1 = int(self.tend/dt)
        self.values = np.mean(cmds[:, t0:t1]*1e-12, axis=1)  # express values in amps
        self.commandLevels = self.values
        info = [{'units': 'A', 'values': self.values, 'name': 'Command'},
                    {'name': 'Time', 'units': 's', 'values': self.time_base},
                    {'ClampState': {'primaryGain': 10.0, 'ClampParams': {'OutputZeroEnable': 0, 'PipetteOffset': 0.0,
                        'Holding': 0, 'PrimarySignalHPF': 0.0, 'BridgeBalResist': 0.0, 'PrimarySignalLPF': 20000.0, 'RsCompBandwidth':
                        0.0, 'WholeCellCompResist': 0.0, 'WholeCellCompEnable': 6004, 'LeakSubResist': 0.0,
                        'HoldingEnable': 1, 'FastCompTau': 0.0, 'SlowCompCap': 0.0, 'WholeCellCompCap': 0.,
                        'LeakSubEnable': 6004, 'NeutralizationCap': 0., 'BridgeBalEnable': 0, 'RsCompCorrection': 0.0,
                        'NeutralizationEnable': 1, 'RsCompEnable': 6004, 'OutputZeroAmplitude': 0., 'FastCompCap': 0.,
                        'SlowCompTau': 0.0}, 'secondarySignal': 'Command Current', 'secondaryGain': 1.0, 'secondaryScaleFactor': 2e-09,
                        'primarySignal': 'Membrane Potential', 'extCmdScale': 4e-10, 'mode': 'IC', 'holding': 0.0, 'primaryUnits': 'V', 'LPFCutoff': 20000.0,
                        'secondaryUnits': 'A', 'primaryScaleFactor': 0.1, 'membraneCapacitance': 0.0}, 'Protocol': {'recordState': True, 'secondary': None,
                        'primary': None, 'mode': 'IC'}, 'DAQ': {'command': {'numPts': 10000, 'rate': 10000.0, 'type': 'ao', 'startTime': 0.},
                        'primary': {'numPts': 10000, 'rate': 10000.0, 'type': 'ai', 'startTime': 0.}, 'secondary': {'numPts': 1000, 'rate':
                        10000.0, 'type': 'ai', 'startTime': 0.}}, 'startTime': 0.}]
        self.traces = np.squeeze(self.traces[:,1,:])
        if fname == '06nov09d.mat':
            print self.traces.shape
            print np.max(self.traces[0,:])
            print np.min(self.traces[0,:])
            print np.min(datac_traces[0, 0, :])
            for i in range(self.traces.shape[0]):
                print i
                self.traces[i,:] = self.traces[i,:]-50e6*1e-12*datac_traces[i, 0, :]  # IR drop
            print np.min(self.traces[0,:])
            print 'bridge corrected!!!!!'

        self.traces = metaarray.MetaArray(self.traces, info=info)
        self.spikecount = np.zeros(len(recs))
        self.rgnrmp = [0, 0.005]
        print('getclamps got it all')


class DirTreeWidget(QtGui.QTreeWidget):

    sigSelectionChanged = QtCore.Signal(object)

    def __init__(self, parent=None, baseDir=None):
        QtGui.QTreeWidget.__init__(self, parent)
        self.baseDir = os.path.normpath(os.path.abspath(baseDir))
        self.currentDir = None
        self.currentItemChanged.connect(self.selectionChanged)
        self.buildTree()

    def selectionChanged(self, item=None, _=None):
        self.sigSelectionChanged.emit(self)

    def selectedFile(self):
        """Return the name of the currently selected file.
        If no items are selected, return None.
        If multiple items are selected, raise an exception."""
        items = self.selectedItems()
        if len(items) == 0:
            return None
        if len(items) > 1:
            raise Exception('Multiple items selected. Use selectedFiles instead.')
        return items[0].fname

    def buildTree(self, root=None):
        if root is None:
            root = self.invisibleRootItem()
            path = self.baseDir
        else:
            path = root.fname
            
        ls = self.ls(path)
        for fname in ls:
            item = QtGui.QTreeWidgetItem([os.path.split(fname)[1]])
            root.addChild(item)
            item.fname = fname
            if os.path.isdir(fname):
                self.buildTree(item)
            
    def ls(self, path):
        """List files ordered by date indicated in file names generated by acq3
        (like 05apr11a.mat).
        """
        ls = [os.path.join(path, ch) for ch in os.listdir(path)]
        def fileval(path):
            if os.path.isdir(path):
                return (0, 0, 0, 0)
            m = re.match('(\d\d)([a-z]{3})(\d\d)([a-z]).mat', os.path.split(path)[1])
            if m is None:
                return (1e12, 0, 0, 0)
            day, mon, year, cell = m.groups()
            mon = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'].index(mon)
            return (int(year), mon, int(day), cell)
        ls.sort(key=fileval)
        return ls
                

class DatacFile(object):
    def __init__(self, fname):
        self.fname = fname
        self.data = matlab.loadmat(fname)
        
        # read block records
        self.items = []
        blockid = None
        blockrecs = []
        for rec in self.data['Index'][0]:
            typ = rec['type'][0]
            if typ == 'HFILE':
                self.items.append(HFile(self, rec))
            elif typ == 'NOTE':
                self.items.append(Note(self, rec))
            elif typ in ('SFILE', 'DFILE', 'DATA'):
                bid = rec['block'][0][0]
                if blockid != bid:
                    if blockid is not None:
                        self.items.append(Block(self, blockrecs))
                        blockrecs = []
                    blockid = bid
                blockrecs.append(rec)
        
        if len(blockrecs) > 0:
            self.items.append(Block(self, blockrecs))
            
    def summary(self):
        summ = 'FILE: ' + self.fname
        typs = set([i.type for i in self.items if i.type not in ('HFILE', 'NOTE')])
        summ += '\nBLOCKS: ' + '  '.join(list(typs))
        
        for i in self.items:
            if i.type == 'HFILE':
                summ += '\n' + i.summary() + '\n'
            if i.type == 'NOTE':
                summ += '\n* ' + i.note 
        return summ


class HFile(object):
    def __init__(self, datac, rec):
        self.datac = datac
        self.rec = rec
        self.data = flatten(self.datac.data[self.rec['MatName'][0]])
        self.type = 'HFILE'
        
    def summary(self):
        exp = self.data['Experiment']['v']
        
        # What does it mean to have multiple values in the Experiment field??
        if isinstance(exp, np.ndarray):
            exp = '\n'.join(exp)
        
        summ = exp + '\n'
        summ += 'Species: %s  Age: %s' % (self.data['Species']['v'],
                                          self.data['Age']['v'])
        return summ

    
class Note(object):
    def __init__(self, datac, rec):
        self.datac = datac
        self.rec = rec
        self.type = 'NOTE'
        self.note = flatten(datac.data[rec['MatName'][0]])['v']
        
    def summary(self):
        return self.note
    

class Block(object):
    def __init__(self, datac, block_recs):
        self.datac = datac
        self._data = None
        self.recs = block_recs
        self.notes = []
        self.type = None
        self.stims = []
        self.recordings = []
        for rec in self.recs:
            typ = rec['type'][0]
            if rec['MatName'][0] not in self.datac.data.keys():
                print 'missing data block: ', rec['MatName'][0]
                continue
            data = flatten(self.datac.data[rec['MatName'][0]])
            if typ == 'DATA':
                self.type = rec['type2'][0]
                self.recordings.append(data)
            elif typ == 'SFILE':
                self.type = rec['type2'][0]
                self.stims.append(data)
            elif typ == 'DFILE':
                self.dfile = data
            else:
                raise ValueError("Unrecognized data type: %s" % typ)
        
    def summary(self, printnow=True):
        summ = '%d data, %d stims\n' % (len(self.recordings[0]), len(self.stims[0]))
        summ += '---- DFILE: ----\n'
        # summ = ''
        # for k in self.dfile.keys():
        #     if k in ['Actual_Rate', 'ChannelList', 'Data_Mode', 'Points', 'Sample_Rate', 'TraceDur']:
        #          fmtstr = '%15s : ' + ('%s\n' % self.dfile[k]['f'])
        #          summ += fmtstr % (k, self.dfile[k]['v'])
        # summ += '\n'
        # summ += '%15s : %3d - %3d' % ('Records', self.dfile['Record'], len(self.recordings[0]))
        summ += pprint.pformat(self.dfile)

        return summ
        
    def data(self):
        """Return (sample rate, [traces]) for this block, where [traces] is
        a list of (primary_channel, secondary_channel, stim) tuples.
        """
        if self._data is not None:
            return self._data
        
        d = self.recordings[0]
        mode = d[0].values()[0]['status']['Data']['mode']
        if isinstance(mode, np.ndarray):
            mode = mode[0]
        x = []
        reps = self.stims[0]['Repeats']['v']
        # print 'reps: ', reps
        for i in range(len(d)):
            data = d[i].values()[0]['data']
            stim = self.stims[0]['waveform']
            # print stim
            # print 'stim shape: ', stim.shape
            # print 'i, reps, index: ', i, reps, i//reps
            if isinstance(stim, np.ndarray):
                # Is this right? when len(stims) < len(recordings), how should
                # they be matched together?
                stim = stim[0]['v2'] # stim[i//reps]['v2']
            else:
                stim = stim['v2']
            
            if mode == 'V-Clamp':
                x.append((pg.gaussianFilter(data[:, 1], sigma=(2,)), data[:, 0], stim))
            else:
                x.append((data[:,0], data[:,1], stim))
                # x.append((pg.gaussianFilter(data[:, 0], sigma=(2,)), data[:, 1], stim))
        
        rate = 1e6 / self.dfile['Sample_Rate']['v']
                
        self._data = rate, x
        return self._data


def flatten(data):
    """Flatten highly nested structures returned by loadmat
    """
    if isinstance(data, np.ndarray):
        while data.ndim > 0 and data.shape[0] == 1:
            data = data[0]
        if data.ndim == 0 and data.dtype.names is not None:
            data = dict(zip(data.dtype.names, data))
        elif data.ndim > 0:
            if data.dtype.kind == 'O':
                for i in ndinds(data.shape):
                    data[i] = flatten(data[i])
            elif data.dtype.names is not None:
                for fname in data.dtype.names:
                    if data.dtype.fields[fname][0].kind == 'O':
                        field = data[fname]
                        for i in ndinds(field.shape):
                            field[i] = flatten(field[i])
                
    if isinstance(data, dict):
        for k,v in data.items():
            data[k] = flatten(v)
    return data


def ndinds(shape):
    ind = np.indices(shape)
    while len(ind.shape) > 2:
        ind = ind.reshape(ind.shape[:-2] + (ind.shape[-2]*ind.shape[-1],))
    return [tuple(x) for x in ind.T]


class IVAnalyzer(DatacBrowser):
    """Browser with additional features for analyzing IVCurves.
    """

    def __init__(self, path, infiles=None, outfile=None):
        DatacBrowser.__init__(self, path)
#        self.lr1 = pg.LinearRegionItem([0, 0.01])  # maybe use one of these to set up the baseline...
#        self.lr2 = pg.LinearRegionItem([0, 0.01])
#        self.plt1.addItem(self.lr1)
#        self.plt3.addItem(self.lr2)
        
        self.split3 = QtGui.QSplitter(QtCore.Qt.Vertical, self)
        self.addWidget(self.split3)
        self.allAnalysis = {}
        self.IV_plot = pg.PlotWidget()
        self.split3.addWidget(self.IV_plot)
        self.fiPlot = pg.PlotWidget()
        self.split3.addWidget(self.fiPlot)
        self.fslPlot = pg.PlotWidget()
        self.split3.addWidget(self.fslPlot)
        
        # self.histogram = pg.PlotWidget()
        # self.histData = self.histogram.plot(fillLevel=0, fillBrush=(200, 0, 0, 255))
        # self.split3.addWidget(self.histogram)
        # Things that IVCurve does that we need to duplicate.
        self.headerPrint = True
        self.analysis_summary={}
        self.ic_modes = ['cc', 'CC', 'IC', 'ic', 'I-Clamp']
        self.vc_modes = ['vc', 'VC', 'V-Clamp']
        self.dataMarkers = {}
        self.fitMarkers = {}
        self.data_plot = self.plt1
        self.colors = ['w', 'g', 'b', 'r', 'y', 'c']
        self.symbols = ['o', 's', 't', 'd', '+']
        self.color_list = itertools.cycle(self.colors)
        self.symbol_list = itertools.cycle(self.symbols)
        self.lrss_show = True
        self.lrpk_show = True
        self.showFISI = True
        self.keep_analysis = False
        self.keep_analysis_count = 0
        self.make_map_symbols()
        
        self.threshold = 0
        self.sub_baseline = False
        self.sub_leak = False

        self.data_template = (
          OrderedDict([('Cell', (14, '{:>12s}')), ('FirstRec', (4, '{:>4d}')), ('Genotype', (12, '{:>12s}')), ('Protocol', (10, '{:>10s}')), ('Species', (12, '{:>12s}')), ('Age', (5, '{:>5s}')), ('Sex', (3, '{:>3s}')), ('Weight', (6, '{:>6s}')),
                       ('Temperature', (10, '{:>10s}')), ('ElapsedTime', (11, '{:>11.2f}')), 
                       ('RMP', (5, '{:>5.1f}')), ('Rin', (5, '{:>5.1f}')), ('Bridge', (5, '{:>5.1f}')),
                       ('tau', (5, '{:>5.1f}')), ('AdaptRatio', (9, '{:>9.3f}')),
                       ('tauh', (5, '{:>5.1f}')), ('Gh', (6, '{:>6.2f}')),
                       ('FiringRate_1p5T', (12, '{:>9.1f}')), 
                       ('AP1_HalfWidth', (13, '{:>13.2f}')), ('AP1_Latency', (11, '{:>11.1f}')), 
                       ('AP2_HalfWidth', (13, '{:>13.2f}')), ('AP2_Latency', (11, '{:>11.1f}')), 
                       ('AHP_Depth', (9, '{:>9.2f}')),
                       ('Fzero', (9, '{:>9.2f}')),
                       ('Ibreak', (9, '{:>9.3f}')),
                       ('F1amp', (9, '{:>9.3f}')),
                       ('F2amp', (9, '{:>9.3f}')),
                       ('Ibreak', (9, '{:>9.3f}')),
                       ('Ibreak1', (12, '{:>12.6f}')),
                       ('Irate1', (12, '{:>12.6f}')),
                       ('Irate2', (12, '{:>12.6f}')),
                       ('Irate3', (12, '{:>12.6f}')),
                       ('Rate0', (12, '{:>12.6f}')),
                       ('FIFitFunc', (12, '{:>12s}')),
                       ('Description', (11, '{:s}')),
                      ]))
        self.datafile = outfile
#        self.lr1.sigRegionChangeFinished.connect(self.updateAnalysis)
#        self.lr2.sigRegionChangeFinished.connect(self.updateAnalysis)
        #self.lr3.sigRegionChangeFinished.connect(self.updateAnalysis)
        if infiles is None:
            return
        for fn in infiles.keys():
            if fn.endswith('.mat'):  # handle the .mat files
                protocol = infiles[fn]
                fnfull = os.path.join(path, fn)
                DatacBrowser.loadFile(self, selfile=fnfull)  # select this file, load the browser lists
                for i in range(self.blockList.count()):
                    p = self.blockList.item(i).text()
                    if p == protocol:
                        self.analysis_summary={}
                        self.blockList.setCurrentItem(self.blockList.item(i))
                        self.show()
            else:  # acq4 files
                pass
                # self.dataModel = PatchEPhys
                # for i in range(len(infiles[fn])):
                #     if infiles[fn][i][0] is None:  # fiel to skip
                #         continue
                #     fnfull = os.path.join(path, fn, infiles[fn][i][0])
                #     dh = DM.getDirHandle(fnfull, create=False)
                #     self.analysis_summary={}
                #     self.loadFileRequested(dh)
                #

    def blockSelected(self, item):
        DatacBrowser.blockSelected(self, item)
        if not isinstance(self.currentBlock, Block):
            return
        if item is None:  # wait for selection of a block 
            return

        self.clearDecorators()
        self.clearFits()
        rate, recs = self.currentBlock.data()
        stimInd = recs[0][2].argmax()
        stimTime = float(stimInd) / rate
        # self.lr1.sigRegionChangeFinished.disconnect(self.updateAnalysis)
        # self.lr2.sigRegionChangeFinished.disconnect(self.updateAnalysis)
        
        # try:
        #     self.lr1.setRegion([stimTime + 0.5e-3, stimTime + 5e-3])
        #     self.lr2.setRegion([stimTime, stimTime+20e-6])
        #
        # finally:
        #     self.lr1.sigRegionChangeFinished.connect(self.updateAnalysis)
        #     self.lr2.sigRegionChangeFinished.connect(self.updateAnalysis)
        (ddir, filename) = os.path.split(self.currentBlock.datac.fname)
        self.analysis_summary['Cell'] = filename
        
        self.analysis_summary['Genotype'] = 'CBA/CaJ'
        if item is not None:
            self.analysis_summary['Protocol'] = item.text()
            self.analysis_summary['EndRec'] = self.currentBlock.recs[2][3][0][0]
            self.analysis_summary['Block'] = self.currentBlock.recs[0][2][0][0]
            self.analysis_summary['FirstRec'] = self.currentBlock.recs[0][3][0][0]
        self.analysis_summary['Species'] = 'Mouse'
        self.analysis_summary['Age'] = 'p18'
        self.analysis_summary['Sex'] = 'U'
        self.analysis_summary['Weight'] = '28'
        self.analysis_summary['Temperature'] = '34.0'
        
        
        self.updateAnalysis(script_header=self.headerPrint)
        self.headerPrint = False  # reset after the first time
        if filename in self.allAnalysis.keys():
            self.allAnalysis[filename].append(self.analysis_summary)
        else:
            self.allAnalysis[filename] = [self.analysis_summary]  # make a list, keep it ordered

    def loadFileRequested(self, dh, analyze=True, bridge=None):
        """
        loadFileRequested is called by "file loader" when a file is requested.
            FileLoader is provided by the AnalysisModule class
            dh is the handle to the currently selected directory (or directories)

        This function loads all of the successive records from the specified protocol.
        Ancillary information from the protocol is stored in class variables.
        Extracts information about the commands, sometimes using a rather
        simplified set of assumptions. Much of the work for reading the data is
        performed in the GetClamps class in PatchEPhys.
        :param dh: the directory handle (or list of handles) representing the selected
        entitites from the FileLoader in the Analysis Module
        :modifies: plots, sequence, data arrays, data mode, etc.
        :return: True if successful; otherwise raises an exception
        """

        #dh = dh[0]  # just get the first one
        self.filename = dh.name()
        self.current_dirhandle = dh  # this is critical!
        self.loaded = dh
        self.analysis_summary = self.dataModel.cell_summary(dh)  # get other info as needed for the protocol
        self.Clamps = self.dataModel.GetClamps()
        ci = self.Clamps.getClampData(dh)
        self.Clamps.rgnrmp = [0., 0.005]
        if ci is None:
            return False
        if bridge is not None:
            self.bridgeCorrection = bridge
            #for i in range(self.Clamps.traces.shape[0]):
            print '******** Doing bridge correction: ', self.bridgeCorrection
            self.Clamps.traces = self.Clamps.traces - (self.bridgeCorrection * self.Clamps.cmd_wave)
        else:
            br = 20e6 # self.ctrl.IVCurve_bridge.value()*1e6
            # print 'br: ', br
            if br != 0.0:
                self.bridgeCorrection = br
                self.Clamps.traces = self.Clamps.traces - (self.bridgeCorrection * self.Clamps.cmd_wave)
            else:
                self.bridgeCorrection = None
        y = self.splitall(self.filename)
        shortname = os.path.join(y[-4], y[-3], y[-2])
        self.analysis_summary['Cell'] = shortname
        
        self.analysis_summary['Genotype'] = 'CBA/CaJ'
        self.analysis_summary['EndRec'] = 0
        self.analysis_summary['Block'] = 0
        (p, run) = os.path.split(self.filename)
        self.analysis_summary['FirstRec'] = int(run[-3:])
        # self.analysis_summary['Species'] = 'Mouse'
        # self.analysis_summary['Age'] = 'p18'
        # self.analysis_summary['Sex'] = 'U'
        # self.analysis_summary['Weight'] = '28'
        # self.analysis_summary['Temperature'] = '34.0'
        # now plot the data 
        self.plot_traces()
        self.updateAnalysis(script_header=self.headerPrint)
        self.headerPrint = False  # reset after the first time
        if shortname in self.allAnalysis.keys():
            self.allAnalysis[shortname].append(self.analysis_summary)
        else:
            self.allAnalysis[shortname] = [self.analysis_summary]  # make a list, keep it ordered
    
    def splitall(self, path):
        allparts = []
        while 1:
            parts = os.path.split(path)
            if parts[0] == path:  # sentinel for absolute paths
                allparts.insert(0, parts[0])
                break
            elif parts[1] == path: # sentinel for relative paths
                allparts.insert(0, parts[1])
                break
            else:
                path = parts[0]
                allparts.insert(0, parts[1])
        return allparts
    
        
    def plot_traces(self):
        """
        Plot the current data traces.
        :param multimode: try using "multiline plot routine" to speed up plots (no color though)
        :return: nothing
        """
        self.clearDecorators()
        self.make_map_symbols()
        ntr = self.Clamps.traces.shape[0]
        self.data_plot.setDownsampling(auto=False, mode='mean')
        self.data_plot.setClipToView(False)  # setting True deletes some points used for decoration of spikes by shape
        self.cmd_plot = self.plt2
        self.cmd_plot.setDownsampling(auto=False, mode='mean')
        self.cmd_plot.setClipToView(True)  # can leave this true since we do not put symbols on the plot
        self.data_plot.disableAutoRange()
        self.cmd_plot.disableAutoRange()
        cmdindxs = np.unique(self.Clamps.commandLevels)  # find the unique voltages
        colindxs = [int(np.where(cmdindxs == self.Clamps.commandLevels[i])[0]) for i in range(len(self.Clamps.commandLevels))]  # make a list to use

        for i in range(ntr):
            atrace = self.Clamps.traces[i]
            acmdwave = self.Clamps.cmd_wave[i]
            self.plt1.plot(x=self.Clamps.time_base, y=atrace, downSample=10, downSampleMethod='mean',
                                pen=pg.intColor(colindxs[i], len(cmdindxs), maxValue=255))
            self.plt2.plot(x=self.Clamps.time_base, y=acmdwave, downSample=10, downSampleMethod='mean',
                               pen=pg.intColor(colindxs[i], len(cmdindxs), maxValue=255))
        self.plt1.autoRange()
        self.plt2.autoRange()

    def updateAnalysis(self, printnow=True, script_header=False):
        self.SA = SpikeAnalysis.SpikeAnalysis()  # create instances of our analysis routines
        threshold = 0
        self.SA.setup(self.Clamps, threshold)
        self.SA.FIGrowth = 1
        self.RmTau = RmTauAnalysis.RmTauAnalysis()
        self.RmTau.setup(clamps=self.Clamps, spikes=self.SA, 
                dataplot=None, baseline = self.Clamps.rgnrmp)
        
        self.analyzeSpikes()
        rgnpk = [self.Clamps.tstart, self.Clamps.tstart+0.15]  # times in seconds
        self.RmTau.tau_membrane(peak_time=None, printWindow=False, whichTau=1,
                vrange=[-0.005, -0.020], region=rgnpk)
        self.update_ssAnalysis()
        self.update_pkAnalysis()
        self.analyzeSpikeShape()
        # print '-'*80
        # print 'Cell: %14s  ' % (self.currentBlock.datac.fname)
        # print self.currentBlock.summary()
        # print 'RMP: %8.3f   Taum: %8.3f  AR: %8.3f   HW: %8.3f' % (
        #     self.RmTau.rmp, self.analysis_summary['tau'],
        #     self.analysis_summary['AdaptRatio'], self.analysis_summary['AP1_HalfWidth'])
        FIKeys = self.SA.FIKeys # ['Ibreak', 'Irate', 'IPower']# ['Ibreak', 'Rate0', 'Ibreak1', 'Irate1', 'Irate2', 'Irate3']
        for i, k in enumerate(FIKeys):
            print '%-15s\t' % k,
        print
        for i, k in enumerate(FIKeys):
            print '%12.6f\t' % (self.analysis_summary[k]),
        print '\n'
        print '-'*80
        ltxt = ''
        # summary table header is written anew for each cell
        htxt = ''
        if script_header:
        #    htxt = '{:14s}\t{:12s}\t{:12s}\t'.format("Cell", "Genotype", "Protocol")
            for k in self.data_template.keys():
                cnv = '{:<%ds}' % (self.data_template[k][0])
                # print 'cnv: ', cnv
                htxt += (cnv + '\t').format(k)
            script_header = False
            htxt += '\n'
        for a in self.data_template.keys():
            if a in self.analysis_summary.keys():
                txt = self.analysis_summary[a]
                if a in ['Description', 'Notes']:
                    txt = txt.replace('\n', ' ').replace('\r', '')  # remove line breaks from output, replace \n with space
                #print a, data_template[a]
                ltxt += (self.data_template[a][1]).format(txt) + ' \t'
            else:
                ltxt += ('{:>%ds}' % (self.data_template[a][0]) + '\t').format('NaN')
        ltxt = ltxt.replace('\n', ' ').replace('\r', '')  # remove line breaks
        ltxt = htxt + ltxt
        if self.datafile is not None:
            outf = open(self.datafile, 'a+')
            outf.write(ltxt + '\n')
            outf.close()

        if printnow:
            print ltxt
        #print pprint.pformat(self.analysis_summary)
        return

    def analyzeSpikes(self):
        """
        analyzeSpikes: Using the threshold set in the control panel, count the
        number of spikes in the stimulation window (self.Clamps.tstart, self.Clamps.tend)
        Updates the spike plot(s).

        The following variables are set:
        self.SA.spikecount: a 1-D numpy array of spike counts, aligned with the
            current (command)
        self.adapt_ratio: the adaptation ratio of the spike train
        self.fsl: a numpy array of first spike latency for each command level
        self.fisi: a numpy array of first interspike intervals for each
            command level
        self.nospk: the indices of command levels where no spike was detected
        self.spk: the indices of command levels were at least one spike
            was detected
        """
        clearFlag = True
        self.analysis_summary['FI_Curve'] = None
        # print '***** analyzing Spikes'
        if self.Clamps.data_mode not in self.ic_modes or self.Clamps.time_base is None:
            print ('IVCurve::analyzeSpikes: Cannot count spikes, ' +
                   'and dataMode is ', self.Clamps.data_mode, 'and ICModes are: ', self.ic_modes, 'tx is: ', self.tx)
            self.SA.spikecount = []
            self.fiPlot.plot(x=[], y=[], clear=clearFlag, pen='w',
                             symbolSize=6, symbolPen='b',
                             symbolBrush=(0, 0, 255, 200), symbol='s')
            self.fslPlot.plot(x=[], y=[], pen='w', clear=clearFlag,
                              symbolSize=6, symbolPen='g',
                              symbolBrush=(0, 255, 0, 200), symbol='t')
            self.fslPlot.plot(x=[], y=[], pen='w', symbolSize=6,
                              symbolPen='y',
                              symbolBrush=(255, 255, 0, 200), symbol='s')
            return
        threshold = 0 # self.ctrl.IVCurve_SpikeThreshold.value() * 1e-3
        self.analysis_summary['SpikeThreshold'] = threshold # self.ctrl.IVCurve_SpikeThreshold.value()
        self.update_rmpAnalysis()
                
        self.SA.setup(self.Clamps, threshold)
        self.SA.FIGrowth = 1 # Choose fit type
        self.SA.analyzeSpikes()
        
        self.adapt_ratio = self.SA.adapt_ratio # np.mean(ar[iAR])  # only where we made the measurement
        self.analysis_summary['AdaptRatio'] = self.SA.adapt_ratio
#        self.ctrl.IVCurve_AR.setText(u'%7.3f' % self.SA.adapt_ratio)
        self.nospk = np.where(self.SA.spikecount == 0)
        self.spk = np.where(self.SA.spikecount > 0)[0]
        self.analysis_summary['FI_Curve'] = np.array([self.Clamps.values, self.SA.spikecount])
        self.spikes_counted = True
        self.update_SpikePlots()

    def _timeindex(self, t):
        return np.argmin(self.Clamps.time_base-t)

    def analyzeSpikeShape(self, printSpikeInfo=False):
        self.SA.analyzeSpikeShape()
        self.spikeShape = self.SA.spikeShape
        self.analysis_summary['spikes'] = self.SA.spikeShape  # save in the summary dictionary too       
        self.analysis_summary['iHold'] = np.mean(self.SA.iHold)
        self.analysis_summary['pulseDuration'] = self.Clamps.tend - self.Clamps.tstart
        # self.getClassifyingInfo()  # build analysis summary here as well.
        # copy all the analysis summary from the SA to here.
        for k in self.SA.analysis_summary.keys():
            self.analysis_summary[k] = self.SA.analysis_summary[k]
#        self.clearDecorators()
        self.spikeDecorator()
#        self.update_IVPlot()

    def spikeDecorator(self):
        """
        Put markers on the spikes to visually confirm the analysis of thresholds, etc.
        """
        # get colors
        cmdindxs = np.unique(self.Clamps.commandLevels)  # find the unique voltages
        colindxs = [int(np.where(cmdindxs == self.Clamps.commandLevels[i])[0]) for i in range(len(self.Clamps.commandLevels))]  # make a list to use
        alllats = []
        allpeakt = []
        allpeakv = []
        for i, trace in enumerate(self.spikeShape):
            aps = []
            tps = []
            paps = []
            ptps = []
            taps = []
            ttps = []
            hwv = []
            tups = []
            tdps = []

            for j, spk in enumerate(self.spikeShape[trace]):
                aps.append(self.spikeShape[trace][spk]['AP_beginV'])
                alllats.append(self.spikeShape[trace][spk]['AP_Latency'])
                tps.append(self.spikeShape[trace][spk]['AP_Latency'])
            if len(tps) == 1:
                tps.append(tps[0]+1e-12)
                aps.append(aps[0]+1e-10)
            u = self.data_plot.plot(np.array(tps), np.array(aps), pen=None, symbol='o', brush=pg.mkBrush('g'), symbolSize=4)
            dn = 'd%03d' % i
            self.dataMarkers[dn] = u
            for j, spk in enumerate(self.spikeShape[trace]):
                paps.append(self.spikeShape[trace][spk]['peak_V'])
                ptps.append(self.spikeShape[trace][spk]['peak_T'])
                allpeakt.append(self.spikeShape[trace][spk]['peak_T']+0.01)
                allpeakv.append(self.spikeShape[trace][spk]['peak_V'])
            # u = self.data_plot.plot(allpeakt, allpeakv, pen=None, symbol='o', brush=pg.mkBrush('r'), size=2)
            # self.dataMarkers.append(u)
            if len(ptps) == 1:
                ptps.append(ptps[0]+1e-10)
                paps.append(paps[0]+1e-10)

            u = self.data_plot.plot(ptps, paps, pen=None, symbol='t', brush=pg.mkBrush('w'), symbolSize=4)
            dn = 'p%03d' % i
            self.dataMarkers[dn] = u

            for j, spk in enumerate(self.spikeShape[trace]):
                taps.append(self.spikeShape[trace][spk]['trough_V'])
                ttps.append(self.spikeShape[trace][spk]['trough_T'])
            if len(ttps) == 1:
                ttps.append(ttps[0]+1e-10)
                taps.append(taps[0]+1e-10)
            u = self.data_plot.plot(ttps, taps, pen=None, symbol='+', brush=pg.mkBrush('r'), symbolSize=4)
            dn = 't%03d' % i
            self.dataMarkers[dn] = u
            
            for j, spk in enumerate(self.spikeShape[trace]):
                tups.append(self.spikeShape[trace][spk]['hw_up'])
                tdps.append(self.spikeShape[trace][spk]['hw_down'])
                hwv.append(self.spikeShape[trace][spk]['hw_v'])
            if len(tups) == 1:
                tups.append(tups[0]+1e-10)
                tdps.append(tdps[0]+1e-10)
                hwv.append(hwv[0])
            u =self.data_plot.plot(tups, hwv, pen=None, symbol='d', brush=pg.mkBrush('c'), symbolSize=4)
            dn = 'du%03d' % i
            self.dataMarkers[dn] = u
            
            u = self.data_plot.plot(tdps, hwv, pen=None, symbol='s', brush=pg.mkBrush('c'), symbolSize=4)
            dn = 'dd%03d' % i
            self.dataMarkers[dn] = u

    def clearDecorators(self):
        if len(self.dataMarkers) > 0:
            [self.dataMarkers[m].clear() for m in self.dataMarkers.keys()]
        self.dataMarkers = {}        

    def clearFits(self):
        if len(self.fitMarkers) > 0:
            [self.fitMarkers[m].clear() for m in self.fitMarkers.keys()]
        self.fitMarkers = {}  
        try:
            if len(self.RmTau.taum_fitted) > 0:
                [self.RmTau.taum_fitted[m].clear() for m in self.RmTau.taum_fitted.keys()]
            self.self.RmTau.taum_fitted = {}  
        except:
            pass

    def update_SpikePlots(self):
        """
        Draw the spike counts to the FI and FSL windows
        Note: x axis can be I, T, or  # spikes
        """
        if self.Clamps.data_mode in self.vc_modes:
            self.fiPlot.clear()  # no plots of spikes in VC
            self.fslPlot.clear()
            return
        (pen, filledbrush, emptybrush, symbol, n, clearFlag) = self.map_symbol()
        fitpen = pg.mkPen('r')
        mode = 1
        #mode = self.ctrl.IVCurve_RMPMode.currentIndex()  # get x axis mode
        self.spcmd = self.Clamps.commandLevels[self.spk]  # get command levels iwth spikes
        iscale = 1.0e12  # convert to pA
        yfslsc = 1.0  # convert to msec
        if mode == 0:  # plot with time as x axis
            xfi = self.Clamps.trace_StartTimes
            xfsl = self.Clamps.trace_StartTimes
            select = range(len(self.Clamps.trace_StartTimes))
            xlabel = 'T (s)'
        elif mode == 1:  # plot with current as x
            select = self.spk
            xfi = self.Clamps.commandLevels * iscale
            xfsl = self.spcmd * iscale
            xlabel = 'I (pA)'
        elif mode == 2:  # plot with spike counts as x
            xfi = self.SA.spikecount
            xfsl = self.SA.spikecount
            select = range(len(self.SA.spikecount))
            xlabel = 'Spikes (N)'
        else:
            return  # mode not in available list
        self.fiPlot.plot(x=xfi, y=self.SA.spikecount, clear=clearFlag,
                         symbolSize=6,
                         symbol=symbol, pen=pen,
                         symbolPen=pen, symbolBrush=filledbrush)
        # also fit the data and compute FI values
        # xdata must be current levels
        xfit = self.Clamps.commandLevels * iscale
        if np.max(self.SA.spikecount) > 0:
            print('update spike plots calling fitOne, with func set to: ', self.SA.FIGrowth)
            (fpar, xf, yf, names, error, f, func) = self.SA.fitOne(xfit, self.SA.spikecount, info='FI',
                    fixNonMonotonic=True, excludeNonMonotonic=False)
            f = self.fiPlot.plot(x=xf[0], y=yf[0], clear=False, pen=fitpen)
            fn = 'fi'
            self.fitMarkers[fn] = f
            self.analysis_summary['FIPlot_fit'] = [xf[0], yf[0]]
            print('func = %s' % func)
            self.analysis_summary['FIFitFunc'] = func
            
        fslmax = 0.
        if self.showFISI:
            self.fslPlot.plot(x=xfsl, y=self.SA.fsl[select] * yfslsc, clear=clearFlag,
                          symbolSize=6,
                          symbol=symbol, pen=pen,
                          symbolPen=pen, symbolBrush=filledbrush)
            self.fslPlot.plot(x=xfsl, y=self.SA.fisi[select] * yfslsc, symbolSize=6,
                          symbol=symbol, pen=pen,
                          symbolPen=pen, symbolBrush=emptybrush)
            if len(xfsl) > 0:
                self.fslPlot.setXRange(0.0, np.max(xfsl))
                self.fslPlot.setYRange(0., max(max(self.SA.fsl[select]), max(self.SA.fisi[select])))
            ylabel = 'Fsl/Fisi (ms)'
            xfsllabel = xlabel
            self.fslPlot.setTitle('FSL/FISI')
        else:
            maxspk = 0
            maxisi = 0.
            clear = clearFlag
            for i, k in enumerate(self.allisi.keys()):
                nspk = len(self.allisi[k])
                xisi = np.arange(nspk)
                self.fslPlot.plot(x=xisi, y=self.SA.allisi[k] * yfslsc, clear=clear,
                              symbolSize=6,
                              symbol=symbol, pen=pen,
                              symbolPen=pen, symbolBrush=filledbrush)
                clear = False
                maxspk = max(nspk, maxspk)
                maxisi = max(np.max(self.allisi[k]), maxisi)
            self.fslPlot.setXRange(0.0, maxspk)
            self.fslPlot.setYRange(0.0, maxisi)
            xfsllabel = 'Spike Number'
            ylabel = 'ISI (s)'
            self.fslPlot.setTitle('ISI vs. Spike Number')
        self.fiPlot.setLabel('bottom', xlabel)
        self.fslPlot.setLabel('bottom', xfsllabel)
        self.fslPlot.setLabel('left', ylabel)

    def update_ssAnalysis(self):
        """
        Compute the steady-state IV from the selected time window

        Parameters
        ----------
            None.
        
        Returns
        -------
            nothing.
        
        modifies:
            ivss, yleak, ivss_cmd, cmd.

        The IV curve is only valid when there are no spikes detected in
            the window. The values in the curve are taken as the mean of the
            current and the voltage in the time window, at each command step.
        We also compute the input resistance.
        For voltage clamp data, we can optionally remove the "leak" current.
        The resulting curve is plotted.
        """
        if self.Clamps.traces is None:
            return
        rgnss = [self.Clamps.tend -0.1, self.Clamps.tend]
        #self.regions['lrwin1']['region'].getRegion()
        r1 = rgnss[1]
        if rgnss[1] == rgnss[0]:
            print 'Steady-state regions have no width; using 100 msec. window for ss '
            r1 = rgnss[0] + 0.1
#        self.ctrl.IVCurve_ssTStart.setValue(rgnss[0] * 1.0e3)
#        self.ctrl.IVCurve_ssTStop.setValue(r1 * 1.0e3)
        data1 = self.Clamps.traces['Time': rgnss[0]:r1]
        if data1.shape[1] == 0 or data1.shape[0] == 1:
            return  # skip it
        self.ivss = []

        # check out whether there are spikes in the window that is selected
        threshold = self.threshold # ctrl.IVCurve_SpikeThreshold.value() * 1e-3
        ntr = len(self.Clamps.traces)
        if not self.spikes_counted:
            print 'updatess: spikes not counted yet? '
            self.analyzeSpikes()

        self.ivss = data1.mean(axis=1)  # all traces
        if self.sub_baseline:
            self.ivss = self.ivss - self.RmTau.ivbaseline

        if len(self.nospk) >= 1:
            # Steady-state IV where there are no spikes
            self.ivss = self.ivss[self.nospk]
            self.ivss_cmd = self.Clamps.commandLevels[self.nospk]
#            self.commandLevels = commands[self.nospk]
            # compute Rin from the SS IV:
            # this makes the assumption that:
            # successive trials are in order (as are commands)
            # commands are not repeated...
            if len(self.ivss_cmd) > 0 and len(self.ivss) > 0:
                self.r_in = np.max(np.diff
                                   (self.ivss) / np.diff(self.ivss_cmd))
                #self.ctrl.IVCurve_Rin.setText(u'%9.1f M\u03A9'
                #                              % (self.r_in * 1.0e-6))
                self.analysis_summary['Rin'] = self.r_in*1.0e-6
            else:
                print (u'No Valid Points')
                #self.ctrl.IVCurve_Rin.setText(u'No valid points')
        self.yleak = np.zeros(len(self.ivss))
        if self.sub_leak:
        # self.ctrl.IVCurve_subLeak.isChecked():
            if self.Clamps.data_mode in self.dataModel.ic_modes:
                sf = 1e-12
            elif self.Clamps.data_mode in self.dataModel.vc_modes:
                sf = 1e-3
            else:
                sf = 1.0
            (x, y) = Utility.clipdata(self.ivss, self.ivss_cmd,
                                      self.ctrl.IVCurve_LeakMin.value() * sf,
                                      self.ctrl.IVCurve_LeakMax.value() * sf)
            try:
                p = np.polyfit(x, y, 1)  # linear fit
                self.yleak = np.polyval(p, self.ivss_cmd)
                self.ivss = self.ivss - self.yleak
            except: 
                raise ValueError('IVCurve Leak subtraction: no valid points to correct')
        isort = np.argsort(self.ivss_cmd)
        self.ivss_cmd = self.ivss_cmd[isort]
        self.ivss = self.ivss[isort]
        self.analysis_summary['IV_Curve_ss'] = [self.ivss_cmd, self.ivss]



    @staticmethod
    def label_up(plot, xtext, ytext, title):
        """helper to label up the plot"""
        plot.setLabel('bottom', xtext)
        plot.setLabel('left', ytext)
        plot.setTitle(title)
        
    def update_Tau_membrane(self, peak_time=None, printWindow=True, whichTau=1, vrange=[-5., -20.]):
        """
        Compute time constant (single exponential) from the
        onset of the response
        using lrpk window, and only steps that produce a voltage change between 5 and 20 mV below rest
        or as specified
        """

        if len(self.Clamps.commandLevels) == 0:  # probably not ready yet to do the update.
            return
        if self.Clamps.data_mode not in self.ic_modes:  # only permit in IC
            return
        rgnpk = [self.Clamps.tstart, self.Clamps.tstart + 0.1] # list(self.regions['lrwin0']['region'].getRegion())
        self.RmTau.setup(clamps=self.Clamps, spikes=self.SA, 
                dataplot=self.data_plot)
        self.RmTau.tau_membrane(peak_time=None, printWindow=False, whichTau=1,
                vrange=[-0.005, -0.020], region=rgnpk)
        
        #self.ctrl.IVCurve_Tau.setText(u'%18.1f ms' % (self.RmTau.taum_taum * 1.e3))
        self.analysis_summary['tau'] = self.RmTau.taum_taum*1.e3
        tautext = 'Mean Tau: %8.1f'
        if printWindow:
            print tautext % (self.RmTau.taum_taum * 1e3)
        self.show_taum_plot()

    def show_taum_plot(self):
        Fits = Fitting.Fitting()
        fitPars = self.RmTau.taum_pars
        xFit = np.zeros((len(self.RmTau.taum_pars), 500))
        for i in range(len(self.RmTau.taum_pars)):
          xFit[i,:] = np.arange(0, self.RmTau.taum_win[1]-self.RmTau.taum_win[0],
                   (self.RmTau.taum_win[1]-self.RmTau.taum_win[0])/500.)
        yFit = np.zeros((len(fitPars), xFit.shape[1]))
        fitfunc = Fits.fitfuncmap[self.RmTau.taum_func]

        for k, whichdata in enumerate(self.RmTau.taum_whichdata):
            yFit[k] = fitfunc[0](fitPars[k], xFit[k], C=None)  # +self.ivbaseline[whichdata]
            n = 'r%03d' % k
            self.RmTau.taum_fitted[n] = self.data_plot.plot(xFit[k]+self.RmTau.taum_win[0], yFit[k],
                pen=pg.mkPen('r', width=2, style=QtCore.Qt.DashLine))
#        print 'whichdata: ', self.RmTau.taum_whichdata
        
    def update_pkAnalysis(self, clear=False, pw=False):
        """
        Compute the peak IV (minimum) from the selected window
        mode can be 'min', 'max', or 'abs'

        Parameters
        ----------
        clear : Boolean, False
        pw : Boolean, False
            pw is passed to update_taumembrane to control printing.
        """
        if self.Clamps.traces is None:
            return
        mode = 'Min' # self.ctrl.IVCurve_PeakMode.currentText()
        rgnpk = [self.Clamps.tstart, self.Clamps.tstart+0.1] # self.regions['lrwin0']['region'].getRegion()
        #self.ctrl.IVCurve_pkTStart.setValue(rgnpk[0] * 1.0e3)
        #self.ctrl.IVCurve_pkTStop.setValue(rgnpk[1] * 1.0e3)
        data2 = self.Clamps.traces['Time': rgnpk[0]:rgnpk[1]]
        if data2.shape[1] == 0:
            return  # skip it - window missed the data
        # check out whether there are spikes in the window that is selected
        # but only in current clamp
        nospk = []
        peak_pos = None
        if self.Clamps.data_mode in self.ic_modes:
            threshold = self.threshold # self.ctrl.IVCurve_SpikeThreshold.value() * 1e-3
            ntr = len(self.Clamps.traces)
            if not self.spikes_counted:
                print 'update_pkAnalysis: spikes not counted'
                self.analyzeSpikes()
            spikecount = np.zeros(ntr)

        if mode == 'Min':
            self.ivpk = data2.min(axis=1)
            peak_pos = np.argmin(data2, axis=1)
        elif mode == 'Max':
            self.ivpk = data2.max(axis=1)
            peak_pos = np.argmax(data2, axis=1)
        elif mode == 'Abs':  # find largest regardless of the sign ('minormax')
            x1 = data2.min(axis=1)
            peak_pos1 = np.argmin(data2, axis=1)
            x2 = data2.max(axis=1)
            peak_pos2 = np.argmax(data2, axis=1)
            self.ivpk = np.zeros(data2.shape[0])
            for i in range(data2.shape[0]):
                if -x1[i] > x2[i]:
                    self.ivpk[i] = x1[i]
                    peak_pos = peak_pos1
                else:
                    self.ivpk[i] = x2[i]
                    peak_pos = peak_pos2
                    # self.ivpk = np.array([np.max(x1[i], x2[i]) for i in range(data2.shape[0]])
                    #self.ivpk = np.maximum(np.fabs(data2.min(axis=1)), data2.max(axis=1))
        if self.sub_baseline: # self.ctrl.IVCurve_SubBaseline.isChecked():
            self.ivpk = self.ivpk - self.RmTau.ivbaseline
        if len(self.nospk) >= 1:
            # Peak (min, max or absmax voltage) IV where there are no spikes
            self.ivpk = self.ivpk[self.nospk]
            self.ivpk_cmd = self.Clamps.commandLevels[self.nospk]
        else:
            self.ivpk_cmd = self.Clamps.commandLevels
        self.ivpk = self.ivpk.view(np.ndarray)
        if self.sub_leak: # self.ctrl.IVCurve_subLeak.isChecked():
            self.ivpk = self.ivpk - self.yleak
        # now sort data in ascending command levels
        isort = np.argsort(self.ivpk_cmd)
        self.ivpk_cmd = self.ivpk_cmd[isort]
        self.ivpk = self.ivpk[isort]
        self.analysis_summary['IV_Curve_pk'] = [self.ivpk_cmd, self.ivpk]
        self.update_IVPlot()
        peak_time = self.Clamps.time_base[peak_pos]
        self.update_Tau_membrane(peak_time=peak_time, printWindow=pw)

    def update_rmpAnalysis(self, **kwargs):
        """
        Compute the RMP over time/commands from the selected window
        """
        try:
            if self.Clamps.traces is None:
                return
        except:
            return
#        if self.RmTau.Clamps is None:
        self.RmTau.setup(clamps=self.Clamps, spikes=self.SA,
            dataplot=self.data_plot)
        rgnrmp = [0., 0.005] # self.regions['lrrmp']['region'].getRegion()
        self.RmTau.rmp_analysis(rgnrmp)
        #self.ctrl.IVCurve_rmpTStart.setValue(rgnrmp[0] * 1.0e3)
        #self.ctrl.IVCurve_rmpTStop.setValue(rgnrmp[1] * 1.0e3)
        #self.ctrl.IVCurve_vrmp.setText('%8.2f' % self.RmTau.rmp)
        #self.update_RMPPlot()
        self.analysis_summary['RMP'] = self.RmTau.rmp


    def update_IVPlot(self):
        """
        Draw the peak and steady-sate IV to the I-V window
        Note: x axis is always I or V, y axis V or I
        """
        if not self.keep_analysis:
#            if self.ctrl.IVCurve_KeepAnalysis.isChecked() is False:
            self.IV_plot.clear()
        (pen, filledbrush, emptybrush, symbol, n, clearFlag) = \
            self.map_symbol()
        if self.Clamps.data_mode in self.ic_modes:
            if len(self.ivss) > 0 and self.lrss_show:
                    #self.ctrl.IVCurve_showHide_lrss.isChecked()):
                self.IV_plot.plot(self.ivss_cmd * 1e12, self.ivss * 1e3,
                                  symbol=symbol, pen=pen,
                                  symbolSize=6, symbolPen=pen,
                                  symbolBrush=filledbrush)
            if len(self.ivpk) > 0 and self.lrpk_show: # self.ctrl.IVCurve_showHide_lrpk.isChecked()):
                self.IV_plot.plot(self.ivpk_cmd * 1e12, self.ivpk * 1e3,
                                  symbol=symbol, pen=pen,
                                  symbolSize=6, symbolPen=pen,
                                  symbolBrush=emptybrush)
            self.label_up(self.IV_plot, 'I (pA)', 'V (mV)', 'I-V (CC)')
        if self.Clamps.data_mode in self.vc_modes:
            if len(self.ivss) > 0 and self.lrss_show: # self.ctrl.IVCurve_showHide_lrss.isChecked()):
                self.IV_plot.plot(self.ivss_cmd * 1e3, self.ivss * 1e9,
                                  symbol=symbol, pen=pen,
                                  symbolSize=6, symbolPen=pen,
                                  symbolBrush=filledbrush)
            if len(self.ivpk) > 0 and self.lrpk_show: # self.ctrl.IVCurve_showHide_lrpk.isChecked()):
                self.IV_plot.plot(self.ivpk_cmd * 1e3, self.ivpk * 1e9,
                                  symbol=symbol, pen=pen,
                                  symbolSize=6, symbolPen=pen,
                                  symbolBrush=emptybrush)
            self.label_up(self.IV_plot, 'V (mV)', 'I (nA)', 'I-V (VC)')        
        
    def make_map_symbols(self):
        """
        Given the current state of things, (keeping the analysis, when
        superimposing multiple results, for example),
        sets self.currentSymDict with a dict of pen, fill color, empty color, a symbol from
        our lists, and a clearflag. Used to overplot different data.
        """
        n = self.keep_analysis_count
        pen = self.color_list.next()
        filledbrush = pen
        emptybrush = None
        symbol = self.symbol_list.next()
        if n == 0:
            clearFlag = True
        else:
            clearFlag = False
        self.currentSymDict = {'pen': pen, 'filledbrush': filledbrush,
                               'emptybrush': emptybrush, 'symbol': symbol,
                               'n': n, 'clearFlag': clearFlag}

    def map_symbol(self):
        cd = self.currentSymDict
        if cd['filledbrush'] == 'w':
            cd['filledbrush'] = pg.mkBrush((128, 128, 128))
        if cd['pen'] == 'w':
            cd['pen'] = pg.mkPen((128, 128, 128))
        self.lastSymbol = (cd['pen'], cd['filledbrush'],
                           cd['emptybrush'], cd['symbol'],
                           cd['n'], cd['clearFlag'])
        return self.lastSymbol


class EPSCAnalyzer(DatacBrowser):
    """Browser with additional features for analyzing stim-evoked EPSCs.
    """
    def __init__(self, path):
        DatacBrowser.__init__(self, path)
        self.lr1 = pg.LinearRegionItem([0, 0.01])
        self.lr2 = pg.LinearRegionItem([0, 0.01])
        self.plt1.addItem(self.lr1)
        self.plt3.addItem(self.lr2)
        
        self.split3 = QtGui.QSplitter(QtCore.Qt.Vertical, self)
        self.addWidget(self.split3)
        
        self.ioPlot = pg.PlotWidget()
        self.ioData = self.ioPlot.plot(symbol='o')
        self.split3.addWidget(self.ioPlot)
        self.lr3 = pg.LinearRegionItem([0, 1])
        self.ioPlot.addItem(self.lr3)
        
        self.histogram = pg.PlotWidget()
        self.histData = self.histogram.plot(fillLevel=0, fillBrush=(200, 0, 0, 255))
        self.split3.addWidget(self.histogram)
        
        self.lr1.sigRegionChangeFinished.connect(self.updateAnalysis)
        self.lr2.sigRegionChangeFinished.connect(self.updateAnalysis)
        self.lr3.sigRegionChangeFinished.connect(self.updateAnalysis)

    def blockSelected(self, item):
        DatacBrowser.blockSelected(self, item)
        
        if not isinstance(self.currentBlock, Block):
            return
        
        rate, recs = self.currentBlock.data()
        stimInd = recs[0][2].argmax()
        stimTime = float(stimInd) / rate
        self.lr1.sigRegionChangeFinished.disconnect(self.updateAnalysis)
        self.lr2.sigRegionChangeFinished.disconnect(self.updateAnalysis)
        
        try:
            self.lr1.setRegion([stimTime + 0.5e-3, stimTime + 5e-3])
            self.lr2.setRegion([stimTime, stimTime+20e-6])
            
        finally:
            self.lr1.sigRegionChangeFinished.connect(self.updateAnalysis)
            self.lr2.sigRegionChangeFinished.connect(self.updateAnalysis)
            
        self.updateAnalysis()
        
    def updateAnalysis(self):
        return
        rate, recs = self.currentBlock.data()
        r1 = [int(x*rate) for x in self.lr1.getRegion()]
        r2 = [int(x*rate) for x in self.lr2.getRegion()]
        
        epsc = np.empty(len(recs))
        stim = np.empty(len(recs))
        for i,rec in enumerate(recs):
            d1 = rec[0][r1[0]:r1[1]]
            d2 = rec[2][r2[0]:r2[1]]
            base = np.median(rec[0][:100])
            epsc[i] = d1.min() - base
            stim[i] = d2.mean()
        print 'r1, r2', r1, r2
            
        self.ioData.setData(stim, epsc)
        
        vals, bins = np.histogram(epsc, bins=40)
        self.histData.setData(bins, vals, stepMode=True)
        
        rgn = self.lr3.getRegion()
        print 'rgn0, 1: ', rgn
        mask = (stim > rgn[0]) & (stim < rgn[1])
        epscrgn = epsc[mask]
        #print "start:  %g   end:  %g   mean: %g   stdev: %g" % (rgn[0], rgn[1], epscrgn.mean(), epscrgn.std())
        if self.loaded is not None:
            name = os.path.splitext(os.path.split(self.loaded.fname)[1])[0]
        else:
            name = 'unknown'
        # print ('%s ' + "%g "*5) % (name, rgn[0], rgn[1], epscrgn.mean(), epscrgn.std(), epscrgn.std() / epscrgn.mean())
        


if __name__ == '__main__':
    import sys
    import commands
    app = pg.mkQApp()
    path = ''
    datasummary = '../DataSummaries'
    if len(sys.argv) > 1:
        dataset = sys.argv[1]
    if path == '':
        cname = commands.getoutput('scutil --get ComputerName')
        # if cname in ['Lytle']:
        #     paths = {'Rao': '/Users/pbmanis/Documents/Lab/Manuscripts and Abstracts in Progress/Rao-STDPModulation/' +
        #             'Rao-OriginalAnalysis(2011)/Data and Analysis/POB computer data/A1 STDP/',
        #              'Kratz': '/Volumes/Pegasus/ManisLab_Data3/Rao_Deepti/Additional-Kratz/FI_Experiments/',
        #              'ODonohue': '/Volumes/Pegasus/ManisLab_Data3/Rao_Deepti/Additional-Heather/'}
        # elif cname in ['Tampalpais']:
        #     paths = {'Rao': '/Users/pbmanis/Documents/Lab/Manuscripts-inProgress/Rao-STDPModulation/Data and Analysis/POB computer data/A1 STDP/',
        #              'Kratz': '/Volumes/Backup2B/ManisLab_Data3/Rao_Deepti/Additional_Kratz/FI_Experiments/',
        #              'ODonohue': '/Volumes/Backup2B/ManisLab_Data3/Rao_Deepti/Additional-Heather/'}

    
    # allow datasets to be concatenated (or not)
    # build dictionary containing dictionaries from above
    # Dict structure:
    # each entry has a key by drug condition.
    # the data is a list
    # the first element of the list is a list of tuples
    #   Each elemnt of the list is a tuple with the path and list if infiles.
    #   if the list has multiple elements, the data are concatenated with the other elemtns of the list
    # the second element of the dict's list is the name of the output file
    
    # infiles_dict = OrderedDict([('carb', [[(paths['Rao'], infiles_carb)], 'Rao_IV_summary_carb.p']),  # test against carb veh if no difference combine
    #                       ('pir75', [[(paths['Rao'], infiles_pir75)], 'Rao_IV_summary_pir75.p', '75 nM']),
    #                       ('pir10', [[(paths['Rao'], infiles_pir10)], 'Rao_IV_summary_pir10.p', '10uM']),
    #                       ('carbvu', [[(paths['Kratz'], infiles_mbk_vu), (paths['ODonohue'], infiles_hao_vu)], 'MBK_IV_summary_carbvu.p']),
    #                       ('carbveh', [[(paths['Kratz'], infiles_mbk_veh)], 'MBK_IV_summary_carbveh.p']),
    #                       ('oxo', [[(paths['Kratz'], infiles_mbk_oxo), (paths['ODonohue'], infiles_hao_oxo)], 'MBKHAO_IV_summary_oxo.p']),
    #                       ])
    fn = '/Users/pbmanis/Documents/data/matdatac/12sep07f.mat'

    lf = DatacFile(os.path.join(fn))

    blocksel = 3
    bname = 'db_%d' % blocksel
    dblk = None
    for i, item in enumerate(lf.items):
        if item.type in ['HFILE', 'NOTE']:
            continue
        if item.dfile['Block'] == blocksel:
            dblk = item
            iblk = i

    if dblk is None:
        print('Block %d not found' % blocksel)
        exit(1)

    print dblk.type
    print 'sr: ', dblk.data()[0]
    print len(dblk.data())
    print len(dblk.data()[1])  # 21 records
    print len(dblk.data()[1][0]) # 3 traces
#    print     lf.items[1].data()[1][0][0]  # voltage
#    print     lf.items[1].dfile  # dict of all dfile arguments
#    print     lf.items[1].stims # sfile and stimulus arrays (at low rate), including method code
    

    P = PH.regular_grid(2, 1, order='columns', figsize=(10, 5), showgrid=False,
                verticalspacing=0.08, horizontalspacing=0.08,
                margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.03, 'bottommargin': 0.1},
                labelposition=(-0.12, 0.95))

    sr = dblk.dfile['Actual_Rate']  # sample rate
    df = dblk.dfile['frame']  # data frame
    
    ntr = dblk.data()[1]
    npts = dblk.dfile['Points']['v']
    tb = np.linspace(0., npts/sr, npts)*1000.
    vax = P.axdict['A']
    iax = P.axdict['B']
    for itr, tr in enumerate(dblk.data()[1]):
        vax.plot(tb, tr[0], linewidth=0.75)  # first index, record, second index is current, voltage (in order)
        iax.plot(tb, tr[1], linewidth=0.75)
    PH.calbar(iax, calbar=[0., 50., 50., 50.], unitNames={'x': 'ms', 'y': 'pA'})
    PH.calbar(vax, calbar=[0., -40., 50., 20.], unitNames={'x': 'ms', 'y': 'mV'})
    #PH.reference_line(vax, refline=-65.)
    vax.set_ylim([-100., 50.])
    iax.set_ylim([-300., 300.])
    mpl.show()

    # else:
    #     if dataset in infiles_dict.keys():
    #         info = infiles_dict[dataset] # get the info for that dataset
    #         outdict = {}
    #         for i, f in enumerate(info[0]): # get the file information for that dataset
    #             thispath = info[0][i][0]
    #             thefiles = info[0][i][1]
    #             win = IVAnalyzer(path=thispath, infiles=thefiles, outfile=None) # 'Rao_IV_summary.txt')
    #             if thefiles is not None:
    #                 print win.allAnalysis.keys()
    #             if i == 0:
    #                 outdict = win.allAnalysis  # just copy the first result
    #             else:  # concatenate results in dictionaries. Because the keys are filenames, they should be unique
    #                 for k in win.allAnalysis.keys():
    #                     outdict[k] = win.allAnalysis[k]
    #         print outdict.keys()
    #         p = open(os.path.join(datasummary, info[1]), 'w')
    #         pickle.dump(outdict, p)
    #         p.close()
    #     else:
    #         raise ValueError('Unrecognized data set: ', dataset)
    win.show()
    if sys.flags.interactive == 0:
        app.exec_()
