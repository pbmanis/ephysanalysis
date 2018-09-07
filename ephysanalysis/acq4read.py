from __future__ import print_function
#!/usr/bin/python

"""
Class to read acq4 data blocks in simple manner, as a standalone program.
Does not require acq4 link; bypasses DataManager and PatchEPhys

Requires pyqtgraph to read the .ma files and the .index file

"""
import os
import re
#from pyqtgraph import metaarray
import matplotlib
matplotlib.use('Agg')
import ephysanalysis.metaarray as EM
from pyqtgraph import configfile
import numpy as np
import datetime
import matplotlib.pyplot as mpl
import pprint
import textwrap as WR
import collections
import tifffile as tf

pp = pprint.PrettyPrinter(indent=4)



class Acq4Read():
    """
    Provides methods to read an acq4 protocol directory
    including data and .index files
    """
    
    def __init__(self, pathtoprotocol=None, dataname=None):
        """
        Parameters
        ----------
        pathtoprotocol str (default: None)
            Path to the protocol directoryto set for this instance of the reader
        
        dataname: str (default: None)
            Name of the data file to read (for example, 'MultiClamp1.ma')
        
        Returns
        -------
        Nothing
        """
        
        self.protocol = None
        if pathtoprotocol is not None:
            self.setProtocol(pathtoprotocol)
        if dataname is None:
            self.dataname = 'MultiClamp1.ma'
        else:
            self.dataname = dataname
        self.clampInfo = {}
        self.lb = '\n'
        # establish known clamp devices:
        maxclamps = 4
        clamps = []
        self.clampdevices = []
        for nc in range(maxclamps):  # create a list of valid clamps and multiclamps
            cno = nc + 1 # clamps numbered from 1 typically, not 0
            cname = 'Clamp%d'%cno
            mcname = 'MultiClamp%d'%cno
            clamps.extend([(cname, 'Pulse_amplitude'), (mcname, 'Pulse_amplitude')])
            self.clampdevices.extend([cname, mcname])
        aon = 'AxoPatch200'
        apn = 'AxoProbe' 
        clamps.extend([(aon, 'Pulse_amplitude'), (apn, 'Pulse_amplitude')])
        self.clampdevices.extend([aon, apn])
        self.clamps = clamps
        self.tstamp = re.compile('\s*(__timestamp__: )([\d.\d]*)')
        
    def setProtocol(self, pathtoprotocol):
        """
        Parameters
        ----------
        pathtoprotocol str (default: None)
            Path to the protocol directory to set for this instance of the reader
        """
        self.protocol = pathtoprotocol
    
    def setDataName(self, dataname):
        """
        Set the type (name) of the data metaarray name that will be read
        """
        self.dataname = dataname

    def subDirs(self, p):
        """
        return a list of the subdirectories just below this path
        
        Parameters
        ----------
        p : str (no default)
            path to investigate
        """
        dirs = filter(os.path.isdir, [os.path.join(p, d) for d in os.listdir(p)])
        dirs = sorted(dirs)  # make sure these are in proper order... 
        return dirs

    def listSequenceParams(self, dh):
        """Given a directory handle for a protocol sequence, return the dict of sequence parameters"""
        try:
            return dh.info()['sequenceParams']
        except KeyError:
            if len(dh.info()) == 0:
                print( '****************** Error: Missing .index file? (fails to detect protocol sequence)')
                raise Exception("Directory '%s' does not appear to be a protocol sequence." % dh.name())

    def getIndex(self, currdir='', lineend='\n'):
        self.lb = lineend  # set line break character
        self._readIndex(currdir=currdir)
        if self._index is not None:
            return self._index['.']
        else:
            return None

    def _readIndex(self, currdir=''):
        self._index = None
        indexFile = os.path.join(self.protocol, currdir, '.index')
#        print self.protocol, currdir, indexFile
        if not os.path.isfile(indexFile):
           # print("Directory '%s' is not managed!" % (self.dataname))
            return self._index
        self._index = configfile.readConfigFile(indexFile)
        return self._index

    def readDirIndex(self, currdir=''):
        self._dirindex = None
        indexFile = os.path.join(currdir, '.index')
       # print (indexFile)
        if not os.path.isfile(indexFile):
            print("Directory '%s' is not managed!" % (currdir))
            return self._dirindex
        try:
            self._dirindex = configfile.readConfigFile(indexFile)
        except:
            print('Failed to read index file for %s' % currdir)
            print('Probably bad formatting or broken .index file')
            return self._dirindex
        return self._dirindex

    def _parse_timestamp(self, lstr):
        tstamp = None
        ts = self.tstamp.match(lstr)
        if ts is not None:
            fts = float(ts.group(2))
            tstamp = datetime.datetime.fromtimestamp(fts).strftime('%Y-%m-%d  %H:%M:%S %z')
        return tstamp

    def convert_timestamp(self, fts):
        tstamp = datetime.datetime.fromtimestamp(fts).strftime('%Y-%m-%d  %H:%M:%S %z')
        return tstamp
       
    def _parse_index(self, index):
        """
        Recursive version
        """
        self.indent += 1
        if isinstance(index, list):
            for i in range(len(index)):
                index[i] = self._parse_index(index[i])
                if isinstance(index[i], list):
                    self.textline += ('{0:s}  list, len={1:d}{2:s}'.format(' '*self.indent*4,  len(index[i]), self.lb))
                else:
                    if not isinstance(index[i], tuple):
                        self.textline += ('{0:s}  {1:d}{2:s}',format(' '*self.indent*4, index[i], self.lb))
        
        elif isinstance(index, tuple):
            self.textline += ('{0:s} Device, Sequence : {1:s}, {2:s}{3:s}'.format(' '*self.indent*4, str(index[0]), str(index[1]),
                self.lb))
 
        elif isinstance(index, dict):
            for k in index.keys():
                if k.endswith('.ma') or k.endswith('.tif'):
                    continue
                if k in ['splitter']:
                    continue

                index[k] = self._parse_index(index[k])
                if isinstance(index[k], list) or isinstance(index[k], np.ndarray):
                    self.textline += ('{0:s} {1:3d} : list/array, len= {2:4d}{3:s}'.format(' '*self.indent*4, k, len(index[k]),
                        self.lb))
                elif k not in ['__timestamp__', '.']:
                    indents = ' '*(self.indent*4)
                    indents2 = ' '*(self.indent*4)
                    # do a textwrap on ths string
                    if k in ['description', 'notes']:
                        hdr = ('{0:s} {1:>20s} : '.format(indents, k))
                      #  self.textline += hdr
                        wrapper = WR.TextWrapper(initial_indent='', subsequent_indent=len(hdr)*' ', width=100)
                        for t in wrapper.wrap(hdr + str(index[k])):
                            self.textline += t+self.lb
                    else:
                        if not isinstance(index[k], collections.OrderedDict):
                            self.textline += ('{0:s} {1:>20s} : {2:<s}{3:s}'.format(indents, k, str(index[k]), self.lb))
                        else:
                            break
                elif k in ['__timestamp__']:
                    tstamp = self.convert_timestamp(index[k])
                    if tstamp is not None:
                        self.textline += ('{0:s} {1:>20s} : {2:s}{3:s}'.format(' '*self.indent*4, 'timestamp', tstamp, self.lb))
        
        elif isinstance(index, bytes):  # change all bytestrings to string and remove internal quotes
            index = index.decode('utf-8').replace("\'", '')
            self.textline += ('{0:s}  b: {1:d}{2:s}'.format(' '*self.indent*4, inde, self.lb))
        self.indent -= 1
        return index
        
    def printIndex(self, index):
        """
        Generate a nice printout of the index, about as far down as we can go
        """
        self.indent = 0
        self.textline = ''
        t = self._parse_index(index)
        print('Index: \n', t)
        return

    def getIndex_text(self, index):
        """
        Generate a nice printout of the index, about as far down as we can go
        """
        self.indent = 0
        self.textline = ''
        t = self._parse_index(index)
        return self.textline
        
        # for k in index['.'].keys():
        #     print( '  ', k, ':  ', index['.'][k])
        #     if isinstance(index['.'][k], dict):
        #         for k2 in index['.'][k].keys():
        #             print ('    ', k, ' ', k2, '::  ', index['.'][k][k2])
        #             if isinstance(index['.'][k][k2], dict):
        #                 for k3 in index['.'][k][k2]:
        #                     print ('    ', k, ' ', k2, ' ', k3, ':::  ', index['.'][k][k2][k3])
        #                     if isinstance(index['.'][k][k2][k3], dict):
        #                         for k4 in index['.'][k][k2][k3]:
        #                             print( '    [', k, '][', k2, '][', k3, '][', k4, '] ::::  ', index['.'][k][k2][k3][k4])

    def file_cell_protocol(self, filename):
        """
        file_cell_protocol breaks the current filename down and returns a
        tuple: (date, cell, protocol)
        last argument returned is the rest of the path...
        """
        (p0, proto) = os.path.split(filename)
        (p1, cell) = os.path.split(p0)
        (p2, sliceid) = os.path.split(p1)
        (p3, date) = os.path.split(p2)
        return (date, sliceid, cell, proto, p3)

    def getClampDevices(self, currdir='', verbose=False):
        """
        Search for a known clamp device in the list of devices 
        used in the current protocol directory...
        
        Return
        ------
        list of valid clamp devices found (there may be more than one)
            List will be empty if no recognized device is found.
        """
        info = self.getIndex(currdir=currdir)
        if verbose:
            print('\ngetClampDevices info: ', info['devices'])
        devs = []
        if info is not None and 'devices' in info.keys():
            devices = info['devices']
            for d in devices:
                if d in self.clampdevices:
                    devs.append(d)
        return devs

    def getDataInfo(self, fn):
        """
        Get the index info for a record, without reading the trace data
        """
        info = None
        if (os.path.isfile(fn)):
            try:
                tr = EM.MetaArray(file=fn, readAllData=False)
            except:
                return info
            info = tr[0].infoCopy()
#            print ('info: ', info)
            self.parseClampInfo(info)
        return(info)

    def parseClampInfo(self, info):
        """
        Get important information from the info[1] directory that we can use
        to determine the acquisition type
        """
        self.mode = info[1]['ClampState']['mode']
        self.units = [info[1]['ClampState']['primaryUnits'], info[1]['ClampState']['secondaryUnits']]
        self.samp_rate = info[1]['DAQ']['primary']['rate']
        if self.mode in ['IC', 'I=0']:
            self.tracepos = 1
            self.cmdpos = 0
        elif self.mode in ['VC']:
            self.tracepos = 1
            self.cmdpos = 0
        else:
            raise ValueError('Unable to determine how to map channels')

    def parseClampWCCompSettings(self, info):
        """
        Given the .index file for this protocol dir, try to parse the 
        clamp state and compensation
        """
        d = {}
        if 'ClampState' in info[1].keys() and 'ClampParams' in info[1]['ClampState'].keys():
            par = info[1]['ClampState']['ClampParams']
            d['WCCompValid'] = True
            d['WCEnabled'] = par['WholeCellCompEnable']
            d['WCResistance'] = par['WholeCellCompResist']
            d['WCCellCap'] = par['WholeCellCompCap']
            d['CompEnabled'] = par['RsCompEnable']
            d['CompCorrection'] = par['RsCompCorrection']
            d['CompBW'] = par['RsCompBandwidth']
            return d
        else:
            return {'WCCompValid': False, 'WCEnable': 0, 'WCResistance': 0., 'WholeCellCap': 0.,
                    'CompEnable': 0, 'CompCorrection': 0., 'CompBW': 50000. }

    def parseClampCCCompSettings(self, info):
        d = {}
        if 'ClampState' in info[1].keys() and 'ClampParams' in info[1]['ClampState'].keys():
            par = info[1]['ClampState']['ClampParams']
            d['CCCompValid'] = True
            d['CCBridgeEnable'] = par['BridgeBalEnable']
            d['CCBridgeResistance'] = par['BridgeBalResist']
            d['CCNeutralizationEnable'] = par['NeutralizationEnable']
            d['CCNeutralizationCap'] = par['NeutralizationCap']
            d['CCLPF'] = par['PrimarySignalLPF']
            d['CCPipetteOffset'] = par['PipetteOffset']
            return d
        else:
            return {'CCCompValid': False, 'CCBridgeEnable': 0, 'CCBridgeResistance': 0., 'CCNeutralizationEnable': 0.,
                    'CCNeutralizationCap': 0, 'CCPipetteOffset': 0., 'CCLPF': 10000. }
    
    def parseClampHoldingLevel(self, info):
        """
        Given the .index file for a protocol dir, try to get
        the holding level from the clamp state
        """
        try:
            return info[1]['ClampState']['holding']
        except:
            return 0.

    def getData(self, pos=1, check=False):
        """
        Get the data for the current protocol
        if check is True, we just check that the requested file exists and return
        True if it does and false if it does not
        """ 
        # non threaded
        dirs = self.subDirs(self.protocol)
        index = self._readIndex()
        self.clampInfo['dirs'] = dirs
        self.clampInfo['missingData'] = []
        self.traces = []
        self.data_array = []
        self.commandLevels = []
        self.cmd_wave = []
        self.time_base = []
        self.values = []
        self.trace_StartTimes = np.zeros(0)
        self.sample_rate = []
        info = self.getIndex(self.protocol)
        holdcheck = info['devices']['MultiClamp1']['holdingCheck']
        holdvalue = info['devices']['MultiClamp1']['holdingSpin']
        if holdcheck:
            self.holding = holdvalue
        else:
            self.holding = 0.

        trx = []
        cmd = []
        sequence_values = None
        if 'sequenceParams' in index['.'].keys():
            self.sequence =  index['.']['sequenceParams']
        else:
            self.sequence = []
        # building command voltages or currents - get amplitudes to clamp

        reps = ('protocol', 'repetitions')
        foundclamp = False
        for clamp in self.clamps:
            if clamp in self.sequence:
                foundclamp = True
                self.clampValues = self.sequence[clamp]
                self.nclamp = len(self.clampValues)
                if sequence_values is not None:
                    sequence_values = [x for x in self.clampValues for y in sequence_values]
                else:
                    sequence_values = [x for x in self.clampValues]
        self.mode = None
        for i, d in enumerate(dirs):
            fn = os.path.join(d, self.dataname)
            if not os.path.isfile(fn):
                print(' acq4read.getData: File not found: ', fn)
                if check:
                    return False
                else:
                    continue
            if check:
                return True
            tr = EM.MetaArray(file=fn)
            info = tr[0].infoCopy()
            self.parseClampInfo(info)
            self.WCComp = self.parseClampWCCompSettings(info)
            self.CCComp = self.parseClampCCCompSettings(info)
            
            # if i == 0:
            #     pp.pprint(info)
            cmd = self.getClampCommand(tr)
            self.traces.append(tr)
            trx.append(tr.view(np.ndarray))
            self.data_array.append(tr.view(np.ndarray)[self.tracepos])
            self.cmd_wave.append(tr.view(np.ndarray)[self.cmdpos])
            if sequence_values is not None:
                self.values.append(sequence_values[i])
            self.time_base.append(tr.xvals('Time'))
            sr = info[1]['DAQ']['primary']['rate']
            self.sample_rate.append(self.samp_rate)
            #print ('i: %d   cmd: %f' % (i, sequence_values[i]*1e12))
        if self.mode is None:
            print ('no files processed...')
            return False
        if 'v' in self.mode.lower():
            units = 'V'
        else:
            units = 'A'
        self.traces = np.array(trx)
        if len(self.values) == 0:
            ntr = len(self.traces)
            self.traces = self.traces[:ntr]
            self.values = np.zeros(ntr) # fake 
        else:
            ntr = len(self.values)

        self.traces = EM.MetaArray(self.data_array,
            info=[{'name': 'Command', 'units': cmd.axisUnits(-1),
             'values': np.array(self.values)},
             tr.infoCopy('Time'), tr.infoCopy(-1)])
        self.cmd_wave = EM.MetaArray(self.cmd_wave,
             info=[{'name': 'Command', 'units': cmd.axisUnits(-1),
              'values': np.array(self.values)},
              tr.infoCopy('Time'), tr.infoCopy(-1)])
        self.sample_interval = 1./self.sample_rate[0]
        self.data_array = np.array(self.data_array)
        self.time_base = np.array(self.time_base[0])
        protoreps = ('protocol', 'repetitions')
        mclamppulses = ('MultiClamp1', 'Pulse_amplitude')
        seqparams = index['.']['sequenceParams']
        #self.printIndex(index)

        if mclamppulses in seqparams.keys():
            self.repetitions = len(seqparams[mclamppulses])
            self.commandLevels = np.array(seqparams[mclamppulses])
            function = index['.']['devices']['MultiClamp1']['waveGeneratorWidget']['function']
            stimuli = index['.']['devices']['MultiClamp1']['waveGeneratorWidget']['stimuli']
            self.tstart = stimuli['Pulse']['start']['value']
            self.tend = self.tstart + stimuli['Pulse']['length']['value']
            
            # print 'commandlevels: ', self.commandLevels
        elif protoreps in seqparams.keys():
            self.repetitions = seqparams[protoreps][0] + 1
        else:
            protoreps = 1
            self.repetitions = 1
        return True

    def getClampCommand(self, data, generateEmpty=True):    
        """Returns the command data from a clamp MetaArray.
        If there was no command specified, the function will 
        return all zeros if generateEmpty=True (default).
        """

        if data.hasColumn('Channel', 'Command'):  # hascolumn is a metaarray method
            return data['Channel': 'Command']
        elif data.hasColumn('Channel', 'command'):
            return data['Channel': 'command']
        else:
            if generateEmpty:
                tVals = data.xvals('Time')
             #   mode = getClampMode(data)
                print ('Mode: ', self.mode)
                if 'v' in self.mode.lower():
                    units = 'V'
                else:
                    units = 'A'
                return EM.MetaArray(np.zeros(tVals.shape), info=[{'name': 'Time', 'values': tVals, 'units': 's'}, {'units': units}])
        return None

    def getBlueLaserTimes(self):
        """
        Get laser pulse times  - handling multiple possible configurations (ugly)
        """
        supindex = self._readIndex(self.protocol)
        #print(supindex['.']['devices']['PockelCell']['channels']['Switch'].keys())
        try:
            stimuli = supindex['.']['devices']['Laser-Blue-raw']['channels']['pCell']
            # print('GOT PCELL')
        except:
            try: 
                stimuli = supindex['.']['devices']['PockelCell']['channels']['Switch']
                # print('GOT POCKEL CELL/Switch')
            except:
                print(supindex['.']['devices']['PockelCell']['channels'].keys())
                raise
        # print('skeys: ', stimuli.keys())
        stimuli = stimuli['waveGeneratorWidget']['stimuli']
        # print ('skeys2: ', stimuli.keys())
        if 'Pulse' in stimuli.keys():
            times = {}
            times['start'] = [stimuli['Pulse']['start']['value']]
            times['duration'] = stimuli['Pulse']['length']['value']
            times['amplitude'] = stimuli['Pulse']['amplitude']['value']
            times['type'] = stimuli['Pulse']['type']
            return times
        elif 'Pulse3' in stimuli.keys():
            times = {}
            times['start'] = [stimuli['Pulse3']['start']['value']]
            times['duration'] = stimuli['Pulse3']['length']['value']
            times['amplitude'] = stimuli['Pulse3']['amplitude']['value']
            times['type'] = stimuli['Pulse3']['type']
            return times
        elif 'PulseTrain' in stimuli.keys():
#            print('PTR: ', stimuli['PulseTrain'])
            times = {}
            times['start'] = []
            tstart = [stimuli['PulseTrain']['start']['value']]
            times['duration'] = stimuli['PulseTrain']['length']['value']
            times['amplitude'] = stimuli['PulseTrain']['amplitude']['value']
            times['npulses'] = stimuli['PulseTrain']['pulse_number']['value']
            times['period'] = stimuli['PulseTrain']['period']['value']
            times['type'] = stimuli['PulseTrain']['type']
#            print('times: ', times)
            for n in range(times['npulses']):
                times['start'].append(tstart[0] + n*times['period'])
            return times
        else:
            raise ValueError('need to find keys for stimulus (might be empty): ' % stimuli)

    def getDeviceData(self, device='Photodiode', devicename='Photodiode'):
        """
        Get the data from a device
        
        Parameters
        ----------
        device : str (default: 'Photodiode')
            The base name of the file holding the data. '.ma' will be appended
            to the name
        devicename : str (default: 'Photodiode')
            The name of the device as set in the config (might be 'pCell', etc)
            This might or might not be the same as the device
        
        Returns
        -------
        Success : boolean
        
        The results are stored data for the current protocol
        """ 
        # non threaded
        dirs = self.subDirs(self.protocol)
        index = self._readIndex()
        trx = []
        cmd = []
        sequence_values = None
        if 'sequenceParams' in index['.'].keys():
            self.sequence =  index['.']['sequenceParams']
        else:
            self.sequence = []
        # building command voltages or currents - get amplitudes to clamp

        reps = ('protocol', 'repetitions')
        foundLaser = False
        self.Device_data = []
        self.Device_sample_rate = []
        self.Device_time_base = []
        for i, d in enumerate(dirs):
            fn = os.path.join(d, device + '.ma')
            if not os.path.isfile(fn):
                print(' acq4read.getDeviceData: File not found: ', fn)
                return None
            lbr = EM.MetaArray(file=fn)
            info = lbr[0].infoCopy()
            self.Device_data.append(lbr.view(np.ndarray)[0])
            self.Device_time_base.append(lbr.xvals('Time'))
            sr = info[1]['DAQ'][devicename]['rate']
            self.Device_sample_rate.append(sr)
        self.Device_data = np.array(self.Device_data)
        self.Device_sample_rate = np.array(self.Device_sample_rate)
        self.Device_time_base = np.array(self.Device_time_base)
        return {'data': self.Device_data, 'time_base': self.Device_time_base, 'sample_rate': self.Device_sample_rate}

    def getLaserBlueCommand(self):
        """
        Get the command waveform for the blue laser
        data for the current protocol
        """ 
        # non threaded
        dirs = self.subDirs(self.protocol)
        index = self._readIndex()
        trx = []
        cmd = []
        sequence_values = None
        if 'sequenceParams' in index['.'].keys():
            self.sequence =  index['.']['sequenceParams']
        else:
            self.sequence = []
        # building command voltages or currents - get amplitudes to clamp

        reps = ('protocol', 'repetitions')
        foundLaser = False
        self.LaserBlueRaw = []
        self.LBR_sample_rate = []
        self.LBR_time_base = []
        for i, d in enumerate(dirs):
            fn = os.path.join(d, 'Laser-Blue-raw.ma')
            if not os.path.isfile(fn):
                print(' acq4read.getLaserBlueCommand: File not found: ', fn)
                return False
            lbr = EM.MetaArray(file=fn)
            info = lbr[0].infoCopy()
            self.LaserBlueRaw.append(lbr.view(np.ndarray)[0])
            self.LBR_time_base.append(lbr.xvals('Time'))
            try:
                sr = info[1]['DAQ']['Shutter']['rate']
            except:
                print(info[1]['DAQ'].keys())
                exit(1)
            self.LBR_sample_rate.append(sr)
        self.LaserBlueRaw = np.array(self.LaserBlueRaw)
        self.LBR_sample_rate = np.array(self.LBR_sample_rate)
        self.LBR_time_base = np.array(self.LBR_time_base)
        return True

    def getPhotodiode(self):
        """
        Get the command waveform for the blue laser
        data for the current protocol
        """ 
        # non threaded
        dirs = self.subDirs(self.protocol)
        index = self._readIndex()
        trx = []
        cmd = []
        sequence_values = None
        if 'sequenceParams' in index['.'].keys():
            self.sequence =  index['.']['sequenceParams']
        else:
            self.sequence = []
        # building command voltages or currents - get amplitudes to clamp
        reps = ('protocol', 'repetitions')
        foundPhotodiode = False
        self.Photodiode = []
        self.Photodiode_time_base = []
        self.Photodiode_sample_rate = []
        for i, d in enumerate(dirs):
            fn = os.path.join(d, 'Photodiode.ma')
            if not os.path.isfile(fn):
                print(' acq4read.getPhotodiode: File not found: ', fn)
                return False
            pdr = EM.MetaArray(file=fn)
            info = pdr[0].infoCopy()
            self.Photodiode.append(pdr.view(np.ndarray)[0])
            self.Photodiode_time_base.append(pdr.xvals('Time'))
            sr = info[1]['DAQ']['Photodiode']['rate']
            self.Photodiode_sample_rate.append(sr)
        self.Photodiode = np.array(self.Photodiode)
        self.Photodiode_sample_rate = np.array(self.Photodiode_sample_rate)
        self.Photodiode_time_base = np.array(self.Photodiode_time_base)
        return True

    def getBlueLaserShutter(self):
        supindex = self._readIndex(self.protocol)
        stimuli = supindex['.']['devices']['Laser-Blue-raw']['channels']['Shutter']['waveGeneratorWidget']['stimuli']
        times = []
        shutter = {}
        shutter['start'] = stimuli['Pulse']['start']['value']
        shutter['duration'] = stimuli['Pulse']['length']['value']
        shutter['type'] = stimuli['Pulse']['type']
        return shutter
        
    def getScannerPositions(self, dataname='Laser-Blue-raw.ma'):
        dirs = self.subDirs(self.protocol)
        self.scannerpositions = np.zeros((len(dirs), 2))
        self.targets = [[]]*len(dirs)
        self.spotsize = 0.
        rep = 0
        tar = 0
        supindex = self._readIndex(self.protocol)
        if not 'sequenceParams' in supindex['.'].keys():
            return(False)
        try:
            ntargets = len(supindex['.']['sequenceParams'][('Scanner', 'targets')])
        except:
            return(False)
        pars={}
        pars['sequence1'] = {}
        pars['sequence2'] = {}
        try:
            reps = supindex['.']['sequenceParams'][('protocol', 'repetitions')]
        except:
            reps = [0]  # just fill in one rep. SOme files may be missing the protocol/repetitions entry for some reason
        pars['sequence1']['index'] = reps
        pars['sequence2']['index'] = ntargets
        self.sequenceparams = pars
        self.scannerCamera = {}
        self.scannerinfo = {}
        for i, d in enumerate(dirs):
            index = self._readIndex(d)
            if 'Scanner' in index['.'].keys():
                self.scannerpositions[i] = index['.']['Scanner']['position']
                self.targets[i] = index['.'][('Scanner', 'targets')]
                self.spotsize = index['.']['Scanner']['spotSize']
                self.scannerinfo[(rep, tar)] = {'directory': d, 'rep': rep, 'pos': self.scannerpositions[i]}
            # elif ('Scanner', 'targets') in index['.']:
            #         print ('scanner targets: ', index['.'][('Scanner', 'targets')])
            #         self.scannerpositions[i] = index['.'][('Scanner', 'targets')]['position']
            #         self.targets[i] = index['.'][('Scanner', 'targets')]
            #         self.spotsize = index['.']['Scanner']['spotSize']
            else:
#                print('Scanner information not found in index: ', d, '\n', index['.'].keys())
                return False # protocol is short... 
#                self.scannerinfo[(rep, tar)] = {'directory': d, 'rep': rep, 'pos': self.scannerpositions[i]}
            if 'Camera' in supindex['.']['devices'].keys() and len(self.scannerCamera) == 0:  # read the camera outline
                cindex = self._readIndex(os.path.join(d, 'Camera'))
                self.scannerCamera = cindex
            else:
                pass
                
                
            tar = tar + 1
            if tar > ntargets:
                tar = 0
                rep = rep + 1
        return True # indicate protocol is all ok

    def getImage(self, filename):
        """
        getImage
        Returns the image file in the dataname
        Requires full path to the data
        """
        self.imageData = tf.imread(filename)
        # imageframe = EM.MetaArray(file=dataname)
#         img = imageframe.view(np.ndarray)
        return(self.imageData)

    def getAverageScannerImages(self, dataname='Camera/frames.ma', mode='average', firstonly=False, limit=None):
        """
        Average the images across the scanner camera files
        the images are collected into a stack prior to any operation
        
        Parameters
        ----------
        dataname : str (default: 'Camera/frames.ma')
            Name of the camera data file (metaarray format)
        
        mode : str (default: 'average')
            Operation to do on the collected images
            average : compute the average image
            max : compute the max projection across the stack
            std : compute the standard deviation across the stack
        
        Returns
        -------
            a single image frame that is the result of the specified operation

        """
        assert mode in ['average', 'max', 'std']
        print('average scanner images')
        dirs = self.subDirs(self.protocol)

        rep = 0
        tar = 0
        supindex = self._readIndex(self.protocol)
        ntargets = len(supindex['.']['sequenceParams'][('Scanner', 'targets')])
        pars={}
        pars['sequence1'] = {}
        pars['sequence2'] = {}
        try:
            reps = supindex['.']['sequenceParams'][('protocol', 'repetitions')]
        except:
            reps = [0]
        pars['sequence1']['index'] = reps
        pars['sequence2']['index'] = ntargets
        scannerImages = []
        self.sequenceparams = pars
        self.scannerinfo = {}
        if limit is None:
            nmax = len(dirs)
        else:
            nmax = min(limit, len(dirs))
        for i, d in enumerate(dirs):
            if i == nmax:  # check limit here first
                break
            index = self._readIndex(d)
            imageframe = EM.MetaArray(file=os.path.join(d, dataname))
            cindex = self._readIndex(os.path.join(d, 'Camera'))
            frsize = cindex['frames.ma']['region']
            binning = cindex['frames.ma']['binning']
           # print ('image shape: ', imageframe.shape)
            if imageframe.ndim == 3 and imageframe.shape[0] > 1:
                imageframed = imageframe[1]
            if imageframe.ndim == 3 and imageframe.shape[0] == 1:
                imageframed = imageframe[0]
            if firstonly:
                resultframe = imageframed.view(np.ndarray)
                return resultframe
            if i == 0:
                scannerImages = np.zeros((limit, int(frsize[2]/binning[0]), int(frsize[3]/binning[1])))
            scannerImages[i] = imageframed.view(np.ndarray)

        resultframe = np.zeros((scannerImages.shape[1], scannerImages.shape[2]))
        # simple maximum projection
        print('mode: %s' % mode)
        print('scanner images: ', scannerImages.shape)
        if mode == 'max':
            for i in range(scannerImages.shape[0]):
                resultframe = np.maximum(resultframe, scannerImages[i])
        elif mode == 'average':
            resultframe = np.mean(scannerImages, axis=0)
        elif mode == 'std':
            resultframe = np.std(scannerImages, axis=0)
        return resultframe

    def plotClampData(self, all=True):
        f, ax = mpl.subplots(2)
        if all:
            for i in range(len(self.data_array)):
                ax[0].plot(self.time_base, self.data_array[i])
                ax[1].plot(self.time_base, self.cmd_wave[i])
        else:
            ax[0].plot(self.time_base, np.array(self.data_array).mean(axis=0))
        mpl.show()
        
if __name__ == '__main__':
    # test on a big file
    a = Acq4Read()
    BRI = BR.BoundRect()
#    a.setProtocol('/Users/pbmanis/Documents/data/MRK_Pyramidal/2018.01.26_000/slice_000/cell_000/CCIV_1nA_max_000/')
    cell = '/Users/pbmanis/Documents/data/mrk/2017.09.12_000/slice_000/cell_001'
    datasets = os.listdir(cell)
    imageplotted = False
    imagetimes = []
    imagename = []
    maptimes = []
    mapname = []
    supindex = a.readDirIndex(currdir=cell)
    for k in supindex:
        if k.startswith('image_'):
            print('Found Image: ', k)
            imagetimes.append(supindex[k]['__timestamp__'])
            imagename.append(k)
        if k.startswith('Map_'):
            maptimes.append(supindex[k]['__timestamp__'])
            mapname.append(k)
    print(maptimes)
    print(imagetimes)
    maptoimage = {}
    for im, m in enumerate(maptimes):
        u = np.argmin(maptimes[im] - np.array(imagetimes))
        maptoimage[mapname[im]] = imagename[u]
            
    print (maptoimage)

    for i, d in enumerate(datasets):
        pa, da = os.path.split(d)
        if 'Map'  not in da:
            continue
        print('d: ', d)
        a.setProtocol(os.path.join(cell, d))
    #    a.setProtocol('/Volumes/Pegasus/ManisLab_Data3/Kasten, Michael/2017.11.20_000/slice_000/cell_000/CCIV_4nA_max_000')
        if not a.getScannerPositions():
           continue
        
    
        print( a.scannerCamera['frames.ma']['transform'])
        pos = a.scannerCamera['frames.ma']['transform']['pos']
        scale = a.scannerCamera['frames.ma']['transform']['scale']
        region = a.scannerCamera['frames.ma']['region']
        binning = a.scannerCamera['frames.ma']['binning']
        print('bining: ', binning)
        if a.spotsize is not None:
            print ('Spot Size: {0:0.3f} microns'.format(a.spotsize*1e6))
        else:
            a.spotsize=50.

        camerabox = [[pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[1]],
               [pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[3]],
               [pos[0] + scale[0]*region[2], pos[1] + scale[1]*region[3]],
               [pos[0] + scale[0]*region[2], pos[1] + scale[1]*region[1]],
               [pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[1]]
           ]
        scannerbox = BRI.getRectangle(a.scannerpositions)
        print(scannerbox)
        print(scannerbox.shape)
        fp = np.array([scannerbox[0][0], scannerbox[1][1]]).reshape(2,1)
        print(fp.shape)
        scannerbox = np.append(scannerbox, fp, axis=1)
        print(scannerbox)
    
        boxw = np.swapaxes(np.array(camerabox), 0, 1)
        print('camera box: ', boxw)
        scboxw = np.array(scannerbox)
        print('scanner box: ', scboxw)
        mpl.plot(boxw[0,:], boxw[1,:], linewidth=1.5)
        avgfr = a.getAverageScannerImages(firstonly=True, mode='average')
        if not imageplotted:
            imgd = a.getImage(os.path.join(cell, 'image_001.tif'))
           # mpl.imshow(np.flipud(np.rot90(avgfr), aspect='equal', extent=[np.min(boxw[0]), np.max(boxw[0]), np.min(boxw[1]), np.max(boxw[1])])
            mpl.imshow(imgd, aspect='equal', cmap='gist_gray',
                extent=[np.min(boxw[0]), np.max(boxw[0]), np.min(boxw[1]), np.max(boxw[1])])
            imageplotted = True
        mpl.plot(a.scannerpositions[:,0], a.scannerpositions[:,1], 'ro', alpha=0.2, markeredgecolor='w')
        mpl.plot(boxw[0,:], boxw[1,:], 'g-', linewidth=5)
        mpl.plot(scboxw[0,:], scboxw[1,:], linewidth=1.5, label=d.replace('_', '\_'))
    
    # a.getData()
    # a.plotClampData(all=True)
    #print a.clampInfo
   # print a.traces[0]
    pos = mpl.ginput(-1, show_clicks=True)
    print(pos)
    
    mpl.legend()
    mpl.show()
            

