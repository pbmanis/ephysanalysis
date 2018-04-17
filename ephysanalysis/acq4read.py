#!/usr/bin/python

"""
Class to read acq4 data blocks in simple manner. Does not requre
acq4 link; bypasses DataManager and PatchEPhys

Requires pyqtgraph to read the .ma files and the .index file

"""
import os
from pyqtgraph import metaarray
from pyqtgraph import configfile
import numpy as np
import matplotlib.pyplot as mpl
import pprint

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
        
    def setProtocol(self, pathtoprotocol):
        """
        Parameters
        ----------
        pathtoprotocol str (default: None)
            Path to the protocol directoryto set for this instance of the reader
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
                print '****************** Error: Missing .index file? (fails to detect protocol sequence)'
                raise Exception("Directory '%s' does not appear to be a protocol sequence." % dh.name())

    def getIndex(self, currdir=''):
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

    def printIndex(self, index):
        """
        Generate a nice printout of the index, about as far down as we can go
        """
        for k in index['.'].keys():
            print '  ', k, ':  ', index['.'][k]
            if isinstance(index['.'][k], dict):
                for k2 in index['.'][k].keys():
                    print '    ', k, ' ', k2, '::  ', index['.'][k][k2]
                    if isinstance(index['.'][k][k2], dict):
                        for k3 in index['.'][k][k2]:
                            print '    ', k, ' ', k2, ' ', k3, ':::  ', index['.'][k][k2][k3]
                            if isinstance(index['.'][k][k2][k3], dict):
                                for k4 in index['.'][k][k2][k3]:
                                    print '    [', k, '][', k2, '][', k3, '][', k4, '] ::::  ', index['.'][k][k2][k3][k4]

    def getClampDevices(self, currdir=''):
        """
        Search for a known clamp device in the list of devices 
        used in the current protocol directory...
        
        Return
        ------
        list of valid clamp devices found (there may be more than one)
            List will be empty if no recognized device is found.
        """
        info = self.getIndex(currdir=currdir)
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
                tr = metaarray.MetaArray(file=fn, readAllData=False)
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
      #  print('MODE: ', self.mode)
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
        trx = []
        cmd = []
        sequence_values = None
        self.sequence =  index['.']['sequenceParams']
        # building command voltages or currents - get amplitudes to clamp

        reps = ('protocol', 'repetitions')
        foundclamp = False
        sequence_values = None
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
            tr = metaarray.MetaArray(file=fn)
            info = tr[0].infoCopy()
            self.parseClampInfo(info)
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
            exit(1)
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

        self.traces = metaarray.MetaArray(self.data_array,
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
                return MetaArray(np.zeros(tVals.shape), info=[{'name': 'Time', 'values': tVals, 'units': 's'}, {'units': units}])
        return None 

    def getScannerPositions(self, dataname='Laser-Blue-raw.ma'):
        dirs = self.subDirs(self.protocol)
        self.scannerpositions = np.zeros((len(dirs), 2))
        self.targets = [[]]*len(dirs)
        self.spotsize = 0.
        rep = 0
        tar = 0
        supindex = self._readIndex(self.protocol)
        ntargets = len(supindex['.']['sequenceParams'][('Scanner', 'targets')])
        pars={}
        pars['sequence1'] = {}
        pars['sequence2'] = {}
        reps = supindex['.']['sequenceParams'][('protocol', 'repetitions')]
        pars['sequence1']['index'] = reps
        pars['sequence2']['index'] = ntargets
        self.sequenceparams = pars
        self.scannerinfo = {}
        for i, d in enumerate(dirs):
            index = self._readIndex(d)
            if 'Scanner' in index['.'].keys():
                self.scannerpositions[i] = index['.']['Scanner']['position']
                self.targets[i] = index['.'][('Scanner', 'targets')]
                self.spotsize = index['.']['Scanner']['spotSize']
                self.scannerinfo[(rep, tar)] = {'directory': d, 'rep': rep, 'pos': self.scannerpositions[i]}
            else:
                self.scannerpositions[i] = [0., 0.]
                self.targets[i] = None
                self.spotsize = None
                self.scannerinfo[(rep, tar)] = {'directory': d, 'rep': rep, 'pos': self.scannerpositions[i]}
            tar = tar + 1
            if tar > ntargets:
                tar = 0
                rep = rep + 1

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
    a.setProtocol('/Users/pbmanis/Documents/data/MRK_Pyramidal/2018.01.26_000/slice_000/cell_000/CCIV_1nA_max_000/')
#    a.setProtocol('/Volumes/Pegasus/ManisLab_Data3/Kasten, Michael/2017.11.20_000/slice_000/cell_000/CCIV_4nA_max_000')
    # a.getScannerPositions()
    # print a.scannerpositions
    # print (a.spotsize)
#    mpl.plot(a.scannerpositions[:,0], a.scannerpositions[:,1], 'ro')
    a.getData()
    a.plotClampData(all=True)
    #print a.clampInfo
   # print a.traces[0]
    mpl.show()
            

