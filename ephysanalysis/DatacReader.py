#!/usr/bin/env python
# encoding: utf-8
"""
readdatac.py
original C code by Daniel Bertrand and Paul Manis (1985-1999).
Matlab version by Patrick Kanold (1996).

Matlab: 4/30/99 P. Manis. Modified to handle modes 0 and 1 also, as well as get the timing information.
now we should be able to read legacy files, and files with just one channel in modes 2 and 3
Python: Converted to 5/13/2014 PBManis.
Tested on old data
Created by Paul Manis on 2014-05-14.
Copyright (c) 2014 Paul B. Manis, Ph.D.. All rights reserved.
"""

import sys
import os
import platform
import unittest
import struct
import re
import numpy as np
import glob
import matplotlib.pylab as MP
import getcomputer


datapath = {'Lytle': '/Volumes/ManisLab_Data3/Rothman_Jason/DATA/',
            'Tamalpais': '/Users/pbmanis/Documents/data/Rothman_Jason/DATA',
            'Tamalpais2': '/Users/pbmanis/Documents/data',}

# computername = platform.node()
# print computername
# sysname = None
# for k in datapath.keys():
#     if computername.find(k) >= 0: # just needs to be part of the name
#         sysname = k
bdir, sysname = getcomputer.getcomputer()
print sysname
if sysname is None:
    print ("Don't know path for system {:s}".format(computername))
    exit()
else:
    print("Found system: {:s}".format(sysname))
    
class ReadDatac:
    def __init__(self):
        # initialize a DATAC file object, used to read data and notefiles
        # dfile - the variables that are in dfile structure, as created by
        # the various versions of DATAC, are defined here
        
        nrl = 16  # 
        self.nrl = nrl
        RL = np.zeros((nrl,1))  # a list of records
        self.fullfile = ''  # full filename with path
        self.filename = ''  # filename, with extension
        self.path = ''  # separated path
        self.ext = ''  # extension (often this is the initials of the person doing the recording)
        self.records_in_file = 0  # number of records stored in the file
        self.nr_points = 0  # number of points per record
        self.comment = ' '  # the "comment" field (string)
        self.mode = -1  # file mode (-1 is placeholder). Mode varies between versions of the file structure
        self.dmode = 'CC'  # default data collection mode if not otherwise defined...
        self.rate = np.zeros(nrl)  # there is a separate sample rate for every record!
        self.ftime = 0.0  # time variable/sampling information
        self.record = []  # the records
        self.nr_channel = 0  # number of channels collected (usually 2 or 3; 3 is max)
        self.junctionpot = 0.0  # junction potential in mV for corrections (nominally 0 mV)
        self.vgain = 1.0  # voltage gain
        self.wgain = 1.0  # channel 3 gain
        self.igain = 1.0  # current gain
        self.gain = np.zeros((nrl,8))  # gain for each channel
        self.low_pass = np.zeros((nrl,8))  # low pass filter settings for each channel
        self.slow = np.zeros((nrl,1))  # slow mode flag (for slow acquisition rates on old AD board)
        self.ztime = np.zeros((nrl,1))  # time since file opened
        self.refvg = 1.0
        self.refig = 1.0
        self.refwg = 1.0
        self.frec = 1.0  # first record of data that is pulled
        self.lrec = 1.0  # last record
        self.steps = 0  # nubmer of steps in the stimulus
        self.fileheaderlen = 128  # fixed length file header, 128 bytes
        self.recheaderlen = 256  # fixed length of header for each record, 256 bytes
        self.datawidth = 2  # standard width of int variable, in bytes
        self.data = None  # the data that is returned
        self.fid = None  # the file id of an open file

    def openfile(self, filename):
        """
        open the requested file for reading
        :params: filename - the full name of the file including the path
        :returns: error  - integer error flag:
             1: requested file was not found or could not be opened
             3: requested file was a matlab file; should be handled in matlab
             0: normal return, file is open and fid is the open file handle.
        """
        
        if not os.path.exists(filename):
            print('ReadDatac:readfile File %s not found' % filename)
            self.err = 1
            return self.err

        # decide what kind of file we are reading
        p, ext = os.path.splitext(filename)
        path, fname = os.path.split(filename)
        if ext in  ['.mat']: # if its a matlab file, use different routine to parse through here
            #[dfile, data, err] = read_acq_mat(filename, RL)
           self.err = 3
           return self.err

        # otherwise, we proceed normally.

        self.fid = open(filename,'rb')  # open for readonly, as binary file
        if self.fid is None:
            print('datac2mat: Unable to open file: %s ',filename)
            self.err == 1
            return self.err

        self.fullfile = filename  # store file info
        self.filename = fname
        self.path = path
        self.ext = ext
        self.err = 0
        return self.err


    def readheader(self, filename, indent=0):
        """
        Read the file header information. Each file in the DATAC format has a header
        that contains basic information about the file, including a comment field,
        the file mode flag, and the number of channels and points per channel. 
        Depending on the mode, other information may also be set. 
        A mode of 0 is 2 channels, 1 is 3 channels. 
        Note that early DATAC files held data in integer format, as received from the data acquistion board.
        Later files (mode = 9) stored the data in a floating point format.

        :params: filename: name of the file to open and read the header.
        :return: nothing (self.err is set)
        """
        data=[]
        err=1
        endian='ieee-le' # data from little endian (i.e. PC's) => integers flipped, floats in ieee
#        dfile = init_dfile # create the current version of the dfile structure.
        #if self.nrl > 0:
        #    print('Opening %s' % filename)

        err = self.openfile(filename)
        if err > 0:
            self.err = err
            return       
        # before reading, check the length of the file
        self.fid.seek(0, 2)
        flen = self.fid.tell()
        if flen < 129:  # must have at least header and one record
            self.fid.close()
            # print('datac2mat: Empty file: %s ',filename)
            self.err = 2
            return self.err
        #read file header
        self.fid.seek(0)
        
        c_np = self.fid.read(2)
        self.nr_points = struct.unpack('H', c_np)[0]
        c_nothing = self.fid.read(27)
        c_comment = self.fid.read(99)
        cp = struct.unpack('99c', c_comment)
        self.comment = ''
        for c in cp:
            if c in ['\n', '\x00']:
                continue
            self.comment += c

        self.fid.seek(self.fileheaderlen, 0)
        c_mode = self.fid.read(1)
        self.mode = struct.unpack('b', c_mode)[0]
        c_time = self.fid.read(1)
        self.ftime = struct.unpack('b', c_time)[0]
        
        self.datawidth = 2  # default
        if(self.mode >= 2):
             c_record = self.fid.read(2)
             c_nrch = self.fid.read(2)

             self.record = struct.unpack('H', c_record)[0] #should be 1 for first record
             self.nr_channel = struct.unpack('H', c_nrch)[0]
        else:
            self.record = 1
            self.nr_channel = self.mode+2  # mode 0 is 2 channels; mode 1 is 3 channels
            self.nbytes = self.nr_points*8   # original form */
            self.nheader = 6       # original header was 6 bytes long */
            self.ndbytes = self.nbytes - self.nheader
            
        if self.mode == 9:
            self.datawidth = 4
            self.vgain = 0.1

        # check the file length now
        self.records_in_file = int(np.floor((flen-128)/(self.nr_points*self.nr_channel*self.datawidth+256.)))
        if self.records_in_file == 0:
            print('No records found in file!')
            self.err = 2
            self.fid.close()
            return self.err
        # print some info on the file that was just opened.
        spc = ' '*indent
        print('%sFile %s   Mode:         %5d       ftime: %6d   Channels:%6d  Data format: %d bytes' %
             (spc, self.filename, self.mode, self.ftime, self.nr_channel, self.datawidth))
        print('%s  First Record: %5d     Num Records: %6d' % (spc, self.record, self.records_in_file))
        
        self.err = 0
        return self.err

    def readnote(self):
        """
        Read the note file associated with the current open datafile
        The note information (at least, some of it) is parsed and stored in 
        the self.expInfo dictionary, numbered by the acquisition blocks
        """
        nf, ext = os.path.splitext(self.fullfile)
        self.notefile = nf + '.NOT'
        print self.notefile
        
        self.expInfo = {}
        with open(self.notefile) as f:
            content = f.readlines()
        f.close()
        k = 0
        protocol = None
        acqflag = False

        for lineno, c in enumerate(content):
            if lineno < 6:
                continue # skip over the header
            if c[0:3] == 'PRI':
                continue  # older ones, some lines have sequence information at start of line.
#            if c[20] == 'L':
#                continue
            t = c[0:8]  # get the time stamp for each line
            try:
                rec = int(c[12:16])  # read the record information
            except:
                continue # some lines in some versions can't be parsed this way; just skip.
            if acqflag:  # retrieve the end record after an acquistion, and store the information.
                caret = c.find('>')   # any command
                if caret != -1:
                    self.expInfo[k] = {'R': firstrec, 'Rend': rec, 'time': t, 'Protocol': protocol}
                    acqflag = False  # no longer acquiring
                    k += 1
            for j, stm in enumerate(['>do STM:','>g STM:']):
                cmd = c.find(stm)  # get the "do" command (get and seq)
                if cmd != -1:
                    if not acqflag:
                        protocol = c[cmd+8:cmd+18].strip()  # get the current protocol.
                    firstrec = rec
                    if acqflag is False and k >= 1:
                        self.expInfo[k-1]['Rend'] = rec  # set last record
                    if j == 0:
                        acqflag = True
#            if not acqflag:
            seq = c.find('>seq')  # check for an seq command from a protocol that is loaded
            if seq != -1:
                acqflag = True
            
    def printNoteInfo(self, indent):
        """
        Print the note file information in a table format.
        """
        spc = ' ' * indent
        k = self.expInfo.keys()
        k.sort()
        print('{:s}  Block  Protocol   Recs   Time'.format(spc))
        for i in k:
            print('{:s}  {:3d}  {:9s} {:3d}-{:3d} {:8s} '.format(spc, i, self.expInfo[i]['Protocol'], 
                self.expInfo[i]['R'], self.expInfo[i]['Rend'],
                self.expInfo[i]['time']))
        
    def readrecords(self, RL):
        """
        Read the selected data records from the currently open file. 
        :params: RL - the list of records (does not need to be contiguous)
        Returns: Error (self.err)
            Errors:
            0: success
            3: last record is past end of data in file
            
        """ 
        if len(RL) == 0: # catch when we just read the header for information
            self.err = 0; # this is ok...
            return self.err

        if RL[0] == 0:
            RL = [x+1 for x in RL]
        maxRL = np.max(RL)
        if maxRL > self.records_in_file:
            print('Last record (%d) greater than length of file (%d recs in file)' % (maxRL,self.records_in_file))
            #fid.close(fid)
            self.err = 3
            return self.err
        block_head=0; # flag for block header reading ONLY.
        #read data sets according to the records in the RL vector
        self.nrl = len(RL)
        self.frec = RL[0]
        self.lrec = RL[-1]
        
        if maxRL == 0: # special mode just to read ZTIME and other header information from WHOLE FILE
            RL = range(1,self.records_in_file)
            block_head = 1 # block header read
            self.data = None
            self.rate = None
        else:
            self.data = np.zeros((self.nrl, self.nr_channel, self.nr_points)); # go ahead and allocate memory.
            self.rate = np.zeros(self.nrl); # separate for every record!

#        data=np.zeros((self.nrl,4)) # just initialize to make rest of program happy - adjust later when we know
        self.record = [None]*self.records_in_file
        self.channels = [None]*self.records_in_file
        self.rate = [None]*self.records_in_file
        self.slow = [None]*self.records_in_file
        self.ztime = [None]*self.records_in_file
        self.gain = np.ones((self.records_in_file, 8))
        self.low_pass = np.ones((self.records_in_file, 8))
        
        for i, rec in enumerate(RL):
            if self.dmode >= 2:
                offset = self.fileheaderlen+(rec-1) * self.nr_points*self.datawidth*self.nr_channel+(rec-1)*self.recheaderlen; #jump to data start
            else:
                offset = self.fileheaderlen+(rec-1)*self.ndbytes
            self.fid.seek(offset, 0)      # skip over earlier records
            # read record header
            c_time = self.fid.read(1)  # reads the mode byte, but we don't need it.
            self.ftime = struct.unpack('b', c_time)[0]
            if self.mode == 9:
                datawid = 4
            else:
                datawid = 2
            if self.mode < 2:
                # print 'self.mode (< 2): ', self.mode
                self.gain = np.ones((self.nrl,8)) # set the gains all to 1
                offset = self.fileheaderlen+(rec-1)*self.ndbytes+rec*self.nheader   # position it correctly...
                self.fid.seek(offset, 0)  # skip earlier records
                if block_head == 0: # read the data, not just the block header..
                    c_data_in = self.fid.read(self.ndbytes*datawid)  # access the data itself
                    try:
                        data_in = struct.unpack('%dh' % int(self.ndbytes), c_data_in)
                    except:
                        print 'Failed to read record %d' % (rec)
                        print 'self.ndbytes: ', self.ndbytes
                        print 'len cdata: ', len(c_data_in)
                        print 'nr points: ', self.nr_points
                        break
                        
                    if self.mode == 0: # 3 channels of information
                        self.data[i,0,:] = data_in[1:self.nr_points*3:self.nr_channel+1] 
                        self.data[i,1,:] = data_in[2:self.nr_points*3:self.nr_channel+1] 
                        self.rate[i] = self.oldtime(data_in[0:self.nr_points*3:3])
                    else:
                        self.data[i,0,:] = datain[1::self.nr_channel+1]
                        self.data[i,1,:] = data_in[2::self.nr_channel+1]
                        self.rate[i] = self.oldtime(data_in[0::self.nr_channel+1])
        
            else:
                # print 'data mode > 2 : %d' % self.mode
                c_rec = self.fid.read(2)
                self.record[i] = struct.unpack('h', c_rec)[0]
                c_chan = self.fid.read(2)
                self.channels[i] = struct.unpack('H', c_chan)[0]  # fread(fid,1,'int16')
                c_rate = self.fid.read(4)
                self.rate[i] = struct.unpack('f', c_rate)[0]  # fread(fid,1,'float32')
                c_gain = self.fid.read(8*4)
                self.gain[i,:] = struct.unpack('8f', c_gain)  # fread(fid,8,'float32')
                c_lpf = self.fid.read(8*4)
                self.low_pass[i,:] = struct.unpack('8f', c_lpf)  #fread(fid,8,'float32')
                c_slow = self.fid.read(2)
                self.slow[i] = struct.unpack('H', c_slow)  # fread(fid,1,'int16')
                c_ztime = self.fid.read(4)
                self.ztime[i] = struct.unpack('I', c_ztime)[0]  # fread(fid,1,'long')
                # now read data
                offset = 128+(rec-1)*datawid*self.nr_points*self.nr_channel+rec*256; #jump to data start
        
                self.fid.seek(offset, 0) #skip earlier records
                if block_head == 0: # we don't actually read the data...
                    if self.mode != 9:
                        c_data = self.fid.read(self.nr_channel*self.nr_points*datawid) 
                        data_in = struct.unpack('%dh' % (self.nr_channel*self.nr_points), c_data) 
                    else:  # floating point format
                        c_data = self.fid.read(self.nr_channel*self.nr_points*datawid)  
                        data_in = struct.unpack('%dh' % (self.nr_channel*self.nr_points), c_data)
                    
                    skip = self.nr_channel # set the skip counter
                    # read the first entry
                    self.data[i,0,:] = data_in[0::2]     # self.nr_points - 1 ??? #voltage
                    if skip > 1:  # get second channel
                        self.data[i,1,:] = data_in[1::2]   #current
                    if skip > 2:  # get third channel ('w')
                        self.data[i,2,:] = data_in[2::2]  #wchannel...
        if self.nrl == 0:  # we sometimes override these with the ctl structure...
            self.refvgain = self.gain(1,1)
            self.refigain = self.gain(1,2)
            self.refwgain = self.gain(1,3)
        self.data = np.flip(self.data, axis=2)
        #self.fid.close()
        self.err = 0
        return self.err

    def close(self):
        """
        close the current open data file
        """
        self.fid.close()

    def oldtime(self, tbuf):
        """
        compute mean rate from the old timey (mode 0 and 1) modes, associated with 
        the old switching amplifier from Prof. Daniel Bertrand, CMU Geneva.
        """
        CTE = 65536
        FTIME = 0.00061
        tx = 0
        dt = FTIME
        tdt = dt*10
        t = tbuf[0]
        ta = dt * (CTE-t)
        if ta < -tdt:
            tx = tx + CTE
            ta = dt * (tx-t)
        rate = ta
        #print('Rate: %12.5f' % rate)
        return(rate)


    def plotcurrentrecs(self, fn, title=None):
        f = MP.figure(fn)
        if title is not None:
            f.suptitle(title)
        ch = 0
        sp = {}
        sp[0] = MP.subplot(211)
        sp[1] = MP.subplot(212)
        if self.data is None:
            return
        for i in range(self.data.shape[0]):
            if self.rate[i] == 0.:
                continue
            tb = np.arange(0., self.nr_points/self.rate[i], 1./self.rate[i])
 
            sp[0].plot(tb, (self.data[i,ch,:]-2047)*0.1339, 'k-')
            sp[1].plot(tb, (self.data[i,ch+1,:]-2047)*-1., 'r-')
        #MP.show()

class GetClamps():
    
    def __init__(self, datac, path):
        self.datac = datac

    def getClampData(self, chmap, block):
        """
        Translates fields as best as we can from the original DATAC structure
        create a Clamp structure for use in SpikeAnalysis and RMTauAnalysis.
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
        protocol = item.text()
        
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

        self.traces = MetaArray(self.traces, info=info)
        self.spikecount = np.zeros(len(recs))
        self.rgnrmp = [0, 0.005]

class untitledTests(unittest.TestCase):
    def setUp(self):
        pass

def printNotes(filename, indent=0):
    dfile = ReadDatac()
    #dfile.readheader('/Volumes/ManisLab_Data3/VCN1/SOM/6DEC88B.SOM')
    dfile.readheader(filename, indent)
    if dfile.err == 2:
        return 1
    dfile.readnote()
    dfile.printNoteInfo(indent)
    return 0
    
def plotonefile(filename):
    dfile = ReadDatac()
    #dfile.readheader('/Volumes/ManisLab_Data3/VCN1/SOM/6DEC88B.SOM')
    dfile.readheader(filename)
    dfile.readnote()
    dfile.printNoteInfo(3)
    if dfile.err == 0:
        k = dfile.expInfo.keys()
        for i in k:
            dfile.readrecords(range(dfile.expInfo[i]['R'],dfile.expInfo[i]['Rend']))
            if dfile.err == 0:
                txt = dfile.expInfo[i]['Protocol'] + ' ' '[R:%d-%d] ' % (dfile.expInfo[i]['R'],dfile.expInfo[i]['Rend'])
                txt += dfile.expInfo[i]['time']
                dfile.plotcurrentrecs(i, title=txt)
                
    dfile.close()
    MP.show()

def listAllFiles():
    n = 0
    nempty = 0
    for d in os.listdir(datapath):
        #f = glob.glob(os.listdir(os.path.join(datapath, d) + '.JSR'))
        p = os.path.join(datapath, d)
        if not os.path.isdir(p):
            continue
        print '\nDirectory: %s' % p
        for fx in os.listdir(p): 
            (f, e) = os.path.splitext(fx)
            ff = os.path.join(p, fx)
            if os.path.isfile(ff) and e == '.JSR':
                nempty += printNotes(ff, indent=3)
                n += 1
    print('\nFiles read: {:d}, Empty file accessed: {:d}'.format(n, nempty))

if __name__ == '__main__':

   # plotonefile(os.path.join('/Users/pbmanis/Documents/data/HWF0001B/VCN', '11SEP96H.HWF'))
    plotonefile(os.path.join('/Users/pbmanis/Documents/data/datac', '07FEB96H.JSR'))