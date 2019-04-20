"""
Analyze EPSCs or IPSCs
Or EPSPs and IPSPs...

This module provides the following analyses:

1. Amplitudes from a train
2. Paired pulse facilitation for pulse pairs, and the first pair in a train.
3. Current-voltage relationship in voltage clamp measured over a time window

The results of the analysis are stored in the class variable analysis_summary


"""

import sys
from pathlib import Path
import os  #legacy
import matplotlib
import scipy.signal
import pandas as pd
matplotlib.use('Qt4Agg')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#rcParams['font.sans-serif'] = ['Arial']
#rcParams['font.family'] = 'sans-serif'
rc('text', usetex=True)
rcParams = matplotlib.rcParams
rcParams['svg.fonttype'] = 'none' # No text as paths. Assume font installed.
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['text.latex.unicode'] = True

import numpy as np

import ephysanalysis as EP
import ephysanalysis.metaarray as EM  # need to use this version for Python 3
import ephysanalysis.cursor_plot as CP
import matplotlib.pyplot as mpl
import matplotlib.colors
import pylibrary.PlotHelpers as PH

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph.Point import Point

def make_key(pathname):
    """
    Make a key string using the date, slice, cell and protocol from the path name
    """
    p = pathname.parts

    return(str('~'.join([p[i] for i in range(-4, 0)])))


class PSCSummary():
    def __init__(self, datapath, plot=True, update_regions=False):
        self.datapath = datapath
        self.AR = EP.acq4read.Acq4Read()  # make our own private cersion of the analysis and reader
        self.plot = plot
        self.db = None
        self.db_filename = None
        self.update_regions = update_regions

    def setup(self, clamps=None, spikes=None, baseline=[0, 0.001]):
        """
        Set up for the fitting
        
        Parameters
        ----------
        clamps : A datamodel structure (required)
            Brings the data to the module. This usually will be a PatchEphys object.
        
        spikes : A spikeAnalysis structure (optional)
            Has information about which traces have spikes
            Use this when analyzing events that may be contaminated by spikes
        
        baseline : list (2 elements)
            times over which baseline is measured (in seconds)
        
        """
        
        if clamps is None:
            raise ValueError("VC analysis requires defined clamps ")
        self.Clamps = clamps
        self.spikes = spikes
        self.set_baseline_times(baseline)

        self.analysis_summary = {}  # init the result structure
    
    def read_database(self, filename):
        self.db_filename = Path(filename)
        if self.db_filename.is_file():
            with(open(self.db_filename, 'rb')) as fh:
                self.db = pd.read_pickle(fh)
        else:
            self.db = pd.DataFrame(columns=['date', 'protocol', 'T0', 'T1'])

    def update_database(self):
        if self.db is not None:
            self.db.to_pickle(self.db_filename)

    def compute_PSC_IV(self, protocolName=None, plot=True, savetimes=False):
        """
        Simple plot voltage clamp traces
        """
        #print('path: ', self.datapath)
        self.AR.setProtocol(self.datapath)  # define the protocol path where the data is
        self.setup(clamps=self.AR)
        self.read_database(f"{protocolName:s}.p")
        # print('database has:')
        # print(self.db['date'])
        # print('protocol: ', protocolName)
        if self.AR.getData():  # get that data.
            ok = False
            if protocolName.startswith('Stim_IO'):
                ok = self.analyze_IO()
            elif protocolName.startswith('VC-EPSC_3'):
                ok = self.analyze_VDEP()
            elif protocolName.startswith('PPF'):
                ok = self.analyze_PPF()
            if not ok:
                print('Failed on protocol: ', self.datapath, protocolName)
                return False
            if plot:
                self.plot_vciv()
            if savetimes:
                date = make_key(self.datapath)
                # print('date: ', date)
                # print('db: ', self.db)
                # print('date list: \n', self.db['date'].tolist())
                if date not in self.db['date'].tolist():
                    self.db.loc[len(self.db)] = [date, protocolName, self.T0, self.T1]
                    print('new date added')
                else:
                    self.db.loc[date, 'date'] = date
                    self.db.loc[date, 'protocol'] = protocolName
                    self.db.loc[date, 'T0'] = self.T0
                    self.db.loc[date, 'T1'] = self.T1
                    print('old date updated')
                self.update_database()
                # print('db head: ', self.db.head())
            return True
        else:
            return False

    def get_stimtimes(self):
        pass

    def set_baseline_times(self, baseline):
        """
        baseline: 2-element list or numpy array
        """
        if len(baseline) != 2:
            raise ValueError('Baseline must be a 2-element array')
        if isinstance(baseline, list):
            baseline = np.array(baseline)
        self.baseline = np.sort(baseline)

    def get_baseline(self):
        """ return the mean values in the data
        """
        bl = self.mean_I_analysis(region=self.baseline, reps=[0])
        return bl
        
    def analyze_IO(self, rmpregion=[0., 0.05], protocolName=None):
 
        ptrain = self.AR.getStim('Stim0')
        # stim dict in ptrain will look like:
        # {'start': [0.05, 0.1], 'duration': [0.0001, 0.0001],
        # 'amplitude': [0.00025, 0.00025], 'npulses': [2], 'period': [0.05], 'type': ['pulseTrain']}
        # try:
        dd  = self.AR.getDeviceData(device='Stim0', devicename='command')
        if dd is None:
            return False
        # except:
        #     print('Failed to get device data from Stim0')
        #     return False
        filekey = Path(make_key(self.datapath))
        # check the db to see if we have parameters already
        dfiles = self.db['date'].tolist()
        if filekey in dfiles:
            delay = self.db.loc[filekey, 'date']['T0']
            t1    = self.db.loc[filekey, 'date']['T1']
            width = t1-delay
        else:
            delay = 1.0*1e-3
            width = 15.0*1e-3            
        reps = self.AR.sequence[('protocol', 'repetitions')]
        stim_dt = np.diff(ptrain['start'])
        stim_I = [ptrain['amplitude'][0]]
        mode = '?'
        if not ('Stim0', 'command.PulseTrain_amplitude') in self.AR.sequence.keys():
            raise ValueError('Cannot find PulseTrain_amplitude in stimulus command')
            
        stim_I = self.AR.sequence[('Stim0', 'command.PulseTrain_amplitude')]

        baseline = []
        meani = []
        stimamp = []
        stimintvl = []
        cmdv = []
        self.sign = -1
        self.i_mean = []
        self.set_baseline_times(rmpregion)
        self.analysis_summary['iHold'] = []
        self.i_mean = []
        self.analysis_summary[f'PSP_IO'] = [[]]*len(ptrain['start'])
        bl = self.mean_I_analysis(region=self.baseline, mode='baseline', reps=[0])
        for i in range(len(ptrain['start'])):
            pdelay = ptrain['start'][i] + delay
            if i == 0 and self.update_regions:
                rgn = self.set_region([ptrain['start'][i], ptrain['start'][i]+0.050], baseline=bl)
            else:
                rgn = [delay, delay+width]
            self.T0 = rgn[0]
            self.T1 = rgn[1]
            print('pdelay: ', pdelay, self.baseline)
            i_mean = self.mean_I_analysis(region=np.array(rgn)+ptrain['start'][i], mode='min', baseline=bl, reps=reps)
            if i_mean is None:
                return False
            self.analysis_summary[f'PSP_IO'][i] = self.sign*1e12*i_mean
            cmdv.append(self.i_mean_cmd)
            stimamp.append(ptrain['amplitude'][i])
            stimintvl.append(ptrain['period'][0])

        self.analysis_summary['psc_stim_amplitudes'] = 1e6*np.array(stim_I)
        self.analysis_summary['psc_intervals'] = np.array(stimintvl)
        self.analysis_summary['ppf_dt'] = np.array(stim_dt)
        self.analysis_summary['stim_times'] = ptrain['start']
        self.analysis_summary['window'] = [self.T0, self.T1]
        return True

    def analyze_VDEP(self, rmpregion=[0., 0.05], protocolName=None):
        """
        Analyze the voltage-dependence of EPSCs 
        
        When selecting the analysis window, choose a window that encompases
        the peak of the inward EPSC in the negative voltage range.
        """
        ptrain = self.AR.getStim('Stim0')
        # stim dict in ptrain will look like:
        # {'start': [0.05, 0.1], 'duration': [0.0001, 0.0001],
        # 'amplitude': [0.00025, 0.00025], 'npulses': [2], 'period': [0.05], 'type': ['pulseTrain']}
        dd  = self.AR.getDeviceData(device='Stim0', devicename='command')
        # print('ptrain: ', ptrain)
        # print('dd: ', dd)
        reps = self.AR.sequence[('protocol', 'repetitions')]
        stim_dt = np.diff(ptrain['start'])
        stim_I = [ptrain['amplitude'][0]]
        mode = '?'
        if not ('MultiClamp1', 'Pulse_amplitude') in self.AR.sequence.keys():
            raise ValueError('Cannot find (MultiClamp1, Pulse_amplitude) in stimulus command')
            
        stim_V = self.AR.sequence[('MultiClamp1', 'Pulse_amplitude')]

        delay = 6*1e-3
        ndelay = 25*1e-3  # 25 msec after stimulus seems a fair time.
        nwidth = 2.5*1e-3
        width = 0.25*1e-3
        filekey = Path(make_key(self.datapath))
        # check the db to see if we have parameters already
        dfiles = self.db['date'].tolist()
        if filekey in dfiles:
            delay = self.db.loc[filekey, 'date']['T0']
            t1    = self.db.loc[filekey, 'date']['T1']
            width = t1-delay
        else:
            delay = 1.0*1e-3
            width = 20.0*1e-3

        bl_region = [ptrain['start'][0]-0.060, ptrain['start'][0]-0.010]  # time just before stimulus
        baseline = []
        self.baseline = bl_region
        meani = []
        stimamp = []
        stimintvl = []
        cmdv = []
        self.sign = 1
        self.i_mean = []
        # self.set_baseline_times(rmpregion)
        self.analysis_summary['iHold'] = []
        self.i_mean = []
        self.analysis_summary[f'PSP_VDEP'] = [[]]*len(ptrain['start'])
        self.analysis_summary[f'PSP_VDEP_NMDA'] = [[]]*len(ptrain['start'])
        bl = self.mean_I_analysis(region=bl_region, mode='baseline', reps=[0])
        print('bl_region: ', bl_region)
        for i in range(len(ptrain['start'])):
            pdelay = ptrain['start'][i] + delay
            if i == 0 and self.update_regions:
                rgn = self.set_region([ptrain['start'][i]+0.002, ptrain['start'][i]+0.025], baseline=bl, slope=False)
            else:
                rgn = [delay, delay+width]
            self.T0 = rgn[0]
            self.T1 = rgn[1]
            pdelay = ptrain['start'][i] + self.T0
            pdelay2 = ptrain['start'][i] + self.T1
            nmdelay = ptrain['start'][i] + ndelay
            i_mean = self.mean_I_analysis(region=[pdelay, pdelay2], mode='mean',
                                        baseline=bl, reps=reps, slope=False)
            if i_mean is None:
                return False
            # values for nmda analysis are currently fixed
            i_nmda_mean = self.mean_I_analysis(region=[nmdelay-nwidth, nmdelay+nwidth], mode='mean',
                                        baseline=bl, reps=reps, slope=False)

            self.analysis_summary[f'PSP_VDEP'][i] = self.sign*i_mean
            self.analysis_summary[f'PSP_VDEP_NMDA'][i] = self.sign*i_nmda_mean
            cmdv.extend(self.i_mean_cmd)
            stimamp.append(ptrain['amplitude'][i])
            stimintvl.append(ptrain['period'][0])
        
        data1, tb = self.get_traces(region=[ptrain['start'][0], ptrain['start'][0]+0.050],
            trlist=None, baseline=bl, intno=0, nint=1, reps=reps, slope=False)
        # find -80 and +30 voltage indices (so we can save them and save the data)
        cmdv = np.array(cmdv)+self.AR.holding
        i80 = np.argmin(np.fabs(cmdv+0.080))
        i30 = np.argmin(np.fabs(cmdv-0.030))
        if i30 >= data1.shape[0]:

            self.analysis_summary['Vindices'] = {'-80': np.nan, '30': np.nan}
            self.analysis_summary['AMPA_NMDA_traces'] = {'T': None, 'VN80': None, 'VP30': None}
        else:
            print('data1 shape: ', data1.shape, i80, i30, cmdv[i80], cmdv[i30])

            self.analysis_summary['Vindices'] = {'-80': i80, '30': i30}
            self.analysis_summary['AMPA_NMDA_traces'] = {'T': tb, 'VN80': data1[i80], 'VP30': data1[i30]}
        self.analysis_summary['psc_stim_amplitudes'] = np.array(stim_I)
        self.analysis_summary['psc_intervals'] = np.array(stimintvl)
        self.analysis_summary['stim_times'] = ptrain['start']
        self.analysis_summary['Vcmd'] = cmdv
        self.analysis_summary['window'] = [self.T0, self.T1]
        return True
        
    def analyze_PPF(self, rmpregion=[0., 0.05], protocolName=None):
        #self.rmp_analysis(region=rmpregion)
#        self.tau_membrane(region=tauregion)
        # r0 = self.Clamps.tstart + 0.9*(self.Clamps.tend-self.Clamps.tstart) #
        # print(dir(self.Clamps))
        ptrain = self.AR.getStim('Stim0')
        # stim dict in ptrain will look like:
        # {'start': [0.05, 0.1], 'duration': [0.0001, 0.0001], 'amplitude': [0.00025, 0.00025], 'npulses': [2], 'period': [0.05], 'type': ['pulseTrain']}
        dd  = self.AR.getDeviceData(device='Stim0', devicename='command')

        reps = self.AR.sequence[('protocol', 'repetitions')]
        stim_I = [ptrain['amplitude'][0]]
        if not ('Stim0', 'command.PulseTrain_period') in self.AR.sequence.keys():
            raise ValueError('Cannot find PulseTrain_period in stimulus command')
        
        stim_dt = self.AR.sequence[('Stim0', 'command.PulseTrain_period')]
        mode = 'PPF'
        print('PPF')
        delay = 1.0*1e-3
        width = 20.0*1e-3
        ndelay = 1.0*1e-3

        filekey = Path(make_key(self.datapath))
        # check the db to see if we have parameters already
        dfiles = self.db['date'].tolist()
        if filekey in dfiles:
            delay = self.db.loc[filekey, 'date']['T0']
            t1    = self.db.loc[filekey, 'date']['T1']
            width = t1-delay
        else:
            delay = 1.0*1e-3
            width = 20.0*1e-3
        nwidth = width
        baseline = []
        meani = []
        stimamp = []
        stimintvl = []
        cmdv = []
        self.sign = 1
        self.set_baseline_times(rmpregion)
        self.i_mean = []
        self.analysis_summary['iHold'] = []
        self.i_mean = []
        self.analysis_summary[f'PPF'] = [[]]*len(stim_dt)
        # print(ptrain)
        self.i_mean = []
        bl = self.mean_I_analysis(region=self.baseline, mode='baseline', reps=[0])
        for i in range(len(stim_dt)):
            if i == 0:
                rgn = [delay, delay+width]

                if self.update_regions:
                    rgn = self.set_region([ptrain['start'][0], ptrain['start'][0]+0.020], baseline=bl)
                self.T0 = float(rgn[0])
                self.T1 = float(rgn[1])
            p2delay = self.T0 + ptrain['start'][0] + stim_dt[i] # second stim of pair
            p1delay = self.T0 + ptrain['start'][0] # first stim of pair
            # print('p1, p2, t1: ', p1delay, p2delay, self.T1)
            # print('T1: ', float(self.T1))
            i_mean_ref = self.mean_I_analysis(region=[p1delay, p1delay+self.T1], mode='min', baseline=bl,
                        intno=i, nint=len(stim_dt), reps=reps)
            if i_mean_ref is None:
                return False
            i_mean = self.mean_I_analysis(region=[p2delay, p2delay+self.T1], mode='min', baseline=bl,
                intno=i, nint=len(stim_dt), reps=reps)
            # i_mean -= bl.mean()
            # i_mean_ref -= bl.mean()
            cmdv.extend(self.i_mean_cmd)
            stimamp.extend(ptrain['amplitude'])
            stimintvl.append(stim_dt[i])
            self.analysis_summary[f'PPF'][i] = i_mean/i_mean_ref
        self.analysis_summary['psc_stim_amplitudes'] = np.array(stim_I)
        self.analysis_summary['psc_intervals'] = np.array(stimintvl)
        self.analysis_summary['ppf_dt'] = np.array(stim_dt)
        self.analysis_summary['stim_times'] = ptrain['start']
        self.analysis_summary['window'] = [self.T0, self.T1]
        return True
        
        # self.analysis_summmary['psc_amplitudes'] = meani
        # print(self.analysis_summary)
        #self.ivpk_analysis(region=[self.Clamps.tstart, self.Clamps.tstart+0.4*(self.Clamps.tend-self.Clamps.tstart)])

    def set_region(self, region=None, baseline=None, slope=True):
        if region is None:
            raise ValueError("PSPSummary, set_region requires a region beginning and end to measure the current")
    
        data1 = self.Clamps.traces['Time': region[0]:region[1]]
        if baseline is None:
            baseline = [0.]

        tb = np.arange(0, data1.shape[1]*self.Clamps.sample_interval, self.Clamps.sample_interval)
        data1 = data1.view(np.ndarray)
        newCP = CP.CursorPlot(str(self.datapath))
        setline = True

        # for i in range(data1.shape[0]):
        newCP.plotData(x=tb, y=np.array([data1[i]-baseline[i] for i in range(data1.shape[0])])*1e12, setline=setline, slope=slope)
            # setline = False # only do on first plot
    
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
         
        self.T0, self.T1 = newCP.selectedRegion
        if self.T0 == None:
            return(None)
        return(newCP.selectedRegion)
    

    def mean_I_analysis(self, region=None, mode='mean', baseline=None, intno=0, nint=1, reps=[0], slope=True):
        """
        Get the mean current in a window
        
        Parameters
        ----------
        region : tuple, list or numpy array with 2 values (default: None)
            start and end time of a trace used to measure the RMP across
            traces.
        
        Return
        ------
        Nothing
        
        Stores computed mean current in the variable "name".
        """
        if region is None:
            raise ValueError("PSPSummary, mean_I_analysis requires a region beginning and end to measure the current")
        
        data1 = self.Clamps.traces['Time': region[0]:region[1]]

        tb = np.arange(0, data1.shape[1]*self.Clamps.sample_interval, self.Clamps.sample_interval)
        data1 = data1.view(np.ndarray)

        sh = data1.shape
        nreps = len(reps)
        if nint > 1:
            dindx = range(intno, sh[0], nint)
            data1 = data1[dindx]

        # subtract the "baseline" from the beginning of the interval to the end. 
        if slope:
            data1 = self.slope_subtraction(tb, data1, region, mode=mode)
            
        if baseline is not None:
            data1 = np.array([data1[i]-baseline[i] for i in range(data1.shape[0])])
        nx = int(sh[0]/len(reps))
        if mode in ['mean', 'baseline']:
            i_mean = data1.mean(axis=1)  # all traces, average over specified time window
        elif mode == 'min':
            # dfilt = scipy.signal.savgol_filter(data1, 5, 2, axis=1, mode='nearest')
            i_mean = data1.min(axis=1)  # all traces, average over specified time window

        if nint == 1:
            nx = int(sh[0]/len(reps))
            try:
                i_mean = np.reshape(i_mean, (len(reps), nx))  # reshape by repetition
            except:
                return None

        i_mean = i_mean.mean(axis=0)  # average over reps
        self.i_mean_cmd = self.Clamps.commandLevels
        return i_mean

    def slope_subtraction(self, tb, data1, region, mode='mean'):
        """
        Based on """
        dt = tb[1]-tb[0]
        minX = 0 #int((region[0])/dt)
        maxX = int((region[1]-region[0])/dt)
        for i in range(data1.shape[0]):
            x0 = list(range(minX,minX+3))
            ml = maxX
            x0.extend(list(range(ml-10, ml)))
            fdx = np.array([tb[j] for j in x0])
            fdy = np.array([data1[i][j] for j in x0])
            pf = np.polyfit(fdx, fdy, 1)
            bline = np.polyval(pf, tb)
            if bline.shape[0] > data1[i].shape[0]:
                bline = bline[:data1[i].shape[0]]
            if mode != 'baseline':
                data1[i,:] -= bline
        return data1

    def get_traces(self, region=None, trlist=None, baseline=None, intno=0, nint=1, reps=[0], slope=True):
        """
        Get the mean current (averages) in a window
        
        Parameters
        ----------
        region : tuple, list or numpy array with 2 values (default: None)
            start and end time of a trace used to measure the RMP across
            traces.
        
        Return
        ------
        Nothing
        
        Stores computed mean current in the variable "name".
        """
        if region is None:
            raise ValueError("PSPSummary, mean_I_analysis requires a region beginning and end to measure the current")
        
        data1 = self.Clamps.traces['Time': region[0]:region[1]]

        tb = np.arange(0, data1.shape[1]*self.Clamps.sample_interval, self.Clamps.sample_interval)
        data1 = data1.view(np.ndarray)

        sh = data1.shape
        nreps = len(reps)
        if nint > 1:
            dindx = range(intno, sh[0], nint)
            data1 = data1[dindx]

        # subtract the "baseline" from the beginning of the interval to the end. 
        if slope:
            data1 = self.slope_subtraction(tb, data1, region, mode=mode)
            
        if baseline is not None:
            data1 = np.array([data1[i]-baseline[i] for i in range(data1.shape[0])])
        
        nx = int(sh[0]/len(reps))
        # print('unshaped data1: ', data1.shape)
        # prop_cycle = mpl.rcParams['axes.prop_cycle']
        # colors = prop_cycle.by_key()['color']
        #
        # from cycler import cycler
        # from itertools import cycle
        # color_cycle = cycler(c=['r', 'g', 'b', 'c', 'y', 'm', 'k'])
        # f, ax = mpl.subplots(1, 1)
        # k = 0
        # for j in range(len(reps)):
        #     for i, sty in zip(range(nx), cycle(color_cycle)):
        # # for i in range(data1.shape[0]):
        #         ax.plot(tb, data1[k], linewidth=0.5, **sty)
        #         k += 1
        # print('nx: ', nx, '  reps: ', reps, '  len reps: ', len(reps))
        data1 = np.reshape(data1, (len(reps), nx,   -1))
        # print('reshaped data: ', data1.shape)
        data1 = data1.mean(axis=0)
        # print('mean data data: ', data1.shape)
        # for i, sty in zip(range(data1.shape[0]), cycle(color_cycle)):
        #     ax.plot(tb, data1[i], '--', **sty)
        # mpl.show()
        # exit()
        
        # if nint == 1:
        #     nx = int(sh[0]/len(reps))
        #     try:
        #         i_mean = np.reshape(i_mean, (len(reps), nx))  # reshape by repetition
        #     except:
        #         return None

        # i_mean = i_mean.mean(axis=0)  # average over reps
        return data1, tb
                

    def vcss_analysis(self, region=None):
        """
        compute steady-state IV curve - from the mean current 
        across the stimulus set over the defined time region 
        (this usually will be the last half or third of the trace)
        
        Parameters
        ----------
        region : list or tuple
            Start and end times for the analysis
        """
        
        data1 = self.Clamps.traces['Time': region[0]:region[1]]
        icmds = EM.MetaArray(self.Clamps.cmd_wave,  # easiest = turn in to a matching metaarray... 
            info=[{'name': 'Command', 'units': 'A',
             'values': np.array(self.Clamps.clampValues)},
             self.Clamps.traces.infoCopy('Time'), self.Clamps.traces.infoCopy(-1)])
        self.vcss_vcmd = icmds['Time': region[0]:region[1]].mean(axis=1)
        self.r_in = np.nan
        self.analysis_summary['Rin'] = np.nan
        self.vcss_v = []
        if data1.shape[1] == 0 or data1.shape[0] == 1:
            return  # skip it

        ntr = len(self.Clamps.traces)

        self.vcss_Im = data1.mean(axis=1)  # steady-state, all traces
        self.analysis_summary['Rin'] = np.NaN
#        self.Clamps.plotClampData()
        
        isort = np.argsort(self.vcss_vcmd)
        self.vcss_Im= self.vcss_Im[isort]
        self.vcss_vcmd = self.vcss_vcmd[isort]
        bl = self.vcbaseline[isort]
        self.vcss_bl = bl
        # compute Rin from the SS IV:
        # this makes the assumption that:
        # successive trials are in order so we sort above
        # commands are not repeated...
        if len(self.vcss_vcmd) > 1 and len(self.vcss_v) > 1:
            pf = np.polyfit(self.vcss_vcmd, self.vcss_v, 3, rcond=None, full=False, w=None, cov=False)
            pval = np.polyval(pf, self.vcss_vcmd)
            #print('pval: ', pval)
            slope = np.diff(pval) / np.diff(self.vcss_vcmd)  # local slopes
            imids = np.array((self.vcss_vcmd[1:] + self.vcss_vcmd[:-1]) / 2.)
            self.rss_fit ={'I': imids, 'V': np.polyval(pf, imids)}
            #print('fit V: ', self.rss_fit['V'])
            #slope = slope[[slope > 0 ] and [self.vcss_vcmd[:-1] > -0.8] ] # only consider positive slope points
            l = int(len(slope)/2)
            maxloc = np.argmax(slope[l:]) + l
            self.r_in = slope[maxloc]
            self.r_in_loc = [self.vcss_vcmd[maxloc], self.vcss_v[maxloc], maxloc]  # where it was found
            minloc = np.argmin(slope[:l])
            self.r_in_min = slope[minloc]
            self.r_in_minloc = [self.vcss_vcmd[minloc], self.vcss_v[minloc], minloc]  # where it was found
            self.analysis_summary['Rin'] = self.r_in*1.0e-6


    def plot_vciv(self):
        
        P = PH.regular_grid(2 , 2, order='columns', figsize=(8., 6.), showgrid=False,
                        verticalspacing=0.1, horizontalspacing=0.1,
                        margins={'leftmargin': 0.12, 'rightmargin': 0.12, 'topmargin': 0.08, 'bottommargin': 0.1},
                        labelposition=(-0.12, 0.95))
        (date, sliceid, cell, proto, p3) = self.file_cell_protocol(self.datapath)
        
        P.figure_handle.suptitle(os.path.join(date, sliceid, cell, proto).replace('_', '\_'), fontsize=12)
        bl = self.get_baseline()
        if "PPF" in self.analysis_summary.keys():
            maxt = 250.
        else:
            maxt = 150.
        for i in range(self.AR.traces.shape[0]):
            P.axdict['A'].plot(self.AR.time_base*1e3, (self.AR.traces[i,:]-bl[i])*1e12, 'k-', linewidth=0.5)
        if 'PSP_VDEP' in self.analysis_summary.keys():
            P.axdict['A'].set_xlim(self.analysis_summary['stim_times'][0]*1e3-10, self.analysis_summary['stim_times'][0]*1e3+50)
            # P.axdict['A'].set_ylim(-2500, 2500)
        else:
            P.axdict['A'].set_xlim(40, maxt)
            P.axdict['A'].set_ylim(-2000, 500)
        # PH.talbotTicks(P.axdict['A'], tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
        P.axdict['A'].set_xlabel('T (ms)')
        P.axdict['A'].set_ylabel('I (pA)')

        if 'PSP_IO' in self.analysis_summary.keys(): # io function
            for i in range(len(self.analysis_summary['stim_times'])):
                try:
                    P.axdict['C'].plot(self.analysis_summary['psc_stim_amplitudes'],
                        np.array(self.analysis_summary[f'PSP_IO'][i]), linewidth=1, markersize=4)
                except:
                    print('Plot Failed on protocol: ', self.datapath, proto)
            P.axdict['C'].set_xlabel('Istim (microAmps)')
            P.axdict['C'].set_ylabel('EPSC I (pA)')
            PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
        elif 'PSP_VDEP' in self.analysis_summary.keys(): # io function
            for i in range(len(self.analysis_summary['stim_times'])):
                P.axdict['C'].plot(self.analysis_summary['Vcmd']*1e3,
                        self.sign*np.array(self.analysis_summary[f'PSP_VDEP'][i])*1e12, marker='o', linewidth=1, markersize=4)
                P.axdict['C'].plot(self.analysis_summary['Vcmd']*1e3,
                        self.sign*np.array(self.analysis_summary[f'PSP_VDEP_NMDA'][i])*1e12, marker='s', linewidth=1, markersize=4)
            P.axdict['C'].set_xlabel('V (mV)')
            P.axdict['C'].set_ylabel('EPSC I (pA)')
            PH.crossAxes(P.axdict['C'], xyzero=(-60., 0.))
            PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
        elif 'PPF' in self.analysis_summary.keys():
            for i in range(len(self.analysis_summary['stim_times'])):
                P.axdict['C'].plot(self.analysis_summary['ppf_dt']*1e3,
                            self.sign*np.array(self.analysis_summary[f'PPF']), linewidth=1, markersize=4)
                P.axdict['C'].set_xlim(0, 200.)
                P.axdict['C'].set_ylim(0, 2.0)
                PH.referenceline(P.axdict['C'], 1.0)
                P.axdict['C'].set_xlabel('Interval (ms)')
                P.axdict['C'].set_ylabel('PPF (R2/R1)')
                PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 0, 'y': 1}, floatAdd={'x': 0, 'y': 1})


        P.axdict['B'].set_xlabel('I (nA)')
        P.axdict['B'].set_ylabel('V (mV)')
        PH.talbotTicks(P.axdict['B'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 2, 'y': 0})

        P.axdict['D'].set_xlabel('I (pA)')
        P.axdict['D'].set_ylabel('Latency (ms)')

        self.IVFigure = P.figure_handle
    
        if self.plot:
             mpl.show()

    def file_cell_protocol(self, filename):
        """
        file_cell_protocol breaks the current filename down and returns a
        tuple: (date, cell, protocol)
        
        Parameters
        ----------
        filename : str
            Name of the protocol to break down
        
        Returns
        -------
        tuple : (date, sliceid, cell, protocol, any other...)
            last argument returned is the rest of the path...
        """
        (p0, proto) = os.path.split(filename)
        (p1, cell) = os.path.split(p0)
        (p2, sliceid) = os.path.split(p1)
        (p3, date) = os.path.split(p2)
        return (date, sliceid, cell, proto, p3)
        

if __name__ == '__main__':
    
    # disk = '/Volumes/Pegasus/ManisLab_Data3'
    # disk = '/Volumes/PBM_005/data'
    disk = '/Users/pbmanis/Documents/Lab'
    middir = 'Kasten_Michael'
    middir = 'data'
    # middir = ''
    directory = 'Maness_PFC_stim'
    cell = '2019.03.19_000/slice_000/cell_001'
    # cell = '2019.03.19_000/slice_001/cell_000'
    
    ddc = Path(disk, middir, directory, cell)
    # protocol = 'Stim_IO_1_001'
    protocol = 'PPF_2_001'
    # protocol = 'VC-EPSC_3_ver2_003'
    fn = Path(ddc, protocol)
    PSC = PSCSummary(fn)
    PSC.compute_PSC_IV(protocol[:-4], savetimes=True)
    
