"""
Analyze Rm, tau - pulled out of IVCurve 2/29/2016 pbm.
Allows routine to be used to responses to hyperpolarizing pulses independent of acq4's data models.
Create instance, then  call setup to define the "Clamps" structure and analysis parameters. 

Clamps must have the following variables defined:

    commandLevels (current injection levels, list)
    time_base (np.array of times corresponding to traces)
    data_mode ()
    tstart (time for start of looking at spikes; ms)
    tend
    trace
    sample_interval (time between samples, sec)
    values (command waveforms; why it is called this in acq4 is a mystery)

The "Clamps" object can be provided by acq4's PatchEphys module, or by
an instance of acq4read.

RmTauAnalysis requires that the SpikeAnalysis be run first.

Paul B. Manis, 2015-2017
for acq4

"""


from collections import OrderedDict
import os
import os.path
import itertools
import functools
import numpy as np
import scipy
import pyqtgraph as pg
import Fitting  # provided local copy
import Utility
import pprint
import time

class RmTauAnalysis():
    """
    RmTau analysis measures Rm (input resistance) and the membrane time constant for
    the traces in the Clamps structure.
    It also makes some measures related to Ih for hyperpolarizing pulses, including
    the peak and steady-state voltages, and the time constant of Ih activation (e.g.
    the sag itself after the peak).
    
    To use: create an instance of the class, then call RmTauAnalysis.setup(...)
    to initialized the parameters. You can then call any of the methods to get
    the tau, rmp, etc. from the data set
    """
    
    def __init__(self):
        self.Clamps = None
        self.Spikes = None
        self.dataPlot = None
        self.baseline = None
        self.rmp = []
        self.taum_fitted = {}
        self.taum_bounds = []
        self.analysis_summary = {}

        

    def setup(self, clamps=None, spikes=None, dataplot=None, baseline=[0, 0.001], taumbounds = [0.002, 0.050]):
        """
        Set up for the fitting
        
        Parameters
        ----------
        clamps : A datamodel structure (required)
            Brings the data to the module. This usually will be a PatchEphys object.
        
        spikes : A spikeAnalysis structure (required)
            Has information about which traces have spikes
        
        dataplot : pyqtgraph plot object
            pyqtgraph plot to use for plotting the data. No data is plotted if not specified
        
        baseline : list (2 elements)
            times over which baseline is measured (in seconds)
        
        taumbounds : list (2 elements)
            Lower and upper bounds of the allowable taum fitting range (in seconds).
        """
        
        if clamps is None or spikes is None:
            raise ValueError("RmTau analysis requires defined clamps and spike analysis")
        self.Clamps = clamps
        self.Spikes = spikes
        self.dataPlot = dataplot
        self.baseline = baseline
#        self.taum_fitted = {}
        self.tauh_fitted = {}
        self.taum_bounds = taumbounds
    
    def analyze(self, rmpregion=[0., 0.05], tauregion=[0.1, 0.125]):
        self.rmp_analysis(region=rmpregion)
        self.tau_membrane(region=tauregion)
        r0 = self.Clamps.tstart + 0.9*(self.Clamps.tend-self.Clamps.tstart) # 
        self.ivss_analysis(region=[r0, self.Clamps.tend])
        
           
    def tau_membrane(self, peak_time=None, printWindow=False, whichTau=1, vrange=[-0.005, -0.020], region=[]):
        """
        Compute time constant (single exponential) from the onset of the response to a current step
        
        Parameters
        ----------
        
        peak_time : float (ms) (default: None)
            Time to the peak of the data for the fitting. If it is set (not None), the 
            fit will only go until the peak time. Otherwise it is set by the end of the region.
        
        printWindow : Boolean (default: False)
            Flag to allow printing of the fit results in detail during a run
        
        whichTau : int (default: 1)
            Not used
        
        vrange : list (V) (default: [-0.005, -0.020])
            Define the voltage range below RMP for the traces that will be fit to obtain tau_m.
            
        
        region: list (s) (default: [])
            Define the time region for the fitting
        
        Return
        ------
            Nothing
            
        Class variables with a leading taum_ are set by this routine, to return results.
        
        """
        
        assert len(region) > 0.
        
        rgnpk = list(region)
        Func = 'exp1'  # single exponential fit with DC offset.
        if self.rmp == []:
            self.rmp_analysis(region=self.baseline)

        Fits = Fitting()  # get a fitting instance
        initpars = [self.rmp*1e-3, 0.010, 0.01]  # rmp is in units of mV
        icmdneg = np.where(self.Clamps.commandLevels < -20e-12)
        maxcmd = np.min(self.Clamps.commandLevels)
        ineg = np.where(self.Clamps.commandLevels[icmdneg] < 0.0)
        if peak_time is not None and ineg != np.array([]):
            rgnpk[1] = np.max(peak_time[ineg[0]])
        dt = self.Clamps.sample_interval
        rgnpk = sorted(rgnpk)
        vrange = np.sort(vrange) # /1000.
        rgnindx = [int((rgnpk[1]-0.005)/dt), int((rgnpk[1])/dt)]
        rmps = self.ivbaseline
        vmeans = np.mean(self.Clamps.traces[:, rgnindx[0]:rgnindx[1]], axis=1) - self.ivbaseline
        # vrange is in mV
        indxs = np.where(np.logical_and((vmeans[ineg] >= vrange[0]), 
                         (vmeans[ineg] <= vrange[1])))
        # print 'ineg: ', ineg
        # print self.Clamps.commandLevels
        # print 'icmdneg: ', icmdneg
        # print 'vrange: ', vrange
        # print 'vmeans[ineg]: ', vmeans[ineg]
        # print 'indxs: ', indxs
        indxs = list(indxs[0])
        whichdata = ineg[0][indxs]  # restricts to valid values
        itaucmd = self.Clamps.commandLevels[ineg]
        whichaxis = 0
        fpar = []
        names = []
        okdata = []
        if len(self.taum_fitted.keys()) > 0:
            [self.taum_fitted[k].clear() for k in self.taum_fitted.keys()]
        self.taum_fitted = {}

        for j, k in enumerate(whichdata):
            if self.dataPlot is not None:
                n = 'w_%03d' % j
                self.taum_fitted[n] = self.dataPlot.plot(self.Clamps.time_base,
                                     self.Clamps.traces[k], pen=pg.mkPen('w'))
            (fparx, xf, yf, namesx) = Fits.FitRegion([k], whichaxis,
                                               self.Clamps.time_base,
                                               self.Clamps.traces,
                                               dataType='2d',
                                               t0=rgnpk[0], t1=rgnpk[1],
                                               fitFunc=Func,
                                               fitPars=initpars,
                                               method='SLSQP',
                                               bounds=[(-0.1, 0.1), (-0.1, 0.1), 
                                               (self.taum_bounds[0], self.taum_bounds[1])])
        
            if not fparx:
              raise Exception('IVCurve::update_Tau_membrane: Charging tau fitting failed - see log')
            #print 'j: ', j, len(fpar)
            if fparx[0][1] < 2.5e-3:  # amplitude must be > 2.5 mV to be useful
                continue
            fpar.append(fparx[0])
            names.append(namesx[0])
            okdata.append(k)
        self.taum_pars = fpar
        self.taum_win = rgnpk
        self.taum_func = Func
        self.taum_whichdata = okdata
        taus = []
        for j in range(len(fpar)):
            outstr = ""
            taus.append(fpar[j][2])
            for i in range(0, len(names[j])):
                outstr += '%s = %f, ' % (names[j][i], fpar[j][i])
            if printWindow:
                print("FIT(%d, %.1f pA): %s " %
                      (whichdata[j], itaucmd[j] * 1e12, outstr))
        if len(taus) > 0:
            self.taum_taum = np.nanmean(taus)
            self.analysis_summary['taum'] = self.taum_taum
        else:
            self.taum_taum = np.NaN
            self.analysis_summary['taum'] = np.NaN
        self.analysis_summary['taupars'] = self.taum_pars
        self.analysis_summary['taufunc'] = self.taum_func
    
    def rmp_analysis(self, region=None):
        """
        Get the resting membrane potential
        
        Parameters
        ----------
        region : tuple, list or numpy array with 2 values (default: None)
            start and end time of a trace used to measure the RMP across
            traces.
        
        Return
        ------
        Nothing
        
        Stores computed RMP in mV in the class variable rmp.
        """
        if region is None:
            raise ValueError("IVCurve, rmp_analysis requires a region beginning and end to measure the RMP")
        data1 = self.Clamps.traces['Time': region[0]:region[1]]
        data1 = data1.view(np.ndarray)
        self.ivbaseline = data1.mean(axis=1)  # all traces
        self.ivbaseline_cmd = self.Clamps.commandLevels
        self.rmp = np.mean(self.ivbaseline) * 1e3  # convert to mV
        self.analysis_summary['RMP'] = self.rmp

    def ivss_analysis(self, region=None):
        data1 = self.Clamps.traces['Time': region[0]:region[1]]
 #       print 'data shape: ', data1.shape
        if data1.shape[1] == 0 or data1.shape[0] == 1:
            return  # skip it
        self.ivss_v = []

        # check out whether there are spikes in the window that is selected
        threshold = self.Spikes
        ntr = len(self.Clamps.traces)
        if not self.Spikes.spikes_counted:
            print "ivss_analysis: spikes not counted yet? - let's go analyze them..."
            self.analyzeSpikes()

        self.ivss_v = data1.mean(axis=1)  # all traces
        # if self.ctrl.IVCurve_SubBaseline.isChecked():
        #     self.ivss = self.ivss - self.RmTau.ivbaseline
        self.analysis_summary['Rin'] = np.NaN
        if len(self.Spikes.nospk) >= 1:
            # Steady-state IV where there are no spikes
            self.ivss_v = self.ivss_v[self.Spikes.nospk]
            self.ivss_cmd = self.Clamps.commandLevels[self.Spikes.nospk]
#            self.commandLevels = commands[self.nospk]
            # compute Rin from the SS IV:
            # this makes the assumption that:
            # successive trials are in order (as are commands)
            # commands are not repeated...
            if len(self.ivss_cmd) > 1 and len(self.ivss_v) > 1:
                self.r_in = np.max(np.diff
                                   (self.ivss_v) / np.diff(self.ivss_cmd))
                self.analysis_summary['Rin'] = self.r_in*1.0e-6


        isort = np.argsort(self.ivss_cmd)
        self.ivss_cmd = self.ivss_cmd[isort]
        self.ivss_v = self.ivss_v[isort]


    def leak_subtract(self):
        self.yleak = np.zeros(len(self.ivss_v))
        # if self.ctrl.IVCurve_subLeak.isChecked():
        #     if self.Clamps.data_mode in self.dataModel.ic_modes:
        #         sf = 1e-12
        #     elif self.Clamps.data_mode in self.dataModel.vc_modes:
        #         sf = 1e-3
        #     else:
        #         sf = 1.0
        #     (x, y) = Utility.clipdata(self.ivss, self.ivss_cmd,
        #                               self.ctrl.IVCurve_LeakMin.value() * sf,
        #                               self.ctrl.IVCurve_LeakMax.value() * sf)
        #     try:
        #         p = np.polyfit(x, y, 1)  # linear fit
        #         self.yleak = np.polyval(p, self.ivss_cmd)
        #         self.ivss = self.ivss - self.yleak
        #     except:
        #         raise ValueError('IVCurve Leak subtraction: no valid points to correct')
        
    
    def tau_h(self, current, rgn, pkRgn, ssRgn, printWindow=False):
        """
        Measure the time constant associated with activation of the hyperpolarization-
        activated current, Ih. The tau is measured from the peak of the response to the
        steady-state part, at a single current level.
        
        Parameters
        ----------
        current : float (current, pA; no default).
             The current level that will be used as the target for measuring Ih. A single
             current level is given; the closest one in the test set will be used.
        
        rgn : list of floats ([time, time], ms, no default)
            time window over which data will be fit for the tau_h measure
            
        pkRgn : list of floats ([time, time], ms; no default)
            The time window over which the peak voltage will be identified.
        
        ssRgn : list of floats ([time, time], ms; no default)
            The time window over which the steady-state voltage will be identified.
        
        Return
        ------
        Nothing
            
        Class variables with a leading tauh_ are set by this routine, to return the
        results of the measurements.
        
        """
        Func = 'exp1'  # single exponential fit to the whole region
        Fits = Fitting.Fitting()

        initpars = [-80.0 * 1e-3, -10.0 * 1e-3, 50.0 * 1e-3]

        # find the current level that is closest to the target current
        itarget = self.Clamps.values[current]  # retrive actual value from commands
        self.tauh_neg_cmd = itarget
        idiff = np.abs(np.array(self.Clamps.commandLevels) - itarget)
        amin = np.argmin(idiff)  # amin appears to be the same as s_target
        # target trace (as selected in cmd drop-down list):
        target = self.Clamps.traces[amin]
        # get Vrmp -  # rmp approximation.
        vrmp = np.median(target['Time': 0.0:self.Clamps.tstart - 0.005]) * 1000.
        self.tauh_vrmp = vrmp
        # get peak and steady-state voltages
        vpk = target['Time': pkRgn[0]:pkRgn[1]].min() * 1000
        self.tauh_vpk = vpk
        self.tauh_neg_pk = (vpk - vrmp) / 1000.
        vss = np.median(target['Time': ssRgn[0]:ssRgn[1]]) * 1000
        self.tauh_vss = vss
        self.tauh_neg_ss = (vss - vrmp) / 1000.
        whichdata = [int(amin)]
        itaucmd = [self.Clamps.commandLevels[amin]]
        fd = self.Clamps.traces['Time': rgn[0]:rgn[1]][whichdata][0]
        if len(self.tauh_fitted.keys()) > 0:
            [self.tauh_fitted[k].clear() for k in self.tauh_fitted.keys()]
        whichaxis = 0
        (fpar, xf, yf, names) = Fits.FitRegion(whichdata, whichaxis,
                                               self.Clamps.traces.xvals('Time'),
                                               self.Clamps.traces.view(np.ndarray),
                                               dataType='2d',
                                               t0=rgn[0], t1=rgn[1],
                                               fitFunc=Func,
                                               fitPars=initpars)
        if not fpar:
            raise Exception('IVCurve::update_Tauh: tau_h fitting failed - see log')
        s = np.shape(fpar)
        taus = []
        for j in range(0, s[0]):
            outstr = ""
            taus.append(fpar[j][2])
            for i in range(0, len(names[j])):
                outstr += '%s = %f, ' % (names[j][i], fpar[j][i])
            if printWindow:
                print("Ih FIT(%d, %.1f pA): %s " %
                      (whichdata[j], itaucmd[j] * 1e12, outstr))
        self.tauh_xf = xf
        self.tauh_yf = yf
        self.tauh_meantau = np.mean(taus)
        self.tauh_bovera = (vss - vrmp) / (vpk - vrmp)
        Gpk = itarget / self.tauh_neg_pk
        Gss = itarget / self.tauh_neg_ss
        self.tauh_Gh = Gss - Gpk
        
        
