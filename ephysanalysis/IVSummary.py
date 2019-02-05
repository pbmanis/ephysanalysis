from __future__ import print_function

"""
Compute IV

"""

import sys

import matplotlib
import numpy as np

import os.path
from collections import OrderedDict
import ephysanalysis as EP
import matplotlib
#matplotlib.use('MacOSX')
import matplotlib.pyplot as mpl
import matplotlib.colors
import matplotlib
#import colormaps.parula
import pylibrary.PlotHelpers as PH
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#rcParams['font.sans-serif'] = ['Arial']
#rcParams['font.family'] = 'sans-serif'
rc('text', usetex=True)
import matplotlib
rcParams = matplotlib.rcParams
rcParams['svg.fonttype'] = 'none' # No text as paths. Assume font installed.
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['text.latex.unicode'] = True
#import seaborn
#cm_sns = mpl.cm.get_cmap('terrain')  # terrain is not bad
#cm_sns = mpl.cm.get_cmap('parula')  # terrain is not bad
#cm_sns = mpl.cm.get_cmap('jet')  # jet is terrible
color_sequence = ['k', 'r', 'b']
colormap = 'snshelix'

      
class IVSummary():
    def __init__(self, datapath, plot=True):
        self.datapath = datapath

        self.IVFigure = None
        self.AR = EP.acq4read.Acq4Read()  # make our own private version of the analysis and reader
        self.SP = EP.SpikeAnalysis.SpikeAnalysis()
        self.RM = EP.RmTauAnalysis.RmTauAnalysis()
        self.plot = plot

    def compute_iv(self, threshold=-0.010, bridge_offset=0.0, tgap=0.0005, pubmode=False, plotiv=True):
        """
        Simple plot of spikes, FI and subthreshold IV
        
        """
        #print('path: ', self.datapath)
        self.AR.setProtocol(self.datapath)  # define the protocol path where the data is
        if self.AR.getData():  # get that data.
            # print(self.AR.tstart, self.AR.tend)
            self.RM.setup(self.AR, self.SP, bridge_offset=bridge_offset)
            self.SP.setup(clamps=self.AR, threshold=threshold, 
                    refractory=0.0001, peakwidth=0.001, interpolate=True, verify=False, mode='peak')
            self.SP.analyzeSpikes()
            self.SP.analyzeSpikeShape()
            self.SP.analyzeSpikes_brief(mode='baseline')
            self.SP.analyzeSpikes_brief(mode='poststimulus')
            self.SP.fitOne(function='fitOneOriginal')
           # self.SP.fitOne(function='piecewiselinear3')
            self.RM.analyze(rmpregion=[0., self.AR.tstart-0.001],
                            tauregion=[self.AR.tstart,
                                       self.AR.tstart + (self.AR.tend-self.AR.tstart)/5.],
                            to_peak=True, tgap=tgap)
            if plotiv:
                self.plot_iv(pubmode=pubmode)
            return True
        else:
            print('getData did not return data')
            return False

    def plot_iv(self, pubmode=False):
        x = -0.08
        y = 1.02
        sizer = {'A': {'pos': [0.05, 0.50, 0.08, 0.78], 'labelpos': (x,y), 'noaxes': False},
                 'B': {'pos': [0.62, 0.30, 0.64, 0.22], 'labelpos': (x,y), 'noaxes': False},
                 'C': {'pos': [0.62, 0.30, 0.34, 0.22], 'labelpos': (x,y)},
                 'D': {'pos': [0.62, 0.30, 0.08, 0.22], 'labelpos': (x,y)}, 
                }

        # dict pos elements are [left, width, bottom, height] for the axes in the plot.
        gr = [(a, a+1, 0, 1) for a in range(0, 4)]   # just generate subplots - shape does not matter
        axmap = OrderedDict(zip(sizer.keys(), gr))
        P = PH.Plotter((4, 1), axmap=axmap, label=True, figsize=(8., 6.))
        # PH.show_figure_grid(P.figure_handle)
        P.resize(sizer)  # perform positioning magic


        # P = PH.regular_grid(2 , 2, order='columns', figsize=(8., 6.), showgrid=False,
        #                 verticalspacing=0.1, horizontalspacing=0.12,
        #                 margins={'leftmargin': 0.12, 'rightmargin': 0.12, 'topmargin': 0.08, 'bottommargin': 0.1},
        #                 labelposition=(-0.12, 0.95))
        P.figure_handle.suptitle(self.datapath, fontsize=8)
        dv = 50.
        jsp = 0
        for i in range(self.AR.traces.shape[0]):
            if i in list(self.SP.spikeShape.keys()):
                idv = float(jsp)*dv
                jsp += 1
            else:
                idv = 0.
            P.axdict['A'].plot(self.AR.time_base*1e3, idv + self.AR.traces[i,:].view(np.ndarray)*1e3, '-', linewidth=0.35)
            ptps = np.array([])
            paps = np.array([])
            if i in list(self.SP.spikeShape.keys()):
                for j in list(self.SP.spikeShape[i].keys()):
                    paps = np.append(paps, self.SP.spikeShape[i][j]['peak_V']*1e3)
                    ptps = np.append(ptps, self.SP.spikeShape[i][j]['peak_T']*1e3)
                P.axdict['A'].plot(ptps, idv+paps, 'ro', markersize=0.5)
            
            # mark spikes outside the stimlulus window
            ptps = np.array([])
            paps = np.array([])
            for window in ['baseline', 'poststimulus']:
                ptps = np.array(self.SP.analysis_summary[window+'_spikes'][i])
                uindx = [int(u/self.AR.sample_interval)+1 for u in ptps]
                paps = np.array(self.AR.traces[i, uindx])
                P.axdict['A'].plot(ptps*1e3, idv+paps*1e3, 'bo', markersize=0.5)
                    
        for k in self.RM.taum_fitted.keys():
            P.axdict['A'].plot(self.RM.taum_fitted[k][0]*1e3, self.RM.taum_fitted[k][1]*1e3, '--k', linewidth=0.30)
        for k in self.RM.tauh_fitted.keys():
            P.axdict['A'].plot(self.RM.tauh_fitted[k][0]*1e3, self.RM.tauh_fitted[k][1]*1e3, '--r', linewidth=0.50)
        if pubmode:
            PH.calbar(P.axdict['A'], calbar=[0., -90., 25., 25.], axesoff=True, 
                orient='left', unitNames={'x': 'ms', 'y': 'mV'}, fontsize=10, weight='normal', font='Arial')
        P.axdict['B'].plot(self.SP.analysis_summary['FI_Curve'][0]*1e9, self.SP.analysis_summary['FI_Curve'][1]/(self.AR.tend-self.AR.tstart), 'ko-', markersize=4, linewidth=0.5)
        clist = ['r', 'b', 'g', 'c', 'm']  # only 5 possiblities
        linestyle = ['-', '--', '-.', '-', '--']
        n = 0
        for i, figrowth in enumerate(self.SP.analysis_summary['FI_Growth']):
            legstr = '{0:s}\n'.format(figrowth['FunctionName'])
            if len(figrowth['parameters']) == 0: # no valid fit
                P.axdict['B'].plot([np.nan, np.nan], [np.nan, np.nan], label='No valid fit')
            else:
                for j, fna in enumerate(figrowth['names'][0]):
                    legstr += '{0:s}: {1:.3f} '.format(fna, figrowth['parameters'][0][j])
                    if j in [2, 5, 8]:
                        legstr += '\n' 
                P.axdict['B'].plot(figrowth['fit'][0][0]*1e9, figrowth['fit'][1][0],
                    linestyle=linestyle[i], color=clist[i], linewidth=0.5,
                    label=legstr)
            n += 1
        if n > 0:
            P.axdict['B'].legend(fontsize=6)
        
        P.axdict['C'].plot(self.RM.ivss_cmd*1e9, self.RM.ivss_v*1e3, 'ko-', markersize=4, linewidth=1.0)
        if not pubmode:
            if self.RM.analysis_summary['CCComp']['CCBridgeEnable'] == 1:
                enable = 'On'
            else:
                enable = 'Off'
            tstr = (r'RMP: {0:.1f} mV {1:s}${{R_{{in}}}}$: {2:.1f} ${{M\Omega}}${3:s}${{\tau_{{m}}}}$: {4:.2f} ms'.
                    format(self.RM.analysis_summary['RMP'], '\n', 
                           self.RM.analysis_summary['Rin'], '\n', 
                           self.RM.analysis_summary['taum']*1e3))
            tstr += (r'{0:s}Holding: {1:.1f} pA{2:s}Bridge [{3:3s}]: {4:.1f} ${{M\Omega}}$'.
                     format('\n', np.mean(self.RM.analysis_summary['Irmp'])*1e12,
                            '\n', enable, 
                            np.mean(self.RM.analysis_summary['CCComp']['CCBridgeResistance']/1e6)))
            tstr += (r'{0:s}Bridge Adjust: {1:.1f} ${{M\Omega}}$ {2:s}Pipette: {3:.1f} mV'.
                    format('\n', self.RM.analysis_summary['BridgeAdjust']/1e6,
                           '\n', np.mean(self.RM.analysis_summary['CCComp']['CCPipetteOffset']*1e3)))

            P.axdict['C'].text(-0.05, 0.80, tstr,
                transform=P.axdict['C'].transAxes, horizontalalignment='left', verticalalignment='top', fontsize=7)
     #   P.axdict['C'].xyzero=([0., -0.060])
        PH.talbotTicks(P.axdict['A'], tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
        P.axdict['A'].set_xlabel('T (ms)')
        P.axdict['A'].set_ylabel('V (mV)')
        P.axdict['B'].set_xlabel('I (nA)')
        P.axdict['B'].set_ylabel('Spikes/s')
        PH.talbotTicks(P.axdict['B'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 2, 'y': 0})
        try:
            maxv = np.max(self.RM.ivss_v*1e3)
        except:
            maxv = 0.  # sometimes IVs do not have negative voltages for an IVss to be available... 
        ycross = np.around(maxv/5., decimals=0)*5.
        if ycross > maxv:
            ycross = maxv
        PH.crossAxes(P.axdict['C'], xyzero=(0., ycross))
        PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 2, 'y': 0})
        P.axdict['C'].set_xlabel('I (nA)')
        P.axdict['C'].set_ylabel('V (mV)')

        for i in range(len(self.SP.spikes)):
            if len(self.SP.spikes[i]) == 0:
                continue
            spx = np.argwhere((self.SP.spikes[i] > self.SP.Clamps.tstart) & (self.SP.spikes[i] <= self.SP.Clamps.tend)).ravel()
            spkl = (np.array(self.SP.spikes[i][spx])-self.SP.Clamps.tstart )*1e3 # just shorten...
            if len(spkl) == 1:
                P.axdict['D'].plot(spkl[0], spkl[0], 'or', markersize=4)
            else:
                P.axdict['D'].plot(spkl[:-1], np.diff(spkl), 'o-', markersize=3, linewidth=0.5)
                
        PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 1, 'y': 0})
        P.axdict['D'].set_yscale('log')
        P.axdict['D'].set_ylim((1.0, P.axdict['D'].get_ylim()[1]))
        P.axdict['D'].set_xlabel('Latency (ms)')
        P.axdict['D'].set_ylabel('ISI (ms)')
        P.axdict['D'].text(0.05, 0.05, 'Adapt Ratio: {0:.3f}'.format(self.SP.analysis_summary['AdaptRatio']), fontsize=9,
            transform=P.axdict['D'].transAxes, horizontalalignment='left', verticalalignment='bottom')
        self.IVFigure = P.figure_handle
    
        if self.plot:
             mpl.show()

