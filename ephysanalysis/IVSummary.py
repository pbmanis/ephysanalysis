from __future__ import print_function

"""
Compute IV

"""

import sys

import matplotlib
import numpy as np

import os.path
import ephysanalysis as EP
import matplotlib.pyplot as mpl
import matplotlib.colors
import matplotlib
#import colormaps.parula
import pylibrary.PlotHelpers as PH
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
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
    def __init__(self, datapath):
        self.datapath = datapath
        self.AR = EP.acq4read.Acq4Read()  # make our own private cersion of the analysis and reader
        self.SP = EP.SpikeAnalysis.SpikeAnalysis()
        self.RM = EP.RmTauAnalysis.RmTauAnalysis()
    
    def compute_iv(self):
        """
        Simple plot of spikes, FI and subthreshold IV
        """
        print('path: ', self.datapath)
        self.AR.setProtocol(self.datapath)  # define the protocol path where the data is
        self.AR.getData()  # get that data.
        self.SP.setup(clamps=self.AR, threshold=-0.010,  refractory=0.0001, peakwidth=0.001, interpolate=False, verify=False, mode='peak')
        self.SP.analyzeSpikes()
        self.SP.analyzeSpikeShape()
        self.RM.setup(self.AR, self.SP)
        self.RM.analyze()
        self.plot_iv()

    def plot_iv(self):
        P = PH.regular_grid(2 , 2, order='columns', figsize=(8., 6.), showgrid=False,
                        verticalspacing=0.1, horizontalspacing=0.12,
                        margins={'leftmargin': 0.12, 'rightmargin': 0.12, 'topmargin': 0.08, 'bottommargin': 0.1},
                        labelposition=(-0.12, 0.95))
        P.figure_handle.suptitle(self.datapath, fontsize=6)
        for i in range(self.AR.traces.shape[0]):
            P.axdict['A'].plot(self.AR.time_base*1e3, self.AR.traces[i,:]*1e3, 'k-', linewidth=0.5)
        P.axdict['B'].plot(self.SP.analysis_summary['FI_Curve'][0]*1e9, self.SP.analysis_summary['FI_Curve'][1]/(self.AR.tend-self.AR.tstart), 'ko-', markersize=4, linewidth=0.5)
        P.axdict['C'].plot(self.RM.ivss_cmd*1e9, self.RM.ivss_v*1e3, 'ko-', markersize=4, linewidth=1.0)
        P.axdict['C'].text(-0.05, 0.95, r'RMP: {0:.1f} mV {1:s}${{R_{{in}}}}$: {2:.1f} ${{M\Omega}}${3:s}${{\tau_{{m}}}}$: {4:.2f} ms'.format(
                    self.RM.analysis_summary['RMP'], '\n', self.RM.analysis_summary['Rin'], '\n', self.RM.analysis_summary['taum']*1e3),
                    transform=P.axdict['C'].transAxes, horizontalalignment='left', verticalalignment='top')
     #   P.axdict['C'].xyzero=([0., -0.060])
        PH.talbotTicks(P.axdict['A'], tickPlacesAdd={'x': 0, 'y': 0}, floatAdd={'x': 0, 'y': 0})
        P.axdict['A'].set_xlabel('T (ms)')
        P.axdict['A'].set_ylabel('V (mV)')
        P.axdict['B'].set_xlabel('I (nA)')
        P.axdict['B'].set_ylabel('Spikes/s')
        PH.talbotTicks(P.axdict['B'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 2, 'y': 0})
        PH.crossAxes(P.axdict['C'], xyzero=(0., -60))
        PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 2, 'y': 0})
        P.axdict['C'].set_xlabel('I (nA)')
        P.axdict['C'].set_ylabel('V (mV)')

        flat = []
        fisi = []
        iflatcur = []
        iisicur = []
        for k in self.SP.analysis_summary['spikes'].keys():
            flat.append(self.SP.analysis_summary['spikes'][k][0]['AP_Latency'])
            iflatcur.append(self.SP.analysis_summary['spikes'][k][0]['current'])
            try:
                fisi.append(self.SP.analysis_summary['spikes'][k][1]['AP_Latency']-self.SP.analysis_summary['spikes'][k][0]['AP_Latency'])
                iisicur.append(self.SP.analysis_summary['spikes'][k][0]['current'])
            except:
                pass
        P.axdict['D'].plot(np.array(iflatcur)*1e9, (np.array(flat)-self.AR.tstart)*1000., 'sk-')
        P.axdict['D'].plot(np.array(iisicur)*1e9, (np.array(fisi))*1000., '^b-')
        PH.talbotTicks(P.axdict['C'], tickPlacesAdd={'x': 1, 'y': 0}, floatAdd={'x': 1, 'y': 0})
        P.axdict['D'].set_xlabel('I (nA)')
        P.axdict['D'].set_ylabel('Latency (ms)')
        
        mpl.show()

