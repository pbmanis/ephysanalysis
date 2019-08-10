from __future__ import print_function
#!/usr/bin/python

"""
Simple adjunct routine to plot LSPS/CRACM maps with traces, over cell image if possible.
Reuqires acq4read. 

Takes advantage of acq4read code to access all data and attributes.


"""
import os
import re
import itertools
from collections import OrderedDict
from pathlib import Path
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as mpl
from pylibrary import PlotHelpers as PH
import seaborn as sns
import ephysanalysis.acq4read as ARC
import ephysanalysis.metaarray as EM
from pyqtgraph import configfile
from pylibrary import picker
import scipy.ndimage
import numpy as np
import datetime
import pprint
import textwrap as WR
import collections
import tifffile as tf
import boundrect as BR
import mapanalysistools.digital_filters as FILT
import mahotas as MH

class ScannerInfo(object):
    """
    Get scanner information, compute the scanner box and some additional parameters
    Do this as a class to encapsulate the data and make reusable.
    """
    def __init__(self, AR):
        BRI = BR.BoundRect()
        self.AR = AR  # save the acq4read instance for access to the data
        self.AR.getScannerPositions()
        self.scannerpositions = np.array(AR.scannerpositions)
        pos = self.AR.scannerCamera['frames.ma']['transform']['pos']
        scale = self.AR.scannerCamera['frames.ma']['transform']['scale']
        region = self.AR.scannerCamera['frames.ma']['region']
        self.binning = self.AR.scannerCamera['frames.ma']['binning']
        scale = list(scale)
        scale[0] = scale[0]/self.binning[0]
        scale[1] = scale[1]/self.binning[1]
        self.scale = scale
        if self.AR.spotsize is not None:
            print ('Spot Size: {0:0.3f} microns'.format(self.AR.spotsize*1e6))
        else:
            self.AR.spotsize=50.

        self.camerabox = [[pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[1]],
               [pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[3]],
               [pos[0] + scale[0]*region[2], pos[1] + scale[1]*region[3]],
               [pos[0] + scale[0]*region[2], pos[1] + scale[1]*region[1]],
               [pos[0] + scale[0]*region[0], pos[1] + scale[1]*region[1]]
           ]
        scannerbox = BRI.getRectangle(self.AR.scannerpositions)
        if scannerbox is None:  # likely just one point
            pt = self.AR.scannerpositions
            fp = np.array([[pt[0][0]], [pt[0][1]]])
            scannerbox = fp
        else:
            fp = np.array([scannerbox[0][0], scannerbox[1][1]]).reshape(2,1)
        print('fp: ', fp)
        scannerbox = np.append(scannerbox, fp, axis=1)
        self.scboxw = np.array(scannerbox)        
        self.boxw = np.swapaxes(np.array(self.camerabox), 0, 1)
        

class MapTraces(object):
    def __init__(self):
        self.cell = None
        self.datasets = OrderedDict()
        self.image = None
        self.AR = ARC.Acq4Read()
        self.invert = True
        self.vmax = 20000.
        self.voff = 0.
        self.ioff = 0.050
        self.basezero = True
        self.xscale = 1.0
        self.yscale = 1.0
        self.nspots = 0
        self.overlay = True
        self.indicesplotted = []
        self.twin = [0, 0.6]
        self.averageScannerImages = False # normally, would not do
        self.calbar = [0.1, 5e-7]  # 100 ms, 500 pA
        self.picker = picker.Picker()
        sns.set()
        sns.color_palette("colorblind", 10)
        self.palette = itertools.cycle(sns.color_palette("colorblind", 10))
        sns.set_style("white")
        sns.set_style("ticks")        

    def setProtocol(self, cell, image=None, videos=None):
        self.cell = Path(cell)
        if not self.cell.is_dir():
            print(f"did not find directory: {str(cell):s}")
            raise ValueError
        if image is not None:
            self.image = Path(self.cell, image)
            print(self.image)
            print(self.cell)
            self.image_data = self.AR.getImage(self.image)
        else:
            self.image = None
        self.videos = []
        if videos is not None:
            for v in videos:
                self.videos.append(Path(self.cell, f"video_0{v:02d}"))
        self.AR.setProtocol(self.cell)
        
    
    def setPars(self, pdict):
        for k in list(pdict.keys()):
            if k == 'invert':
                self.invert = pdict[k]
            if k == 'vmax':
                self.vmax = pdict[k]
            if k == 'voff':
                self.voff = pdict[k]
            if k == 'ioff':
                self.ioff = pdict[k]
            if k == 'xscale':
                self.xscale = pdict[k]
            if k == 'yscale':
                self.yscale = pdict[k]
            if k == 'calbar':
                self.calbar[0] = pdict[k][0]
                self.calbar[1] = pdict[k][1]
            if k == 'twin':
                self.twin[0] = pdict[k][0]
                self.twin[1] = pdict[k][1]

    def plot_maps(self, protocols, traces=None):
        """
        Plot map or superimposed maps...
        """
        self.figure = mpl.figure()
        # print(dir(self.figure))
        self.figure.set_size_inches(14., 8.)
        self.ax = self.figure.add_subplot('131')
        self.ax2 = self.figure.add_subplot('132')
        sns.despine(ax=self.ax2, left=True, bottom=True, right=True, top=True)
        # self.ax3 = self.figure.add_subplot('133')
        self.data = dict.fromkeys(list(protocols.keys()))
        cols = ['r', 'b', 'c', 'g']
        self.traces = traces

        for i, p in enumerate(protocols):
            prot = protocols[p]
            self.datasets[p] = []
            if i == 0:
                self.setProtocol(prot, image)
                self.show_traces(self.figure, self.ax, pcolor=cols[i], name=p)
            else:
                self.setProtocol(prot)
                self.show_traces(self.figure, self.ax, pcolor=cols[i], name=p)
        if self.traces is not None:
            for tr in self.traces:
                self.handle_event(index=tr)
            PH.calbar(self.ax2, calbar=[0, -6., 50, 5.], scale=[1.0, 1.0],
                axesoff=True, orient='left', unitNames=None, fontsize=11, weight='normal', color='k', font='Arial')

        else:
            self.picker = picker.Picker()
            self.picker.setData(2, self.scp)
            self.picker.setAction(self.handle_event)

        # self.figure.canvas.mpl_connect('button_press_event', self.picker.pickEvent)
        self.figure.canvas.mpl_connect('pick_event', self.picker.pickEvent)
        # self.figure.canvas.mpl_connect('motion_notify_event', self.picker.onMouseMotion)

        cp = self.cell.parts
        cellname = '/'.join(cp[-4:])
        self.figure.suptitle(cellname, fontsize=11)
        self.fig2 = None
        mpl.show()
        
        
    def show_traces(self, f, ax, pcolor='r', name=None):

        self.cell.glob('*')
        # imageplotted = False
        # imagetimes = []
        # imagename = []
        # maptimes = []
        # mapname = []
        supindex = self.AR.readDirIndex(currdir=self.cell)
    
        self.SI = ScannerInfo(self.AR)
        if self.invert:
            cmap = 'gist_gray_r'
        else:
            cmap = 'gist_gray'
            

        max_camera = None
        if self.averageScannerImages:
            max_camera = self.AR.getAverageScannerImages(dataname='Camera/frames.ma', mode='max', firstonly=False, limit=None)
            ax.imshow(max_camera, aspect='equal', cmap='Reds', alpha=0.7, vmin = 1000, vmax=self.vmax,
                extent=[np.min(self.SI.boxw[0]), np.max(self.SI.boxw[0]), np.min(self.SI.boxw[1]), np.max(self.SI.boxw[1])])
        # max_camera = scipy.ndimage.gaussian_filter(max_camera, sigma=256/(4.*10))
        if len(self.videos) > 0:
            self.process_videos()
            ax.imshow(self.merged_image, aspect='equal', cmap=cmap, alpha=0.75, vmin = 0, vmax=self.vmax,
            extent=[np.min(self.SI.boxw[0]), np.max(self.SI.boxw[0]), np.min(self.SI.boxw[1]), np.max(self.SI.boxw[1])])
            
        else:
            ax.imshow(self.image_data, aspect='equal', cmap=cmap, alpha=0.75, vmin = 0, vmax=self.vmax,
            extent=[np.min(self.SI.boxw[0]), np.max(self.SI.boxw[0]), np.min(self.SI.boxw[1]), np.max(self.SI.boxw[1])])

        self.scp = self.SI.scannerpositions
        scp = self.scp
        xmin = np.min(scp[:,0])
        xmax = np.max(scp[:,0])
        ymin = np.min(scp[:,1])
        ymax = np.max(scp[:,1])
        # print(xmin, ymin, xmax, ymax)
        ax.scatter(scp[:,0], scp[:,1], s=4, c='c', marker='o', alpha=0.3, picker=5)
        print('getdata: ', name, self.datasets)
        d = self.AR.getData()
        if name is not None:
            self.datasets[name] = self.AR.data_array

        if self.AR.mode in ['v', 'V', 'VC']:
            self.vscale = 1e5
            self.off = self.voff
        else:
            self.vscale = 1e-3
            self.off = self.ioff
        # print(len(self.AR.traces))
        tb = self.AR.time_base
        im0 = np.argmin(np.fabs(tb - self.twin[0]))
        im1 = np.argmin(np.fabs(tb - self.twin[1]))
        self.im0 = im0
        self.im1 = im1
        self.tb = tb[im0:im1]-tb[im0]
        # just plot as many as we have!
        dshape = self.datasets[name].shape
        # for p in range(dshape[0]): # scp.shape[0]):
        #     self._plot_one(ax, p, pcolor, name=name)
        xcal = self.xscale*3.5e-5*self.calbar[0]*1.25
        ycal = self.yscale*self.vscale*self.calbar[1]*0.5
        zero = 0

        # ax.plot(self.xscale*3.5e-5*np.array([0., 0., self.calbar[0]])+xmin - xcal,
        #        (self.yscale*self.vscale*(np.array([self.calbar[1],  0., 0. ])-zero))+self.off*self.vscale+ymin - ycal, 'k-', linewidth=1)
        # ax.text(xmin-xcal, ymin-ycal, f"{self.calbar[0]*1e3} ms\n{self.calbar[1]*1e12} pA",
        # verticalalignment='top', horizontalalignment='center', fontsize=8)

    def _plot_one(self, ax, p, pcolor, name=None, yscaleflag=True, tscale=True, offflag=True, ystep = 0., ythick=0.3):
        zero = 0.
        vdat = FILT.SignalFilter_LPFBessel(self.datasets[name][p, :], 2000.,
                        samplefreq=self.AR.sample_rate[0], NPole=8)
        if self.basezero:
            zero = np.mean(vdat[self.im0:self.im0+20])
        if offflag:
            xoff = self.scp[p,0]
            yoff = self.off*self.vscale+self.scp[p, 1]
        else:
            xoff = 0.
            yoff = 0.
        if yscaleflag:
            y_scale = self.yscale*self.vscale
        else:
            y_scale = 1e3
        yoff += ystep
        if tscale:
            ts = self.xscale*3.5e-5
        else:
            ts = 1e3
        ax.plot(ts*self.tb+xoff, (y_scale*(vdat[self.im0:self.im1]-zero))+yoff, color=pcolor, linewidth=ythick)

    def handle_event(self, index):
        # print('handle event index: ', index)
        # print(self.SI.scannerpositions[index,:])
        if index in self.indicesplotted:
            return
        if self.overlay:
            ystep = -self.nspots*10.
        self.palette = itertools.cycle(sns.color_palette("colorblind", 10))
        for i, name in enumerate(list(self.datasets.keys())):
            c = next(self.palette)
            self._plot_one(self.ax2, index, pcolor=c, name=name, yscaleflag=False, offflag=False,
                ystep=ystep, ythick=1.3-(i+1)*0.3, tscale=False)
            # if i == 1:
            #     self._plot_one(self.ax3, index, pcolor=c, name=name, yscaleflag=False, offflag=False,
            #         ystep=ystep, ythick=0.5, tscale=False)
        # self.xscale*3.5e-5*self.tb+xoff, (y_scale*(vdat[self.im0:self.im1]-zero))+yoff, color=pcolor, linewidth=0.3)
        trn = self.traces.index(index)+1
        self.ax.text(self.SI.scannerpositions[index][0], self.SI.scannerpositions[index][1],
                f"{trn:d}", fontsize=9, horizontalalignment='center')
        self.ax2.text(0., ystep,
                f"{trn:d}", fontsize=9, horizontalalignment='right')
        self.nspots += 1
        self.indicesplotted.append(index)
        mpl.draw()
        
if __name__ == '__main__':

    basepath = '/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107ai32Het/'
    MT = MapTraces()
    #    AR.setProtocol('/Users/pbmanis/Documents/data/MRK_Pyramidal/2018.01.26_000/slice_000/cell_000/CCIV_1nA_max_000/')
        # this won't work in the wild, need appropriate data for testing.
    # test on a big file    
    #cell = '/Users/pbmanis/Documents/data/mrk/2017.09.12_000/slice_000/cell_001'
    # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.06_000/slice_002/cell_000/LSPS_dendrite_VC_testmap_MAX_000_001')
    # image = 'image_002.tif'
    # # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.06_000/slice_002/cell_000/LSPS_dendrite_CC_testmap_strongest_000')
    # # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.05_000/slice_002/cell_000/LSPS_dendrite_VC_testmap_MAX_004')
    # # image = Path(cell.parent, 'image_001.tif')
    # cell = Path(basepath, '2017.08.22_000/slice_000/cell_001/Map_NewBlueLaser_VC_10Hz_000')  # pyr
    # image = '../image_008.tif'
    # MT.setPars({'invert': False, 'vmax': 30000, 'xscale': 1.5, 'yscale': 1.5, 'calbar': [0.5, 200.e-12]})  # calbar in ms, pA

    ctl = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_000')
    apv = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_002')
    image = '../image_002.tif'
    # tr = [111,148,28,36,96,117,140,61,85,165,65,173,89,119]
    tr = [111, 148, 140, 65, 89, 85, 36]
    tr = [111, 148, 36, 140, 85, 65, 89]  # reorder more or less top to bottom
    MT.setPars({'invert': True, 'vmax': 30000, 'xscale': 6, 'yscale': 1.5, 'ioff': -0.00,
            'calbar': [0.5, 200.e-12], 'twin': [0.25, 0.5]})  # calbar in ms, pA
    prots = {'ctl': ctl, 'apv': apv}
    # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_003')
    # image = '../image_002.tif'
    # MT.setPars({'invert': True, 'vmax': 30000, 'xscale': 6, 'yscale': 1.5, 'calbar': [0.5, 20e-3], 'twin': [0.25, 0.4], 'voff': -0.0})  # calbar in ms, pA

    # bushy
    # cell = Path(basepath, '2017.03.01_000/slice_000/cell_001/Map_NewBlueLaser_VC_single_test_001')  # pyr
    # image = '../image_001.tif'
    #
    # cell = Path(basepath, '2017.03.22_000/slice_000/cell_000/Map_NewBlueLaser_VC_pt2mW_000')  # pyr
    # image = '../image_000.tif'
    # # videos = [0, 1, 2, 3, 4]
    # cell = Path(basepath, '2017.03.24_000/slice_001/cell_001/Map_NewBlueLaser_VC_2mW_005')  # pyr
    # image = '../image_001.tif'
    # MT.setPars({'invert': False, 'vmax': 30000, 'xscale': 1.5, 'yscale': 0.05, 'calbar': [0.5, 5000.e-12]})  # calbar in ms, pA
    # MT.setPars({'invert': True, 'vmin': 1000, 'vmax': 18000, 'xscale': 6, 'yscale': 1.5, 'calbar': [0.5, 20.e-3], 'twin': [0.25, 0.5],
             # 'ioff': -0.0})  # calbar in ms, pA
    # cell1 = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_000')
    # cell2 = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_001')
    # cell3 = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_CC_testmap_MAX_002')
    # cell_vc = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_VC_testmap_MAX_000')
    # MT.setPars({'invert': True, 'vmin': 1000, 'vmax': 18000, 'xscale': 6, 'yscale': 1.5, 'calbar': [0.5, 200e-12], 'twin': [0.25, 0.5],
    #          })
    # image = '../image_002.tif'
    # prots = [cell1, cell2, cell3]
    # prots = [cell_vc]

    # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.05_000/slice_002/cell_000/LSPS_dendrite_VC_testmap_MAX_001')  # pyr
 #    image = '../image_001.tif'
 #    MT.setPars({'invert': True, 'vmax': 30000, 'xscale': 1.5, 'yscale': 1.5, 'calbar': [0.5, 200.e-12]})  # calbar in ms, pA
 #    prots = [cell]

    # cell1 = Path('/Users/pbmanis/Desktop/Data/CN Glu uncaging CBA/2017.12.01_000/slice_003/cell_000/Map_NewBlueLaser_VC_Single_008')
    # # cell2 = Path('/Users/pbmanis/Desktop/Data/CN Glu uncaging CBA/2017.12.01_000/slice_003/cell_000/Map_NewBlueLaser_VC_Single_009')
    # # cell2 = Path('/Users/pbmanis/Desktop/Data/CN Glu uncaging CBA/2017.12.01_000/slice_003/cell_000/Map_NewBlueLaser_VC_Single_010')
    # # cell3 = Path('/Users/pbmanis/Desktop/Data/CN Glu uncaging CBA/2017.12.01_000/slice_003/cell_000/Map_NewBlueLaser_VC_Single_011')
    # # cell4 = Path('/Users/pbmanis/Desktop/Data/CN Glu uncaging CBA/2017.12.01_000/slice_003/cell_000/Map_NewBlueLaser_VC_Single_012')
    # image = '../../image_001.tif'
    # MT.setPars({'invert': False, 'vmax': 30000, 'xscale': 6, 'yscale': 1.5, 'calbar': [0.5, 200.e-12], 'twin': [0.25, 0.4]})  # calbar in ms, pA
    # prots = [cell1] #, cell2, cell3, cell4]

    MT.plot_maps(prots, traces=tr)
    