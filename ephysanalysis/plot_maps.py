from __future__ import print_function
#!/usr/bin/python

"""
Simple adjunct routine to plot LSPS/CRACM maps with traces, over cell image if possible.
Reuqires acq4read. 

Takes advantage of acq4read code to access all data and attributes.


"""
import os
import re
from pathlib import Path
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as mpl
import ephysanalysis.acq4read as ARC
import ephysanalysis.metaarray as EM
from pyqtgraph import configfile
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
        fp = np.array([scannerbox[0][0], scannerbox[1][1]]).reshape(2,1)
        scannerbox = np.append(scannerbox, fp, axis=1)
        self.scboxw = np.array(scannerbox)        
        self.boxw = np.swapaxes(np.array(self.camerabox), 0, 1)
        

class MapTraces(object):
    def __init__(self):
        self.cell = None
        self.image = None
        self.AR = ARC.Acq4Read()
        self.invert = True
        self.vmax = 20000.
        self.voff = 0.
        self.ioff = 0.050
        self.basezero = True
        self.xscale = 1.0
        self.yscale = 1.0
        self.twin = [0, 0.6]
        self.averageScannerImages = False # normally, would not do
        self.calbar = [0.1, 5e-7]  # 100 ms, 500 pA

    def setProtocol(self, cell, image, videos=None):
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

    def show_traces(self):

        datasets = self.cell.glob('*')
        # imageplotted = False
        # imagetimes = []
        # imagename = []
        # maptimes = []
        # mapname = []
        supindex = self.AR.readDirIndex(currdir=self.cell)
        f, ax = mpl.subplots(1, 1)
    
        SI = ScannerInfo(self.AR)
        if self.invert:
            cmap = 'gist_gray_r'
        else:
            cmap = 'gist_gray'
            

        max_camera = None
        if self.averageScannerImages:
            max_camera = self.AR.getAverageScannerImages(dataname='Camera/frames.ma', mode='max', firstonly=False, limit=None)
            ax.imshow(max_camera, aspect='equal', cmap='Reds', alpha=0.7, vmin = 1000, vmax=self.vmax,
                extent=[np.min(SI.boxw[0]), np.max(SI.boxw[0]), np.min(SI.boxw[1]), np.max(SI.boxw[1])])
        # max_camera = scipy.ndimage.gaussian_filter(max_camera, sigma=256/(4.*10))
        if len(self.videos) > 0:
            self.process_videos()
            ax.imshow(self.merged_image, aspect='equal', cmap=cmap, alpha=0.75, vmin = 0, vmax=self.vmax,
            extent=[np.min(SI.boxw[0]), np.max(SI.boxw[0]), np.min(SI.boxw[1]), np.max(SI.boxw[1])])
            
        else:
            ax.imshow(self.image_data, aspect='equal', cmap=cmap, alpha=0.75, vmin = 0, vmax=self.vmax,
            extent=[np.min(SI.boxw[0]), np.max(SI.boxw[0]), np.min(SI.boxw[1]), np.max(SI.boxw[1])])

        scp = SI.scannerpositions
        xmin = np.min(scp[:,0])
        xmax = np.max(scp[:,0])
        ymin = np.min(scp[:,1])
        ymax = np.max(scp[:,1])
        print(xmin, ymin, xmax, ymax)
        ax.plot(scp[:,0], scp[:,1], 'co', alpha=0.3, markersize=4, markeredgecolor='c')
        d = self.AR.getData()

        if self.AR.mode in ['v', 'V', 'VC']:
            vscale = 1e5
            voff = self.voff
        else:
            vscale = 1e-3
            voff = self.ioff
        print(len(self.AR.traces))
        tb = self.AR.time_base
        im0 = np.argmin(np.fabs(tb - self.twin[0]))
        im1 = np.argmin(np.fabs(tb - self.twin[1]))
        tb = tb[im0:im1]-tb[im0]
        for p in range(scp.shape[0]):
            vdat = FILT.SignalFilter_LPFBessel(self.AR.data_array[p, :], 2000.,
                            samplefreq=self.AR.sample_rate[0], NPole=8)
            zero = 0.
            if self.basezero:
                zero = np.mean(vdat[im0:im0+20])
            ax.plot(self.xscale*3.5e-5*tb+scp[p,0], (self.yscale*vscale*(vdat[im0:im1]-zero))+voff*vscale+scp[p, 1], 'r-', linewidth=0.3)
        xcal = self.xscale*3.5e-5*self.calbar[0]*1.25
        ycal = self.yscale*vscale*self.calbar[1]*0.5
        ax.plot(self.xscale*3.5e-5*np.array([0., 0., self.calbar[0]])+xmin - xcal, 
               (self.yscale*vscale*(np.array([self.calbar[1],  0., 0. ])-zero))+voff*vscale+ymin - ycal, 'k-', linewidth=1)
        ax.text(xmin-xcal, ymin-ycal, f"{self.calbar[0]*1e3} ms\n{self.calbar[1]*1e12} pA", 
        verticalalignment='top', horizontalalignment='center', fontsize=8)
        cp = self.cell.parts
        cellname = '/'.join(cp[-4:])
        f.suptitle(cellname, fontsize=11)
        mpl.show()


if __name__ == '__main__':

    basepath = '/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107ai32Het/'
    MT = MapTraces()
    #    AR.setProtocol('/Users/pbmanis/Documents/data/MRK_Pyramidal/2018.01.26_000/slice_000/cell_000/CCIV_1nA_max_000/')
        # this won't work in the wild, need appropriate data for testing.
    # test on a big file    
    #cell = '/Users/pbmanis/Documents/data/mrk/2017.09.12_000/slice_000/cell_001'
    cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.06_000/slice_002/cell_000/LSPS_dendrite_VC_testmap_MAX_000_001')
    image = 'image_002.tif'
    # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.06_000/slice_002/cell_000/LSPS_dendrite_CC_testmap_strongest_000')
    # cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.05_000/slice_002/cell_000/LSPS_dendrite_VC_testmap_MAX_004')
    # image = Path(cell.parent, 'image_001.tif')
    cell = Path(basepath, '2017.08.22_000/slice_000/cell_001/Map_NewBlueLaser_VC_10Hz_000')  # pyr
    image = '../image_008.tif'
    MT.setPars({'invert': False, 'vmax': 30000, 'xscale': 1.5, 'yscale': 1.5, 'calbar': [0.5, 200.e-12]})  # calbar in ms, pA

    cell = Path('/Users/pbmanis/Desktop/Data/Glutamate_LSPS_DCN/2019.08.08_000/slice_001/cell_000/LSPS_dendrite_VC_testmap_MAX_000')
    image = '../image_002.tif'
    MT.setPars({'invert': True, 'vmax': 30000, 'xscale': 6, 'yscale': 1.5, 'calbar': [0.5, 200.e-12], 'twin': [0.25, 0.4]})  # calbar in ms, pA

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
    
    MT.setProtocol(cell, image)
    MT.show_traces()
    