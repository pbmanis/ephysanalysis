#!/usr/bin/env python3
from __future__ import print_function
"""
Translate metaarray file to tiff file
11 July 2018 PBM

Specify path in this file, and then pass the filename on the command line
The output file will be written to the same directory as the input file.

Note we set bigtiff true because some of the tiff files will be > 4 GB

"""

import os
import sys
python_version = sys.version_info[0]
import fnmatch
import pylibrary.tifffile as tf
import numpy as np


def convertfiles():
# define path to data
    path = "/Users/pbmanis/Desktop/Python/mrk_nf107" # Documents/data/MRK_Pyramidal/2017.10.04_000/slice_000"
    if len(sys.argv) >= 2:
        pathname = sys.argv[1]  # really does nothing right now.

    if python_version == 2:
        # python >3.5 has an automatic way to do this, but this works in Python 2.7
        # the walk will inspect all subdirectories under the top directory in path
        #
        import glob
        import pyqtgraph.metaarray as MA
        matches = []
        for root, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, 'video_*.ma'):
                matches.append(os.path.join(root, filename))

    else:  # python 3
        from ephysanalysis import metaarray as MA
        from pathlib import Path
        p = Path('.')
        matches = list(p.glob('video_*.ma'))

    print('files: ', matches)

    for m in matches:

        fullfile = str(m)
        data = MA.MetaArray(file=fullfile)
        ofile, ext = os.path.splitext(fullfile) # cleaner with pathlib... 
        ofile = ofile + '.tif'
        print('output file: ', ofile)

        fh = open(ofile, 'wb')
        tf.imsave(fh, data.view(np.ndarray).astype('float32'), imagej=True, bigtiff=True, software='acq4')
        fh.close()
        
if __name__ == '__main__':
    convertfiles()



