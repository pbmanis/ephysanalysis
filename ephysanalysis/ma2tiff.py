"""
Translate metaarray file to tiff file
11 July 2018 PBM

Specify path in this file, and then pass the filename on the command line
The output file will be written to the same directory as the input file.

"""

import os
import sys
import fnmatch
import pyqtgraph.metaarray as MA
import tifffile as tf
import numpy as np
import glob

path = "/Users/pbmanis/Documents/data/MRK_Pyramidal/2017.10.04_000/slice_000"
pathname = sys.argv[1]

# python >3.5 has an automatic way to do this, but this works in Python 2.7
# the walk will inspect all subdirectories under the top directory in path
#
matches = []
for root, dirnames, filenames in os.walk(path):
    for filename in fnmatch.filter(filenames, 'video_*.ma'):
        matches.append(os.path.join(root, filename))

for m in matches:
    print m

    fullfile = m
    data = MA.MetaArray(file=fullfile)
    #print data
    ofile, ext = os.path.splitext(fullfile)
    ofile = ofile + '.tif'
    print('output file: ', ofile)

    fh = open(ofile, 'w')
    tf.imsave(fh, data.view(np.ndarray).astype('float32'), imagej=True, software='acq4')
    fh.close()
        



