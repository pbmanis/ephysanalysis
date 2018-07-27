from __future__ import print_function

__author__ = 'pbmanis'
"""
dir_check.py
Check an acq4 directory struture to be sure that the hierarchy is correct.
day
    slice
        cell
            protocols
        pair
            cells
                protocols

The clice may have images and videos, but no protocols
The cells and paris may have both
Pairs should have cells and protocols under cells

Any other organization is flagged.


"""
import sys
import os
import re
import math  # use to check nan value...
import argparse
import os.path
import stat
import gc
import argparse
import datetime
import dateutil.parser as DUP
from termcolor import colored
import numpy as np
import textwrap
from collections import OrderedDict
import pandas as pd
import pandas.compat # for StringIO

from ephysanalysis import acq4read
from pyqtgraph.metaarray import MetaArray

class Printer():
    """Print things to stdout on one line dynamically"""
    def __init__(self, data):
        sys.stdout.write("\033[1;36m\r\x1b[K\033[0;0m"+data.__str__())
        sys.stdout.flush()


class DirCheck():
    def __init__(self, topdir):
        """
        Check directory structure
        """
        self.topdir = topdir
        
        self.img_re = re.compile('^[Ii]mage_(\d{3,3}).tif')  # make case insensitive - for some reason in Xuying's data
        self.s2p_re = re.compile('^2pStack_(\d{3,3}).ma')
        self.i2p_re = re.compile('^2pImage_(\d{3,3}).ma')
        self.video_re = re.compile('^[Vv]ideo_(\d{3,3}).ma')
        self.daytype = re.compile("(\d{4,4}).(\d{2,2}).(\d{2,2})_(\d{3,3})")
        self.tstamp = re.compile('(    __timestamp__: )([\d.\d]*)')
        
        if not os.path.isdir(self.topdir):
            raise ValueError('top dir does not exist: %s' % self.topdir)

        metaarraytypes = []
        protocols = []
        level = 0
        parent = ''
        lastcb = ''
        # walk the directory from top down - but in a highly controlled manner:
        
        topdirs = os.listdir(topdir)
        for d in sorted(topdirs):
            if d in ['.index', '.DS_Store', 'log.txt']:
                continue
            if any([d.endswith(e) for e in ['.tif', '.ma']]):
                continue
            m = self.daytype.match(d)
            tstamp = self.gettimestamp(os.path.join(self.topdir, d))
            print('-'*80)
            if m is not None:
                print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}'.format(d, '', '', tstamp), 'white'))
            else:
                print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}   is not a DAY directory'.format(d, '', '', tstamp), 'red'))
            for s in sorted(os.listdir(os.path.join(topdir, d))):
                if s in ['.index', '.DS_Store', 'log.txt']:
                    continue
                tstamp = self.gettimestamp(os.path.join(self.topdir, d, s))
                if any([s.endswith(e) for e in ['.tif', '.ma']]):
                    st = os.stat(os.path.join(self.topdir, d, s))  # unmanaged (though may be in top index file)
                    tstamp = datetime.datetime.fromtimestamp(st[stat.ST_MTIME]).strftime('%Y-%m-%d  %H:%M:%S %z')
                    print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}  data file not associated with slice or cell'.format('', s, '', tstamp), 'cyan'))
                    continue
                if s.startswith('slice_'):
                    print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}'.format('', s, '', tstamp), 'white'))
                else:
                    print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}   is not a SLICE directory'.format('', s, '', tstamp), 'red'))
                
                for c in sorted(os.listdir(os.path.join(topdir, d, s))):
                    if c in ['.index', '.DS_Store', 'log.txt']:
                        continue
                    if any([c.endswith(e) for e in ['.tif', '.ma']]):
                        continue
                    tstamp = self.gettimestamp(os.path.join(self.topdir, d, s, c))
                    if c.startswith('cell_'):
                        print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}'.format('', '', c, tstamp), 'white'))
                    else:
                        print(colored('{0:>15s} {1:>10s} {2:>10s} {3:>20}  is not a CELL directory'.format('', '', c, tstamp), 'red'))
        print(colored('------', 'white'))

    def gettimestamp(self, path):
        """
        Get the timestamp of an .index file
        """
        tstamp = 'None'
        indexfile = os.path.join(path, '.index')
        if not os.path.isfile(indexfile):
            return tstamp
        with open(indexfile, 'r') as f:
            for l in f:
                ts = self.tstamp.match(l)
                if ts is not None:
                    fts = float(ts.group(2))
                    tstamp = datetime.datetime.fromtimestamp(fts).strftime('%Y-%m-%d  %H:%M:%S %z')
#                            print('time: ', tstamp)
        return tstamp



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Data Summaries from acq4 datasets')
    parser.add_argument('basedir', type=str,
                        help='Base Directory')
    parser.add_argument('-r', '--read', action='store_true', dest='read',
                        help='just read the protocol')
    args = parser.parse_args()
    DirCheck(args.basedir)
    