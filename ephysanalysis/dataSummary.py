from __future__ import print_function

__author__ = 'pbmanis'
"""
dataSummary: This script reads all of the data files in a given directory, and prints out top level information
including notes, protocols run (and whether or not they are complete), and image files associated with a cell.
Currently, this routine makes assumptions about the layout as a hierarchical structure [days, slices, cells, protocols]
and does not print out information if there are no successful protocols run.
June, 2014, Paul B. Manis.

Mar 2015:
added argparse to expand command line options rather than editing file. 
The following options are recognized:
begin (b) (define start date; end is set to current date) default: 1/1/1970
end (e)(define end date: start date set to 1/1/1970; end is set to the end date) default: "today"
mode =  full (f) : do a full investigation of the data files. Makes processing very slow. (reports incomplete protocols)
        partial (p) : do a partial investiagion of protocols: is there anything in every protocol directory? (reports incomplete protocols) - slow
        quick (q) : do a quick scan : does not run through protocols to find incomplete protocols. Default (over full and partial)
debug (d) : debug monitoring of progress
output (o) : define output file (tab delimited file for import to other programs)


Mar 2018: Version 2
Uses acq4read and is independent of acq4 itself.

July 2018: 
Major surgery - to output Pandas (pickled) files as well. UGH.

"""
import sys
import os
import re
import math  # use to check nan value...
import argparse
import os.path
import gc
import argparse
import datetime
import dateutil.parser as DUP
import numpy as np
import textwrap
from collections import OrderedDict
import pandas as pd
import pandas.compat # for StringIO

from ephysanalysis import acq4read
from pyqtgraph.metaarray import MetaArray

def ansi_colors(color):
    colors = {
    'black': '\u001b[30m',
    'red': '\u001b[31m',
    'green': '\u001b[32m',
    'yellow': '\u001b[33m',
    'blue': '\u001b[34m',
    'magenta': '\u001b[35m',
    'cyan': '\u001b[36m',
    'white': '\u001b[37m',
    'reset': '\u001b[0m',
    }
    return colors[color]
    
class Printer():
    """Print things to stdout on one line dynamically"""
    def __init__(self, data, color='white'):
        
        sys.stdout.write(u"\r\u001b[2K%s"%ansi_colors(color)+data.__str__())
        sys.stdout.flush()

class DataSummary():

    def __init__(self, basedir, outputMode='terminal', outputFile=None, daylistfile=None,
                 after=None, before=None, dryrun=False, depth='all', inspect=True):
        """
        Note that the init is just setup - you have to call getDay the object to do anything
    
        Parameters
        ----------
        basedir : str (required)
            base directory to be summarized
    
        outputMode : str (default: 'terminal')
            How to write the output. Current options are 'terminal', which writes to the terminal, and
            'pandas', which will create a pandas dataframe and pickle it to a file
    
        daylistfile : str (default: None)
            A filename that stores a list of days that sould be processed
    
        after : str (default: Jan 1, 1970)
            A starting date for the summary - files before this date will not be processed
            A string in a format that can be parsed by dateutil.parser.parse
            The following will work:
            without string quotes:
            2018.01.01
            2018.1.1
            2018/7/6
            2018-7-6
            With quotes:
            'Jan 1 2017'
            'Jan 1, 2017'
            etc..
    
        before : str (default 2266)
            An ending date for the summary. Files after this date will not be processed. 
            The format is the same as for the starting date
    
        dryrun : bool (default False)
            Causes the output ti be limited and protocols are not fully investigated
    
        depth : str (default: 'all')
            Causes output with dry-run to go to the depth specified
            options are day, slice, cell, prototol
    
        Note that if neither before or after are specified, 
        """

        self.setups()  # build some regex and wrappers
        # gather input parameters
        self.basedir = basedir
        self.outputMode = outputMode  # terminal, tabfile, pandas
        self.outFilename = outputFile
        self.daylistfile = daylistfile
        self.dryrun = dryrun
        self.after = after
        self.before = before
        self.depth = depth
 
        self.daylist = None
        self.index = 0
        # flags - for debugging and verbosity
        self.monitor = False
        self.reportIncompleteProtocols = True  # do include incomplete protocol runs in print
        self.InvestigateProtocols = inspect  # set True to check out the protocols in detail
        if self.dryrun:
            self.monitor = True
            self.reportIncompleteProtocols = False # do include incomplete protocol runs in print
            self.InvestigateProtocols = False  # set True to check out the protocols in detail
        self.panda_string = ''
        
        # initialized dictionary that holds all the stuff
        #self.analysis_summary = {}

        # column definitions - may need to adjust if change data that is pasted into the output
        self.day_defs = ['date', 'description', 'notes', 'species', 'strain', 'genotype', 'age', 
                         'sex', 'weight', 'solution', 'internal', 'temperature', 'important', 'expUnit']
        self.slice_defs = ['slice_slice', 'slice_notes', 'slice_location', 'slice_orientation', 'important']
        self.cell_defs =  ['cell_cell',   'cell_notes',  'cell_type', 'cell_location', 'cell_important']
        self.data_defs = ['data_incomplete', 'data_complete', 'data_images']
        
        # expected keys in various structural levels: days, slices, cells
        self.day_keys = ['description', 'notes', 'species', 'strain', 'genotype', 'age', 'sex', 'weight', 'solution', 'animal identification', '__timestamp__', 
                        'internal', 'temperature', 'expUnit', 'dirType', 'important', 'time']
        self.slice_keys = ['notes', 'location', 'orientation', 'important', '__timestamp__']
        self.cell_keys = [ 'notes', 'type', 'location', 'important', '__timestamp__']
        self.data_dkeys = ['incomplete', 'complete', 'data_images']
        
        self.day_template = (
            OrderedDict([('species', '{:>s}'), ('strain', '{:>s}'),('genotype', '{:>12s}'), ('age', '{:>5s}'), ('sex', '{:>1s}'), ('weight', '{:>5s}'),
                         ('solution', '{:>s}'), ('internal', '{:>s}'), ('temperature', '{:>5s}'), ('important', '{:>s}'), ('elapsedtime', '{:>8.2f}'), 
                         ('expUnit', '{:>s}')])
                        )        
        self.slice_template = (
            OrderedDict([('type', '{:>s}'), ('location', '{:>12s}'), ('orientation', '{:>5s}')])
                        )        
        self.cell_template = (
            OrderedDict([('type', '{:>s}'), ('location', '{:>12s}'), ('important', '{:>s}')])
                        )        
        self.data_template = (
            OrderedDict([('incomplete', '{0:s}'), ('complete', '{1:s}'), ('images', '{2:s}')])
                        )        
        self.AR = acq4read.Acq4Read()  # instance of the reader

#        outputDir = os.path.join(os.path.expanduser("~"), 'Desktop/acq4_scripts')
        if self.outputMode == 'tabfile':
            # self.outFilename = self.basedir.replace('/', '_') + '.tab'
            # self.outFilename = self.outFilename.replace('\\', '_')
            # if self.outFilename[0] == '_':
            #     self.outFilename = self.outFilename[1:]
            # self.outFilename = os.path.join(outputDir, self.outFilename)
            print('Tabfile output: Writing to {:<s}'.format(self.outFilename))
            h = open(self.outFilename, 'w')  # write new file
            h.write(self.basedir+'\n')
            h.write(self.coldefs + '\n')
            h.close()
        elif self.outputMode == 'pandas':  # save output as a pandas data structure, pickled
            # self.outFilename = self.basedir.replace('/', '_') + '.pkl'
            # self.outFilename = self.outFilename.replace('\\', '_')
            # if self.outFilename[0] == '_':
            #     self.outFilename = self.outFilename[1:]
            # self.outFilename = os.path.join(outputDir, self.outFilename)
            print('Pandas output: will write to {:<s}'.format(self.outFilename))
        else:
            pass  # just print if terminal
            
        # figure out the before/after limits
        # read the before and after day strings, parse them and set the min and maxdays
        if self.after is None:
            mindayx = (1970, 1, 1)  # all dates after this - unix start date convention
        else:
            try:
                dt = DUP.parse(self.after)
                mindayx = (dt.year, dt.month, dt.day)
            except:
                raise ValueError('Date for AFTER cannot be parsed : {0:s}'.format(self.after))
        if self.before is None:
            maxdayx = (2266, 1, 1)  # far enough into the future for you? Maybe when the Enterprise started it's journey?
        else:
            try:
                dt = DUP.parse(self.before)
                maxdayx = (dt.year, dt.month, dt.day)
            except:
                raise ValueError('Date for BEFORE cannot be parsed : {0:s}'.format(self.before))
        
        if self.daylistfile is None:  # get from command line
            self.minday = mindayx[0]*1e4+mindayx[1]*1e2+mindayx[2]
            self.maxday = maxdayx[0]*1e4+maxdayx[1]*1e2+maxdayx[2]
        else:
            self.daylist = []
            with open(self.daylistfile, 'r') as f:
                for line in f:
                    if line[0] != '#':
                        self.daylist.append(line[0:10])
            f.close()  # silly - it dhould be closed.
        
    def setups(self):
        self.tw = {}  # for notes
        self.tw['day'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)  # used to say "initial_indent ="Description: ""
        self.tw['slice'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)
        self.tw['cell'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)

        self.twd = {}  # for description
        self.twd['day'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)  # used to ays initial_indent ="Notes: ""
        self.twd['slice'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)
        self.twd['cell'] = textwrap.TextWrapper(initial_indent="", subsequent_indent=" "*2)
        
        self.img_re = re.compile('^[Ii]mage_(\d{3,3}).tif')  # make case insensitive - for some reason in Xuying's data
        self.s2p_re = re.compile('^2pStack_(\d{3,3}).ma')
        self.i2p_re = re.compile('^2pImage_(\d{3,3}).ma')
        self.video_re = re.compile('^[Vv]ideo_(\d{3,3}).ma')

        self.daytype = re.compile("(\d{4,4}).(\d{2,2}).(\d{2,2})_(\d{3,3})")
#        daytype = re.compile("(2011).(06).(08)_(\d{3,3})")  # specify a day
        
    def wrapstring(self, datastring, data, argument, wrapper, default='?'):
        """
        Wrap up text in a generic way, returning the modified string
        """
        if data is not None and argument in data.keys() and len(data[argument]) > 0:
            lstr = wrapper.wrap(data[argument])
            for i in lstr:
                nstr = i.replace('\t', '    ')  # clean out tabs so printed formatting is not confused
                datastring += nstr
            datastring = datastring.lstrip()
            return datastring + '\t'
        else:
            return datastring + default + '\t'

    def getDay(self):
        """
        getDay is the entry point for scanning through all the data files in a given directory,
        returning information about those within the date range, with details as specified by the options
        """
        print('getDay')
        allfiles = os.listdir(self.basedir)
        self.pstring = ''
        days = []
        for thisfile in allfiles:
            m = self.daytype.match(thisfile)
            if m in ['.DS_Store']:
                continue
            if m is None:
                continue  # no match
            if len(m.groups()) >= 3:  # perfect match
                idl = [int(d) for d in m.groups()]
                id = idl[0]*1e4+idl[1]*1e2+idl[2]

                if self.daylist is None:
                    if id >= self.minday and id <= self.maxday:
                        days.append(thisfile)  # was [0:10]
                else:
                    if thisfile[0:10] in self.daylist:
                        days.append(thisfile)
        if self.monitor:
            print ('Days reported: ', days)
        for nd, day in enumerate(days):
            if self.monitor:
                self.pstring = 'Processing day[%3d/%3d]: %s ' % (nd, len(days), day)
            self.AR.setProtocol(os.path.join(self.basedir, day))
            self.day_index = self.AR.getIndex()
            self.day_index['date'] = day.strip()
            #print('\nday index: ', self.day_index.keys())
            # now add the rest of the index information to the daystring
            for k in self.day_defs:
                if k not in self.day_index.keys():
                    #print('\nadded: ', k)
                    self.day_index[k] = ' '
                if isinstance(self.day_index[k], bool):
                    self.day_index[k] = str(self.day_index[k])
                self.day_index[k].replace('\n', ' ')
                if len(self.day_index[k]) == 0:
                    self.day_index[k] = ' '
            # for k in self.day_defs:
            #     print('{:>32s} : {:<40s}'.format(k, self.day_index[k]))

            self.doSlices(os.path.join(self.basedir, day))  # next level
            os.closerange(8, 65535)  # close all files in each iteration
            gc.collect()

    def doSlices(self, day):
        """
        process all of the slices for a given day
        :param day:
        :return nothing:
        """
#        print('\ndo slices')
        allfiles = os.listdir(day)
        slicetype = re.compile("(slice\_)(\d{3,3})")
        slices = []
        for thisfile in allfiles:
            m = slicetype.match(thisfile)
            if m is None:
                continue
            if len(m.groups()) == 2:
                slices.append(thisfile)
        for slicen in slices:
            self.sstring = self.pstring + ' %s' % slicen
            Printer(self.sstring)
            self.slicestring = '%s\t' % (slicen)
#            print('slicestring: ', self.slicestring)
            self.slice_index = self.AR.getIndex(os.path.join(day, slicen))
            if self.slice_index is None:  # directory is not managed and probably empty
                self.slice_index = {}
                continue
            self.slice_index['slice'] = slicen
            for k in self.slice_defs:
                ks = k.replace('slice_', '')
                if ks not in self.slice_index.keys():
                    self.slice_index[ks] = ' '
                if isinstance(self.slice_index[ks], bool):
                    self.slice_index[ks] = str(self.slice_index[ks])
                if len(self.slice_index[ks]) == 0 :
                    self.slice_index[ks] = ' '
                self.slice_index[ks].replace('\n', ' ')
            self.doCells(os.path.join(day, slicen))
            gc.collect()

    def doCells(self, thisslice):
        """
        process all of the cells from a slice
        :param slice:
        :return nothing:
        """
#        print('\ndo cells')
        allfiles = os.listdir(thisslice)
        cell_re = re.compile("(cell_)(\d{3,3})")
        cells = []
        for thisfile in allfiles:
            m = cell_re.match(thisfile)
            if m is None:
                continue
            if len(m.groups()) == 2:
                cells.append(thisfile)
        for cell in cells:
            self.cstring = self.sstring + " %s" % cell
            Printer(self.cstring)
            self.cell_index = self.AR.getIndex(os.path.join(thisslice, cell))  # possible that .index file is missing, so we cannot read
            if self.cell_index is None:
                self.cell_index = {}  # directory is not managed, so skip
                continue
#            print('\nCell Index: ', self.cell_index)
            self.cell_index['cell'] = cell
            for k in self.cell_defs:
                ks = k.replace('cell_', '')
                if ks not in self.cell_index.keys():
                    self.cell_index[ks] = ' '
                if isinstance(self.cell_index[ks], bool):
                    self.cell_index[ks] = str(self.cell_index[ks])
                self.cell_index[ks].replace('\n', ' ')
                if len(self.cell_index[ks]) == 0:
                    self.cell_index[ks] = ' '
#            print('\n cell index: ', self.cell_index)
            self.doProtocols(os.path.join(thisslice, cell))
            gc.collect()

    def doProtocols(self, thiscell):
        """
        process all of the protocols for a given cell
        :param cell:
        :return nothing:
        """
#        print( '\n\nSearching protocols')
        allfiles = os.listdir(thiscell)
        protocols = []
        nonprotocols = []
        anyprotocols = False
        images = []  # tiff
        stacks2p = []
        images2p = []
        videos = []

        endmatch = re.compile("[\_(\d{3,3})]$")  # look for _lmn at end of directory name
        for thisfile in allfiles:
            if os.path.isdir(os.path.join(thiscell, thisfile)):
                protocols.append(thisfile)
            else:
                nonprotocols.append(thisfile)

        self.protocolstring = ''
        self.allprotocols = []
        self.incompleteprotocols = []
        self.completeprotocols = []
        self.compprotstring = ''
        if self.InvestigateProtocols is True:
           # self.thiscell_summarystring = 'NaN\t'*6

            for np, protocol in enumerate(protocols):
                Printer(self.cstring + ' Prot[%2d/%2d]: %s' % (np, len(protocols), protocol))
                self.allprotocols += protocol+', '
                protocolpath = os.path.join(thiscell, protocol)
                # if np == 0:
               #      self.getDay(name=protocolpath)
                dirs = self.AR.subDirs(protocolpath)
                protocolok = True  # assume that protocol is ok
                modes = []
                nexpected = len(dirs)  # acq4 writes dirs before, so this is the expected fill
                ncomplete = 0  # count number actually done
                info = self.AR.getIndex(protocolpath)
                clampDevices = self.AR.getClampDevices(protocolpath)

                # must be able to handle multiple data formats, even in one experiment...
                data_mode = 'Missing mode information'
                if len(clampDevices) > 0 and 'devices' in info.keys():
                    data_mode = info['devices'][clampDevices[0]]['mode']  # get mode from top of protocol information
                else:  # try to set a data mode indirectly
                    protocolok = False  # can't parse protocol device...
                if data_mode not in modes:
                    modes.append(data_mode)
                
                protocol_truncated = False
                for i, directory_name in enumerate(dirs):  # dirs has the names of the runs within the protocol
                    if protocol_truncated or not protocolok:
                        break  # terminate, don't both continuing to search.
                    try:
                        info = self.AR.getIndex(directory_name)
                    except:
                       # print('dir %s seems not managed', directory_name)
                        protocol_truncated = True  # condition: missint .index file in protocol sequence dir.
                        continue
                    if info is not None:  # .index file is found, so proceed
                        datafile = os.path.join(directory_name, clampDevices[0]+'.ma')  # clamp device file name
                        clampInfo = self.AR.getDataInfo(datafile)
                        if clampInfo is None:
                            protocol_truncated = True  # condition: no valid clamp file in sequence dir
                            continue
                        if clampInfo is None:
                            protocol_truncated = True
                            continue
                        self.holding = self.AR.parseClampHoldingLevel(clampInfo)
                        self.amp_settings = self.AR.parseClampWCCompSettings(clampInfo)
                        if protocol not in self.completeprotocols:
                            ncomplete = ncomplete + 1
                            self.completeprotocols.append(protocol)
                if protocol_truncated:
                    self.incompleteprotocols.append(protocol)
                gc.collect()  # and force garbage collection of freed objects inside the loop
                if modes == []:
                    modes = ['Unknown mode']
                if protocolok and ncomplete == nexpected:  # accumulate protocols
                    self.protocolstring += '[{:<s}: {:s} {:d}], '.format(protocol, modes[0][0], ncomplete)
                    anyprotocols = True  # indicate that ANY protocol ran to completion
                else:
                    if self.reportIncompleteProtocols:
                        self.protocolstring += '{0:<s}.{1:s}.{2:d}/{3:d}, '.format(protocol, modes[0][0], ncomplete, nexpected)
                gc.collect()
        else:
            self.protocolstring += ''

            anyprotocols = True
            prots = {}
            for np, protocol in enumerate(protocols):
                Printer(self.cstring + ' Prot[%2d/%2d]: %s' % (np,len(protocols), protocol))
                if protocol not in self.allprotocols:
                    self.allprotocols.append(protocol)
                m = endmatch.search(protocol)
                if m is not None:
                    p = protocol[:-4]
                else:
                    p = protocol
                if p not in prots.keys():
                    prots[p] = 1
                else:
                    prots[p] += 1
            if len(prots.keys()) > 0:
                for p in prots.keys():
                    self.protocolstring += '{:<s}({:<d}), '.format(p, prots[p])
            else:
                self.protocolstring = '<No protocols found>'

        if len(self.completeprotocols) == 0:
            self.completeprotocols = ' '
        else:
            self.compprotstring = ', '.join([str(cp) for cp in self.completeprotocols])
        # print('\n: ', self.completeprotocols)
        # exit(1)
        self.allprotocols = ', '.join(self.allprotocols)
        for thisfile in nonprotocols:
            x = self.img_re.match(thisfile)  # look for image files
            if x is not None:
                images.append(thisfile)
            x = self.s2p_re.match(thisfile)  # two photon stacks
            if x is not None:
                stacks2p.append(thisfile)
            x = self.i2p_re.match(thisfile)  # simple two photon images
            if x is not None:
                images2p.append(thisfile)
            x = self.video_re.match(thisfile)  # video images
            if x is not None:
                videos.append(thisfile)
        self.imagestring = ''
        if len(images) > 0:
            self.imagestring += 'Images: %3d' % len(images)
        if len(stacks2p) > 0:
            self.imagestring += ', 2pStacks: %3d' % len(stacks2p)
        if len(images2p) > 0:
            self.imagestring += ', 2pImages: %3d' % len(images2p)
        if len(videos) > 0:
            self.imagestring += ', Videos: %3d' % len(videos)
        if len(images) + len(stacks2p) + len(images2p) + len(videos) == 0:
            self.imagestring = 'No Images or Videos'
        
        # if anyprotocols:
        #     ostring = self.daystring + self.thiscell_summarystring + self.slicestring + self.cellstring + self.protocolstring + self.completeprotocols + self.imagestring
        # else:
        ostring = OrderedDict([('incomplete', self.protocolstring.rstrip(', ')),  ('complete', self.compprotstring.rstrip(', ')), ('images', self.imagestring)])
        self.outputString(ostring)

    def colprint(self, phdr, ostring):
        ps = phdr.split('\t')
        os = ostring.split('\t')
        for i in range(len(ps)):
            if i > len(os):
                break
            if os[i] == ' ':
                os[i] = '--'
            print('{0:3d}: {1:>20s} : {2:<s}'.format(i, ps[i], os[i]))

    def outputString(self, ostring):
        day_string = ''
        phdr = ''
        for k in self.day_defs:
            day_string += str(self.day_index[k])+'\t'
            phdr += k +'\t'
        
        slice_string = ''
        for k in self.slice_defs:
            ks = k.replace('slice_', '')
            slice_string += str(self.slice_index[ks])+'\t'
            phdr += k +'\t'
        
        cell_string = ''
        for k in self.cell_defs:
            kc = k.replace('cell_', '')
            cell_string += str(self.cell_index[kc])+'\t'
            phdr += k +'\t'
        
        prot_string = ''
        ohdr = ''
        for k in self.data_defs:
            kc = k.replace('data_', '')
            pstx = str(ostring[kc])
            if len(pstx) == 0:
                pstx = ' '
            prot_string += pstx + '\t'
            phdr += k + '\t'
            ohdr += k + '\t'
        # print('\n\nestring: ', ostring)
        # print('\nprot_string: ', prot_string)
        # print('\nohdr: ', ohdr)
            
        ostring = day_string + slice_string + cell_string + prot_string
        ostring = ostring.replace('\n', ' ')
        ostring = ostring.rstrip('\t ')
        ostring += '\n'
        # print('\nOSTRING: \n', ostring.replace('\t', '##').replace('\n', '&&'))
        phdr = phdr.rstrip('\t\n')
        if len(self.panda_string) == 0:  # catch the header
            self.panda_string = phdr.replace('\n', '') + '\n'  # clip the last \t and add a 
        
        if self.outputMode in [None, 'terminal']:
            print('{0:s}'.format(ostring))
        
        elif self.outputMode == 'text':
            h = open(self.outFilename, 'a')  # append mode
            h.write(ostring )
            h.close()
        
        elif self.outputMode == 'pandas':
            # self.colprint(phdr, ostring)            
            # print('\n******* building Pandas string', 'ostring: \n', ostring)
            self.panda_string += ('{0:d}\t{1:s}').format(self.index, ostring)  # -1 needed to remove last tab...
            #self.panda_string += ostring + '\n'
            self.index += 1
        else:
            pass
    
    def write_string(self):
        """
        Write an output string using pandas dataframe
        """
        if self.outputMode == 'pandas':
            print('\nOUTPUTTING VIA PANDAS')
          #  self.colprint()
            df = pd.read_csv(pandas.compat.StringIO(self.panda_string), delimiter='\t')
           # print('Head write: \n', df.head(5), '\n')
            df.to_pickle(self.outFilename)
            print('Wrote pandas dataframe to pickled file: {0:s}'.format(self.outFilename))

    def get_file_information(self, dh=None):
        """
        get_file_information reads the sequence information from the
        currently selected data file

        Two-dimensional sequences are supported.
        :return nothing:
        """
        self.sequence = self.dataModel.listSequenceParams(dh)
        keys = self.sequence.keys()
        leftseq = [str(x) for x in self.sequence[keys[0]]]
        if len(keys) > 1:
            rightseq = [str(x) for x in self.sequence[keys[1]]]
        else:
            rightseq = []
        leftseq.insert(0, 'All')
        rightseq.insert(0, 'All')



    def file_cell_protocol(self, filename):
        """
        file_cell_protocol breaks the current filename down and returns a
        tuple: (date, cell, protocol)
        last argument returned is the rest of the path...
        """
        (p0, proto) = os.path.split(filename)
        (p1, cell) = os.path.split(p0)
        (p2, sliceid) = os.path.split(p1)
        (p3, date) = os.path.split(p2)
        return (date, sliceid, cell, proto, p3)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Generate Data Summaries from acq4 datasets')
    parser.add_argument('basedir', type=str,
                        help='Base Directory')
    parser.add_argument('-o', '--output', type=str, default='terminal', dest='output',
                        choices=['terminal', 'pandas', 'excel', 'tabfile'],
                        help='Specify output dataplan key for one entry to process')
    parser.add_argument('-f', '--filename', type=str, default='terminal', dest='outputFilename',
                        help='Specify output file name (including full path)')
    parser.add_argument('-r', '--read', action='store_true', dest='read',
                        help='just read the summary table')
    parser.add_argument('-w', '--write', action='store_true', dest='write',
                        help='janalyze and write the data summary')
    
    parser.add_argument('--daylist', type=str, default=None, dest='daylist',
                        help='Specify daylistfile')
    parser.add_argument('-a', '--after', type=str, default=None, dest='after',
                        help='only analyze dates on or after a date')
    parser.add_argument('-b', '--before', type=str, default=None, dest='before',
                        help='only analyze dates on or before a date')
    parser.add_argument('--dry-run', action='store_true', dest='dryrun',
                        help='Do a dry run, reporting only directories')
    parser.add_argument('--no-inspect', action='store_false', dest='noinspect',
                        help='Do not inspect protocols, only report directories')
    parser.add_argument('-d', '--depth', type=str, default='all', dest='depth',
                        choices = ['days', 'slices', 'cells', 'protocols', 'all'],
                        help='Specify depth for --dry-run')
    args = parser.parse_args()
    ds = DataSummary(basedir=args.basedir, daylistfile=args.daylist, outputMode=args.output, outputFile=args.outputFilename,
            after=args.after, before=args.before, dryrun=args.dryrun, depth=args.depth, inspect=args.noinspect)
 
    if args.write:
        ds.getDay()
        if args.output in ['pandas']:
            ds.write_string()
        
    if args.read:
        print('args.read')
        print('reading: ', ds.outFilename)
        df = pd.read_pickle(ds.outFilename)
        print(df.head(10))

        df2 = df.set_index("date", drop=False)
        print(df2.head(5))
        maps = []
        IVs = []
        stdIVs = []
        
        print(df2.columns.values)
        for day in range(len(df2.index)):
#       Only returns a dataframe if there is more than one entry
#       Otherwise, it is like a series or dict
            date = df2.iloc[day]['date']
            u = df2.iloc[day]['data_complete'].split(', ')
            prox = sorted(list(set(u)))  # adjust for duplicates (must be a better way in pandas)
            for p in prox:
                print('    protocol: ', p)

                c = date + '/' + df2.iloc[day]['slice_slice'] + '/' + df2.iloc[day]['cell_cell'] + '/' + p
                if 'Map' in c:
                    maps.append(c)
                if 'IV' in c:
                    IVs.append(c)
        # else:
       #      u = prots['data_complete'].split(', ')
       #      prox = sorted(list(set(u)))  # adjust for duplicates (must be a better way in pandas)
       #      for x in prox:  # construct filenames and sort by analysis types
       #          c = day + '/' + prots['slice_slice'] + '/' + prots['cell_cell'] + '/' + x
       #          if 'Map' in c:
       #              maps.append(c)
       #          if 'CCIV' in c in c:
       #              IVs.append(c)
       #          if 'CCIV_1nA_max' in c:
       #              stdIVs.append(c)
        print('maps: ', '   \n'.join(maps))
        print('ivs: ', '    \n'.join(IVs))
        print('std ivs: ', '    \n'.join(stdIVs))

    
