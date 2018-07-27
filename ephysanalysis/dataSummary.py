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
Major surgery - to output Pandas (pickled) files as well. 

Future:
    Provide API/class interface (see all the "set" functions)

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

class Printer():
    """Print things to stdout on one line dynamically"""
    def __init__(self,data):
        sys.stdout.write("\033[1;36m\r\x1b[K\033[0;0m"+data.__str__())
        sys.stdout.flush()

class DataSummary():
    """
    Note that this init is just setup - you have to call getSummary on the object to do anything
    
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
    def __init__(self, basedir, outputMode='terminal', daylistfile=None,
                 after=None, before=None, dryrun=False, depth='all', inspect=True):


        self.setups()  # build some regex and wrappers
        # gather input parameters
        self.basedir = basedir
        self.outputMode = outputMode  # terminal, tabfile, pandas
        self.outFilename = ''
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
        
        # initialized dictionary that holds all the stuff
        self.analysis_summary = {}

        # column definitions - may need to adjust if change data that is pasted into the output
        self.coldefs = 'Date\tDescription\tNotes\tSpecies\tGenotype\tAge\tSex\tWeight\tTemp\tImportant\tElapsedT\tCellType\tSlice\tSliceNotes\t'
        self.coldefs += 'Cell\tCellNotes\tProtocols\tCompleteProtocols\tImages\n'
        self.panda_string = self.coldefs
        
        self.AR = acq4read.Acq4Read()  # instance of the reader

        outputDir = os.path.join(os.path.expanduser("~"), 'Desktop/acq4_scripts')
        if self.outputMode == 'tabfile':
            self.outFilename = self.basedir.replace('/', '_') + '.tab'
            self.outFilename = self.outFilename.replace('\\', '_')
            if self.outFilename[0] == '_':
                self.outFilename = self.outFilename[1:]
            self.outFilename = os.path.join(outputDir, self.outFilename)
            print('Tabfile output: Writing to {:<s}'.format(self.outFilename))
            h = open(self.outFilename, 'w')  # write new file
            h.write(self.basedir+'\n')
            h.write(self.coldefs + '\n')
            h.close()
        elif self.outputMode == 'pandas':  # save output as a pandas data structure, pickled
            self.outFilename = self.basedir.replace('/', '_') + '.pkl'
            self.outFilename = self.outFilename.replace('\\', '_')
            if self.outFilename[0] == '_':
                self.outFilename = self.outFilename[1:]
            self.outFilename = os.path.join(outputDir, self.outFilename)
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
            #maxdayx = datetime.datetime.now().timetuple()[0:3]  # get today
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
        
    def setMonitor(self):
        pass
    
    def setBaseDir(self):
        pass
    
    def setBegin(self):
        """
        begin (b) (define start date; end is set to current date) default: 1/1/1970
        """
        pass
    
    def setEnd(self):
        """
        end (e)(define end date: start date set to 1/1/1970; end is set to the end date) default: "today"
        """
        pass
    
    def setModeFull(self):
        """
        mode =  full (f) : do a full investigation of the data files. Makes processing very slow. (reports incomplete protocols)
        """
        pass
    
    def setModePartial(self):
        """
        partial (p) : do a partial investiagion of protocols: is there anything in every protocol directory? (reports incomplete protocols) - slow
        """
        pass
    
    def setModeQuick(self):
        """
        quick (q) : do a quick scan : does not run through protocols to find incomplete protocols. Default (over full and partial)
        """
        pass
    
    def setDebug(self):
        """
        debug (d) : debug monitoring of progress
        """
        pass
        
    def setOutput(self):
        """
        output (o) : define output file (tab delimited file for import to other programs)
        """
        pass
    
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

    def getSummary(self):
        """
        getSummary is the entry point for scanning through all the data files in a given directory,
        returning information about those within the date range, with details as specified by the options
        """
        allfiles = os.listdir(self.basedir)
        self.pstring = ''
        days = []
        for thisfile in allfiles:
            m = self.daytype.match(thisfile)
            if m == '.DS_Store':
                continue
            if m is None:
               # print 'Top level file %s is incorrectly placed ' % thisfile
                continue  # no match
            if len(m.groups()) >= 3:  # perfect match
                # print m.groups()
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
#            if self.monitor:
            self.pstring = 'Processing day[%3d/%3d]: %s ' % (nd, len(days), day)
            Printer(self.pstring)
            self.daystring = '%s\t' % (day.strip())
            self.AR.setProtocol(os.path.join(self.basedir, day))
            self.day = self.AR.getIndex()
            self.daystring = self.wrapstring(self.daystring, self.day, 'description', self.twd['day'])
            self.daystring = self.wrapstring(self.daystring, self.day, 'notes', self.tw['day'])
            if self.dryrun and self.depth=='days':
                self.outputString(self.daystring + self.day_summarystring)                
            else:
                self.doSlices(os.path.join(self.basedir, day))
            os.closerange(8, 65535)  # close all files in each iteration
            gc.collect()

    def doSlices(self, day):
        """
        process all of the slices for a given day
        :param day:
        :return nothing:
        """

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
            self.slice = self.AR.getIndex(os.path.join(day, slicen))
            self.slicestring = self.wrapstring(self.slicestring, self.slice, 'notes', self.tw['slice'])
            if self.dryrun and self.depth=='slices':
                self.outputString(self.daystring + self.day_summarystring + self.slicestring)
            else:
                self.doCells(os.path.join(day, slicen))
            gc.collect()

    def doCells(self, slice):
        """
        process all of the cells from a slice
        :param slice:
        :return nothing:
        """
        allfiles = os.listdir(slice)
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
            self.cellstring = '%s\t' % (cell)
            try:
                self.cell = self.AR.getIndex(os.path.join(slice, cell))  # possible that .index file is missing, so we cannot read
            except:
                self.cell = None
            self.cellstring = self.wrapstring(self.cellstring, self.cell, 'notes', self.tw['cell'])
            self.day_summary(name=os.path.join(slice, cell))
            if self.dryrun and self.depth=='cells':
                self.outputString(self.daystring + self.day_summarystring + self.slicestring + self.cellstring)
            else:
                self.doProtocols(os.path.join(slice, cell))
            gc.collect()

    def doProtocols(self, cell):
        """
        process all of the protocols for a given cell
        :param cell:
        :return nothing:
        """
        #print( 'doing protocols')
        allfiles = os.listdir(cell)
        #celltype = re.compile("(Cell_)(\d{3,3})")
        protocols = []
        nonprotocols = []
        anyprotocols = False
        images = []  # tiff
        stacks2p = []
        images2p = []
        videos = []

        endmatch = re.compile("[\_(\d{3,3})]$")  # look for _lmn at end of directory name
        for thisfile in allfiles:
            if os.path.isdir(os.path.join(cell, thisfile)):
                protocols.append(thisfile)
            else:
                nonprotocols.append(thisfile)

        self.protocolstring = ''
        self.allprotocols = ''
        self.completeprotocols = ''
        if self.InvestigateProtocols is True:
           # self.day_summarystring = 'NaN\t'*6

            for np, protocol in enumerate(protocols):
                Printer(self.cstring + ' Prot[%2d/%2d]: %s' % (np, len(protocols), protocol))
                self.allprotocols += protocol+', '
                protocolpath = os.path.join(cell, protocol)
                if np == 0:
                    self.day_summary(name=protocolpath)
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
                        ncomplete = ncomplete + 1
                        self.completeprotocols += protocol+', '
                
                gc.collect()  # and force garbage collection of freed objects inside the loop
                if modes == []:
                    modes = ['Unknown mode']
                if protocolok and ncomplete == nexpected:  # accumulate protocols
                    self.protocolstring += '[{:<s}: {:s} {:d}], '.format(protocol, modes[0][0], ncomplete)
                    anyprotocols = True  # indicate that ANY protocol ran to completion
                else:
                    if self.reportIncompleteProtocols:
                        self.protocolstring += '[{:<s}, ({:s}, {:d}/{:d}, Incomplete)], '.format(protocol, modes[0][0], ncomplete, nexpected)
                gc.collect()
        else:
            self.protocolstring += 'Protocols: '

            anyprotocols = True
            prots = {}
            for np, protocol in enumerate(protocols):
                Printer(self.cstring + ' Prot[%2d/%2d]: %s' % (np,len(protocols), protocol))
                self.allprotocols += protocol+', '
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
                self.protocolstring += '['
                for p in prots.keys():
                    self.protocolstring += '{:<s}({:<d}), '.format(p, prots[p])
                self.protocolstring += ']'
            else:
                self.protocolstring = '<No protocols found>'
        self.protocolstring += '\t'
        self.allprotocols += '\t'
        self.completeprotocols += '\t'

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
        
        if anyprotocols:
            ostring = self.daystring + self.day_summarystring + self.slicestring + self.cellstring + self.protocolstring + self.completeprotocols + self.imagestring
        else:
            ostring = self.daystring + self.day_summarystring + self.slicestring + self.cellstring + self.protocolstring + self.completeprotocols + self.imagestring
        self.outputString(ostring)

    def outputString(self, ostring):
        if self.outputMode in [None, 'terminal']:
            print('{0:s}'.format(ostring))
        elif self.outputMode == 'text':
            h = open(self.outFilename, 'a')  # append mode
            h.write(ostring + '\n')
            h.close()
        elif self.outputMode == 'pandas':
            #print('\n******* building pd string', 'ostring: \n', ostring)
            self.panda_string += ('{0:d}\t').format(self.index) + ostring+'\n'  # -1 needed to remove last tab... 
            # this just is useful to check that the data are in the right tabbed fields
            # print('\n', self.panda_string)
            # tstr = ( ostring+'\n').split('\t')
            # for i, c in enumerate(self.coldefs.split('\t')):
            #     #print (c)
            #     print(c, ' : ', tstr[i])
            #
            # exit(1)
            self.index += 1
        else:
            pass
    
    def write_string(self):
        """
        Write an output string using pandas dataframe
        """
        if self.outputMode == 'pandas':
            df = pd.read_csv(pandas.compat.StringIO(self.panda_string), delimiter='\t')
            df.to_pickle(self.outFilename)
            # print('Wrote pandas dataframe to pickled file: {0:s}'.format(self.outFilename))

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

    def day_summary(self, name=None):
        """
        day_summary generates a dictionary of information about the cell
        for the selected directory handle (usually a protocol; could be a file)
        builds a formatted string with some of the information.
        
        :return nothing:
        """
        self.analysis_summary = {}  # always clear the summary.
        self.day_summarystring = ''
        # other info into a dictionary
        self.analysis_summary['Day'] = self.day
        self.analysis_summary['Slice'] = self.slice
        self.analysis_summary['Cell'] = self.cell
        self.analysis_summary['ACSF'] = self.day['solution']
        self.analysis_summary['Internal'] = self.day['internal']
        self.analysis_summary['Temp'] = self.day['temperature']
        self.analysis_summary['Important'] = 'N'
        print('\nself cell: ', self.cell)
        print('analysis summary[cell]: ', self.analysis_summary['Cell'])
        if 'important' in self.day.keys():
#            print('important: ', self.day['important'])
            if self.day['important']:
                self.analysis_summary['Important'] = 'Y'
                self.day['important'] = 'N'
        else:
            self.day['important'] = 'N'
        self.analysis_summary['CellType'] = 'X'
        if self.cell is not None and 'type' in self.cell.keys():
            print ('\n   cell type: ', self.cell['type'], len(self.cell['type']))
            if len(self.cell['type']) > 0:  # only set if there is a value
                self.analysis_summary['CellType'] = self.cell['type']
            if isinstance(self.cell['type'], str):
                print('cell is string, len is: %d', len(self.cell['type']))
            else:
                print('isnan cell type: ', math.isnan(self.cell['type']))
        
        today = self.analysis_summary['Day']
        if today is not None:

            if 'description' in today.keys():
                self.analysis_summary['Description'] = today['description']
            if 'weight' in today.keys():
                self.analysis_summary['Weight'] = today['weight']
            if 'strain' in today.keys():
                self.analysis_summary['Strain'] = today['strain']
            if 'notes' in today.keys():
                self.analysis_summary['Notes'] = today['notes']
            if 'age' in today.keys():
                self.analysis_summary['Age'] = today['age']
            if 'sex' in today.keys():
                self.analysis_summary['Sex'] = today['sex']
            if 'animal identifier' in today.keys():
                self.analysis_summary['AnimalIdentifier'] = today['animal identifier']
            if 'genotype' in today.keys():
                self.analysis_summary['Genotype'] = today['genotype']
            if 'temperature' in today.keys():
                self.analysis_summary['Temp'] = today['temperature']
            if 'species' in today.keys():
                self.analysis_summary['Species'] = today['species']


        if self.analysis_summary['Cell'] is not None:
            ct = self.analysis_summary['Cell']['__timestamp__']
        else:
            ct = 0.
        try:
            pt = dh.info()['__timestamp__']
        except:
            pt = 0.
        self.analysis_summary['ElapsedTime'] = pt-ct  # save elapsed time between cell opening and protocol start
        (date, sliceid, cell, proto, p3) = self.file_cell_protocol(name)
        self.analysis_summary['CellID'] = os.path.join(date, sliceid, cell)  # use this as the "ID" for the cell later on
        data_template = (
            OrderedDict([('Species', '{:>s}'), ('Genotype', '{:>12s}'), ('Age', '{:>5s}'), ('Sex', '{:>1s}'), ('Weight', '{:>5s}'),
                         ('Temp', '{:>5s}'), ('Important', '{:>1s}'), ('ElapsedTime', '{:>8.2f}'), ('CellType', '{:>12s}')])
                        )
        ltxt = ''
        utxt = ''  # to make a little table comparing the template and data retrieved
        for a in data_template.keys():
            if a in self.analysis_summary.keys():
                ltxt += ((data_template[a] + '\t').format(self.analysis_summary[a]))
                utxt += '{0:s} = {1:s}\n'.format(a, str(self.analysis_summary[a]))
            else:
                ltxt += (('NaN\t'))
                utxt += '{0:s} = {1:s}\n'.format(a, 'missing in keys')
        self.day_summarystring = ltxt
        
       # print('cell summary: \n', utxt)

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
    parser.add_argument('-r', '--read', action='store_true', dest='read',
                        help='just read the protocol')
    args = parser.parse_args()
    ds = DataSummary(basedir=args.basedir, daylistfile=args.daylist, outputMode=args.output,
            after=args.after, before=args.before, dryrun=args.dryrun, depth=args.depth, inspect=args.noinspect)
    if not args.read:
        ds.getSummary()
        if args.output in ['pandas']:
            ds.write_string()
        
        # now read the same file and look for some protocols
        df = pd.read_pickle(ds.outFilename)
        # print(df.head(5))
        # print('------')
        # print('Known column names: ', df.columns.values)
        # for c in [ 'Date', 'Description', 'Notes', 'Species', 'Genotype', 'Age', 'Sex', 'Weight',
        #     'Temp', 'Important', 'ElapsedT', 'CellType', 'Slice', 'SliceNotes' ,'Cell',
        #     'CellNotes', 'Protocols', 'AllProtocols', 'Images',]:
        #    print('****  {0:s} ****'.format(c))
        #    print(df[c])
        #allp = df.CompleteProtocols
        df.set_index("Date", inplace=True)
        # print("protocols day: ", df.loc['2017.10.25_000']['AllProtocols'])
        allp = df.loc[['2017.10.25_000'],['Important', 'Slice', 'Cell','CompleteProtocols']]
        print(allp)
    if args.read:
        df = pd.read_pickle(ds.outFilename)
        print(df.head(10))

        df2 = df.set_index("Date", drop=False)
        #print(df2.head(5))
        # select a date:
        
        day= '2017.07.28_000'
        
        # get the complete protocols:
        prots = df2.loc[day, ['Slice', 'Cell', 'CompleteProtocols']]  # to build datapaths 
        
        maps = []
        IVs = []
        stdIVs = []
#        print('isdataframe: ', isinstance(prots, pd.DataFrame))
#       Only returns a dataframe if there is more than one entry
#       Otherwise, it is like a series or dict
        if isinstance(prots, pd.DataFrame):
            for i in range(prots.shape[0]):
                u = prots.iloc[i]['CompleteProtocols'].split(', ')  # get the protocols
                prox = sorted(list(set(u)))  # adjust for duplicates (must be a better way in pandas)
                for x in prox:  # construct filenames and sort by analysis types
                    c = day + '/' + prots.iloc[i]['Slice'] + '/' + prots.iloc[i]['Cell'] + '/' + x
                    if 'Map' in c:
                        maps.append(c)
                    if 'IV' in c:
                        IVs.append(c)
#                    print (c)
        else:
            u = prots['CompleteProtocols'].split(', ')
            prox = sorted(list(set(u)))  # adjust for duplicates (must be a better way in pandas)
            for x in prox:  # construct filenames and sort by analysis types
                c = day + '/' + prots['Slice'] + '/' + prots['Cell'] + '/' + x
                if 'Map' in c:
                    maps.append(c)
                if 'CCIV' in c in c:
                    IVs.append(c)
                if 'CCIV_1nA_max' in c:
                    stdIVs.append(c)
#                print (c)
#            print('u: ', u)
        print('maps: ', '   \n'.join(maps))
        print('ivs: ', '    \n'.join(IVs))
        print('std ivs: ', '    \n'.join(stdIVs))
        #print(allp)
        # prots = allp['CompleteProtocols']
 #        print(prots)
        #
        
    
    
