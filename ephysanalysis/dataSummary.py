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

Future:
    Provide API/class interface (see all the "set" functions)

"""
import sys
import os
import re
import os.path
import gc
import argparse
import datetime
import numpy as np
import textwrap
from collections import OrderedDict

from ephysanalysis import acq4read
from pyqtgraph.metaarray import MetaArray



class DataSummary():
    def __init__(self, basedir=None, daylistfile=None, dayspec=None):
        self.monitor = False
        self.reportIncompleteProtocols = True  # do include incomplete protocol runs in print
        self.InvestigateProtocols = True  # set True to check out the protocols in detail
        
        self.analysis_summary = {}
        self.AR = acq4read.Acq4Read()  # instance of the reader
#        self.dataModel = PatchEPhys
        self.basedir = basedir
        self.dayspec = dayspec

        # column definitions - may need to adjust if change data that is pasted into the output
        self.coldefs = 'Date\tDescription\tNotes\tGenotype\tAge\tSex\tWeight\tTemp\tElapsed T\tSlice\tSlice Notes\t'
        self.coldefs += 'Cell\t Cell Notes\t\tProtocols\tImages\t'

        self.outputMode = 'terminal' # 'tabfile'
        outputDir = os.path.join(os.path.expanduser("~"), 'Desktop/acq4_scripts')
        if self.outputMode == 'tabfile':
            self.outFilename = basedir.replace('/', '_') + '.tab'
            self.outFilename = self.outFilename.replace('\\', '_')
            if self.outFilename[0] == '_':
                self.outFilename = self.outFilename[1:]
            self.outFilename = os.path.join(outputDir, self.outFilename)
            print('Writing to: {:<s}'.format(self.outFilename))
            h = open(self.outFilename, 'w')  # write new file
            h.write(basedir+'\n')
            h.write(self.coldefs + '\n')
            h.close()
        else:
            print ('Base Directory: ', basedir)
            print (self.coldefs)
            
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


        
        # look for names that match the acq4 "day" template:
        # example: 2013.03.28_000
        self.daytype = re.compile("(\d{4,4}).(\d{2,2}).(\d{2,2})_(\d{3,3})")
#        daytype = re.compile("(2011).(06).(08)_(\d{3,3})")  # specify a day
        #2011.10.17_000
        # operate in two modes:
        # second, between two dates
        self.daylist = None
        self.daylistfile = daylistfile
        if daylistfile is None:

            mindayx = (1970, 1, 1)
            #mindayx = (2018, 1, 26)
            self.minday = mindayx[0]*1e4+mindayx[1]*1e2+mindayx[2]
            #maxdayx = datetime.datetime.now().timetuple()[0:3]  # get today
            maxdayx = (2019, 1, 26)
            self.maxday = maxdayx[0]*1e4+maxdayx[1]*1e2+maxdayx[2]
        else:
            self.daylist = []
            with open(self.daylistfile, 'r') as f:
                for line in f:
                    if line[0] != '#':
                        self.daylist.append(line[0:10])
            f.close()

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
        
    def getSummary(self):
        """
        getSummary is the entry point for scanning through all the data files in a given directory,
        returning information about those within the date range, with details as specified by the options
        """
        allfiles = os.listdir(self.basedir)
        
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
        for day in days:
            if self.monitor:
                print ('processing day: %s' % day)
            self.daystring = '%s \t' % (day)
            self.AR.setProtocol(os.path.join(self.basedir, day))
            self.day = self.AR.getIndex()
            if self.day is not None and 'description' in self.day.keys() and len(self.day['description']) > 0:
                l = self.twd['day'].wrap(self.day['description'])
                for i in l:
                    i = i.replace('\t', '    ')  # clean out tabs so printed formatting is not confused
                    self.daystring += i
            else:
                self.daystring += ' [no description]'
            self.daystring += '\t'
            if self.day is not None and 'notes' in self.day.keys() and len(self.day['notes']) > 0:
                l = self.tw['day'].wrap(self.day['notes'])
                for i in l:
                    i = i.replace('\t', '    ')  # clean out tabs so printed formatting is not confused
                    self.daystring += i
            else:
                self.daystring += ' [no notes]'
            self.daystring += '\t'
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
        for slice in slices:
            self.slicestring = '%s\t' % (slice)
            self.slice = self.AR.getIndex(os.path.join(day, slice))
            if self.slice is not None and 'notes' in self.slice.keys() and len(self.slice['notes']) > 0:
                l = self.tw['slice'].wrap(self.slice['notes'])
                for i in l:
                    i = i.replace('\t', '    ')  # clean out tabs so printed formatting is not confused
                    self.slicestring += i
            else:
                self.slicestring += ' No slice notes'
            self.slicestring += '\t'
            self.doCells(os.path.join(day, slice))
            gc.collect()

    def doCells(self, slice):
        """
        process all of the cells from a slice
        :param slice:
        :return nothing:
        """
        allfiles = os.listdir(slice)
        celltype = re.compile("(cell_)(\d{3,3})")
        cells = []
        for thisfile in allfiles:
            m = celltype.match(thisfile)
            if m is None:
                continue
            if len(m.groups()) == 2:
                cells.append(thisfile)
        for cell in cells:
            self.cellstring = '%s\t' % (cell)
            try:
                self.cell = self.AR.getIndex(os.path.join(slice, cell))  # possible that .index file is missing, so we cannot read
            except:
                self.cell = None

            if self.cell is not None and 'notes' in self.cell.keys() and len(self.cell['notes']) > 0:
                l = self.tw['cell'].wrap(self.cell['notes'])
                for i in l:
                    i = i.replace('\t', '    ')  # clean out tabs so printed formatting is not confused
                    self.cellstring += i
            else:
                self.cellstring += ' No cell notes'
            self.cellstring += '\t'
            self.cell_summary(name=os.path.join(slice, cell))
            self.cellstring += '\t'
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
        if self.InvestigateProtocols is True:
            self.summarystring = 'NaN \t'*6

            for np, protocol in enumerate(protocols):

                protocolpath = os.path.join(cell, protocol)
                if np == 0:
                    self.cell_summary(name=protocolpath)
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
            for protocol in protocols:
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

        for thisfile in nonprotocols:
#            if os.path.isdir(os.path.join(cell, thisfile)):  # skip protocols
#                continue
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
            self.imagestring += 'Images: %d ' % len(images)
        if len(stacks2p) > 0:
            self.imagestring += '2pStacks: %d ' % len(stacks2p)
        if len(images2p) > 0:
            self.imagestring += '2pImages: %d ' % len(images2p)
        if len(videos) > 0:
            self.imagestring += 'Videos: %d' % len(videos)
        if len(images) + len(stacks2p) + len(images2p) + len(videos) == 0:
            self.imagestring = 'No Images or Videos'
        
        if anyprotocols:
            ostring =  self.daystring + self.summarystring + self.slicestring + self.cellstring + self.protocolstring + self.imagestring + '\t'
        else:
            ostring = self.daystring + self.summarystring + self.slicestring + self.cellstring + '<No complete protocols> \t' + self.imagestring + '\t'
        self.outputString(ostring)

    def outputString(self, ostring):
        if self.outputMode == 'terminal':
            print (ostring)
        else:
            h = open(self.outFilename, 'a')  # append mode
            h.write(ostring + '\n')
            h.close()
        
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

    def cell_summary(self, name=None):
        """
        cell_summary generates a dictionary of information about the cell
        for the selected directory handle (usually a protocol; could be a file)
        builds a formatted string with some of the information.
        :param dh: the directory handle for the data, as passed to loadFileRequested
        :return nothing:
        """
        self.analysis_summary = {}  # always clear the summary.
        self.summarystring = ''
        # other info into a dictionary
        self.analysis_summary['Day'] = self.day
        self.analysis_summary['Slice'] = self.slice
        self.analysis_summary['Cell'] = self.cell
        self.analysis_summary['ACSF'] = self.day['solution']
        self.analysis_summary['Internal'] = self.day['internal']
        self.analysis_summary['Temp'] = self.day['temperature']
        if self.cell is not None and 'type' in self.cell.keys():
            self.analysis_summary['CellType'] = self.cell[u'type']
        else:
            self.analysis_summary['CellType'] = 'N/A'
        today = self.analysis_summary['Day']
        if today is not None:
            #print today.keys()
            if 'species' in today.keys():
                self.analysis_summary['Species'] = today['species']
            if 'age' in today.keys():
                self.analysis_summary['Age'] = today['age']
            if 'sex' in today.keys():
                self.analysis_summary['Sex'] = today['sex']
            if 'weight' in today.keys():
                self.analysis_summary['Weight'] = today['weight']
            if 'temperature' in today.keys():
                self.analysis_summary['Temperature'] = today['temperature']
            if 'description' in today.keys():
                self.analysis_summary['Description'] = today['description']
            if 'notes' in today.keys():
                self.analysis_summary['Notes'] = today['notes']
        
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
            OrderedDict([('Species', '{:>s}'), ('Age', '{:>5s}'), ('Sex', '{:>1s}'), ('Weight', '{:>5s}'),
                         ('Temperature', '{:>5s}'), ('ElapsedTime', '{:>8.2f}')]))
 
        ltxt = ''
        for a in data_template.keys():
            if a in self.analysis_summary.keys():
                ltxt += ((data_template[a] + '\t').format(self.analysis_summary[a]))
            else:
                ltxt += (('NaN \t'))
        self.summarystring = ltxt

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
    if len(sys.argv) == 2:
        ds = DataSummary(basedir=sys.argv[1])
    if len(sys.argv) == 3:
        ds = DataSummary(basedir=sys.argv[1], daylistfile=sys.argv[2])
    ds.getSummary()
