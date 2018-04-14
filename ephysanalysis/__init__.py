#!/usr/bin/env python

# Use Semantic Versioning, http://semver.org/
version_info = (0, 2, 0, 'a')
__version__ = "%d.%d.%d%s" % version_info

print "apparent version: ", __version__

import Fitting
import Utility

import acq4read
import MatdatacRead
import RmTauAnalysis
import SpikeAnalysis
import dataSummary

