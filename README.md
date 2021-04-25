# ephysanalysis
This repository holds three interlocking sets of tools for analysis of voltage and current clamp data from electrophysiology 
experiments, laser scanning photostimulation and optogenetics. Includes analysis of IV curves, spike shape, 
rate adaptation, miniature syanptic events (multiple algorithms), and laser-scanning photostimulation maps.

Installation
------------

Clone the repository from github, and run "./make_env.sh". This will build a Python 3 environment with the
required dependencies. 

Tests
-----
A few tests are currently implemented, principally in the minianalysis section. Run "python tests.py" at the top level.

Issues
------
This is an evolving set of packages, and have been recently refactored. They are meant to provide support for data analysis and should
be called/imported from scripts that manage the incoming data files.



