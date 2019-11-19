"""
Fix the objective scale factor
Give "stated objective", "actual objective"

('objective', '4x 0.1na ACHROPLAN')

"""
import numpy as np
from pathlib import Path
from pyqtgraph import configfile
import pprint
import datetime
pp = pprint.PrettyPrinter(indent=4)
# retiga:

CineScale = 1.0
refscale = [(1./CineScale)*6.54e-6, -(1./CineScale)*6.54e-6]
print(refscale)
sfactors = [4., 10., 40., 63.]

objlist = {4: '4x 0.1na ACHROPLAN',
        10: '10x 0.3na W N-ACHROPLAN',
        40: '40x 0.8na ACHROPLAN',
        63: '63x 0.9na ACHROPLAN',
    }


bp = '/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107Ai32Het'
fn = '2019.11.15_000/slice_001/cell_001'

p = Path(bp, fn, '.index')

changelist = {'image_001.tif': [4, 10], 'image_002.tif': [4, 10], 'video_000.ma': [4, 10]}

def read_index(p):
    d = datetime.datetime.now()
    dstr = d.strftime('%Y-%m-%d %H:%M:%S')
    index = configfile.readConfigFile(p)
    print(list(index.keys()))
    for k in list(index.keys()):
        if k in list(changelist.keys()):
            print('index: ', k)
            print('old obj: ', index[k]['objective']) # pp.pprint(index[k] )
            oldobj = index[k]['objective']
            print('old trans: ', index[k]['transform'])
            print('old devtrans: ', index[k]['deviceTransform'])
            binning = index[k]['binning']
            newobj = changelist[k][1]
            fnewobj = float(newobj)
            index[k]['transform']['scale'] = (binning[0]*refscale[0]/fnewobj, binning[1]*refscale[1]/fnewobj, 1.0)
            index[k]['deviceTransform']['scale'] = (binning[0]*refscale[0]/fnewobj, binning[1]*refscale[1]/fnewobj, 1.0)
            index[k]['objective'] = objlist[newobj]
            index[k]['note'] = f'Objective scale corrected from {oldobj:s} to {objlist[newobj]:s} on {dstr:s} by PBM'
            
    for k in list(index.keys()):
        if k in list(changelist.keys()):
            print('newindex: ', k)
            print('newobj: ', index[k]['objective']) # pp.pprint(index[k] )
            print('newtrans: ', index[k]['transform'])
            print('newdevtrans: ', index[k]['deviceTransform'])
            print(index[k]['note'])
    
    return index

def write_index(p, index):
    configfile.writeConfigFile(index, p)

def main():
    index = read_index(p)
    write_index(p, index)

if __name__ == '__main__':
    main()




