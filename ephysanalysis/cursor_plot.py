import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph.Point import Point
from pyqtgraph.parametertree import Parameter, ParameterTree
font = QtGui.QFont()
font.setFamily('Arial')
font.setFixedPitch(False)
font.setPointSize(11)


#generate layout
class CursorPlot(object):
    def __init__(self):
        

        self.nextStore = 'T0'
        self.app = QtGui.QApplication([])
        self.win = pg.GraphicsWindow()
        self.layout = QtGui.QGridLayout()
        self.win.setLayout(self.layout)
        self.win.resize(1024,800)
        
        self.win.setWindowTitle('pyqtgraph example: crosshair')
        self.label = pg.LabelItem(justify='right')
        self.win.addItem(self.label)
        # self.plots['Windowed'] = self.win.addPlot(row=1, col=0)
        # self.plots['Full'] = self.win.addPlot(row=2, col=0)

        params = [
                {'name': 'Cursor Plot', 'type': 'group', 'children': [
                    {'name': 'Save T0', 'type': 'action'},
                    {'name': 'Save T1', 'type': 'action'},
                    {'name': 'Done', 'type': 'action'},
                ],
                },
                {'name': 'Data', 'type': 'group', 'children': [
                    {'name': 'NextStore', 'type': 'str', 'value': self.nextStore},
                    {'name': 'T0', 'type': 'float', 'value': 0.},
                    {'name': 'T1', 'type': 'float', 'value': 0.},
                ],
                }
            
            ]
        ptree = ParameterTree()
        ptree.setMaximumWidth(150)
        ptree.setFont(font)
        ptree.sizeHintForRow(24)
#        ptree.setUniformRowHeights(True)
        self.ptreedata = Parameter.create(name='params', type='group', children=params)
        ptree.setParameters(self.ptreedata, showTop=False)
        print(self.ptreedata.param)
        self.T0 = self.ptreedata.param('Data').param('T0')
        print(self.T0.value()) # 'T0')
        self.T1 = self.ptreedata.param('Data').param('T1')
        self.T0.sigValueChanged.connect(self.T0_Changed)
        self.T1.sigValueChanged.connect(self.T1_Changed)
        


        self.layout.addWidget(ptree, 0, 0, 9, 1) # Parameter Tree on left

        # add space for the graphs
        view = pg.GraphicsView()
        layout = pg.GraphicsLayout(border=(50,100,50))
        view.setCentralItem(layout)
        self.layout.addWidget(view, 1, 1, 8, 1)  # data plots on right

        self.plots = {}
        self.plots['Windowed'] = layout.addPlot(title="Windowed")
        layout.nextRow()
        self.plots['Full'] = layout.addPlot(title="Full")

        self.ptreedata.sigTreeStateChanged.connect(self.command_dispatcher)  # connect parameters to their updates

        self.region = pg.LinearRegionItem()
        self.region.setZValue(10)
        # Add the LinearRegionItem to the ViewBox, but tell the ViewBox to exclude this 
        # item when doing auto-range calculations.
        self.plots['Full'].addItem(self.region, ignoreBounds=True)

        #pg.dbg()
        self.plots['Windowed'].setAutoVisible(y=True)
        self.proxy = pg.SignalProxy(self.plots['Windowed'].scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
        self.proxyclick = pg.SignalProxy(self.plots['Windowed'].scene().sigMouseClicked, rateLimit=60, slot=self.mouseClicked)
    
    def command_dispatcher(self, param, changes):
        for param, change, data in changes:
            path = self.ptreedata.childPath(param)
            if path[1] == 'Done':
                # close window
                pass
            elif path[1] == 'Save T0':
                self.nextstore = 'T0'
            elif path[1] == 'Save T1':
                self.nextstore = 'T1'
            
            
            pass
            
    def T0_Changed(self):
        pass
    
    def T1_Changed(self):
        pass
    

    def make_testdata(self): #create numpy arrays
        #make the numbers large to show that the xrange shows data from 10000 to all the way 0
        dataY = 15000 + 15000 * pg.gaussianFilter(np.random.random(size=10000), 10) + 3000 * np.random.random(size=10000)
        # dataX = 10000 + 15000 * pg.gaussianFilter(np.random.random(size=10000), 10) + 3000 * np.random.random(size=10000)
        dataX = np.linspace(0, 100., dataY.shape[0])
        return(dataX, dataY)
        
    def plotData(self, x, y, setline=True):
        self.dataX = x
        self.dataY = y
        self.plots['Windowed'].plot(x=self.dataX, y=self.dataY, pen="g")
        # self.plots['Windowed'].plot(self.dataY, pen="g")
        self.plots['Full'].plot(x=self.dataX, y=self.dataY, pen="w")
        
        if setline:
            self.region.sigRegionChanged.connect(self.update)
            self.plots['Windowed'].sigRangeChanged.connect(self.updateRegion)
            xmin = np.min(self.dataX)
            xmax = np.max(self.dataX)
            xmin = 0.1*(xmax-xmin) + xmin
            xmax = 0.9*(xmax-xmin) + xmin
            self.region.setRegion([xmin, xmax])

            #cross hair
            self.vLine = pg.InfiniteLine(angle=90, movable=False)
            self.hLine = pg.InfiniteLine(angle=0, movable=False)
            self.plots['Windowed'].addItem(self.vLine, ignoreBounds=True)
            self.plots['Windowed'].addItem(self.hLine, ignoreBounds=True)
            self.vb = self.plots['Windowed'].vb

    def update(self):
        self.region.setZValue(10)
        minX, maxX = self.region.getRegion()
        self.plots['Windowed'].setXRange(minX, maxX, padding=0)    

    def updateRegion(self, window, viewRange):
        rgn = viewRange[0]
        self.region.setRegion(rgn)
        self.rgn = rgn

    def mouseMoved(self, evt):
        # evt is a PyQt5.QtCore.QPointF
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.plots['Windowed'].sceneBoundingRect().contains(pos):
            mousePoint = self.vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            if index > 0 and index < len(self.dataX):
                self.label.setText("<span style='font-size: 12pt'>x=%0.3f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), mousePoint.y(), self.dataY[index]))
            self.vLine.setPos(mousePoint.x())
            self.hLine.setPos(mousePoint.y())
            self.lastx = mousePoint.x()
            self.lasty = mousePoint.y()

    def mouseClicked(self, evt):
        # print(evt)
        # print(f"<{self.nextStore:s}>")
        # print(dir(self.T0))
        if self.nextStore == 'T0':
            self.T0.setValue(self.lastx)
            self.nextStore = 'T1'
        elif self.nextStore == 'T1':
            self.T1.setValue(self.lastx)
            self.nextStore = 'T0'
        else:
            pass
        

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    CP = CursorPlot()
    x, y = CP.make_testdata()
    CP.plotData(x, y)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()