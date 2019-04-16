import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph.Point import Point

#generate layout
class CursorPlot(object):
    def __init__(self):
        self.app = QtGui.QApplication([])
        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('pyqtgraph example: crosshair')
        self.label = pg.LabelItem(justify='right')
        self.win.addItem(self.label)
        self.p1 = self.win.addPlot(row=1, col=0)
        self.p2 = self.win.addPlot(row=2, col=0)

        self.region = pg.LinearRegionItem()
        self.region.setZValue(10)
        # Add the LinearRegionItem to the ViewBox, but tell the ViewBox to exclude this 
        # item when doing auto-range calculations.
        self.p2.addItem(self.region, ignoreBounds=True)

        #pg.dbg()
        self.p1.setAutoVisible(y=True)
        self.proxy = pg.SignalProxy(self.p1.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)

    def make_testdata(self): #create numpy arrays
        #make the numbers large to show that the xrange shows data from 10000 to all the way 0
        dataY = 15000 + 15000 * pg.gaussianFilter(np.random.random(size=10000), 10) + 3000 * np.random.random(size=10000)
        # dataX = 10000 + 15000 * pg.gaussianFilter(np.random.random(size=10000), 10) + 3000 * np.random.random(size=10000)
        dataX = np.linspace(0, 100., dataY.shape[0])
        return(dataX, dataY)
        
    def plotData(self, x, y, setline=True):
        self.dataX = x
        self.dataY = y
        self.p1.plot(x=self.dataX, y=self.dataY, pen="g")
        # self.p1.plot(self.dataY, pen="g")
        self.p2.plot(x=self.dataX, y=self.dataY, pen="w")
        
        if setline:
            self.region.sigRegionChanged.connect(self.update)
            self.p1.sigRangeChanged.connect(self.updateRegion)
            xmin = np.min(self.dataX)
            xmax = np.max(self.dataX)
            xmin = 0.1*(xmax-xmin) + xmin
            xmax = 0.9*(xmax-xmin) + xmin
            self.region.setRegion([xmin, xmax])

            #cross hair
            self.vLine = pg.InfiniteLine(angle=90, movable=False)
            self.hLine = pg.InfiniteLine(angle=0, movable=False)
            self.p1.addItem(self.vLine, ignoreBounds=True)
            self.p1.addItem(self.hLine, ignoreBounds=True)
            self.vb = self.p1.vb

    def update(self):
        self.region.setZValue(10)
        minX, maxX = self.region.getRegion()
        self.p1.setXRange(minX, maxX, padding=0)    

    def updateRegion(self, window, viewRange):
        rgn = viewRange[0]
        self.region.setRegion(rgn)
        self.rgn = rgn

    def mouseMoved(self, evt):
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.p1.sceneBoundingRect().contains(pos):
            mousePoint = self.vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            if index > 0 and index < len(self.dataX):
                self.label.setText("<span style='font-size: 12pt'>x=%0.3f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), mousePoint.y(), self.dataY[index]))
            self.vLine.setPos(mousePoint.x())
            self.hLine.setPos(mousePoint.y())
    #p1.scene().sigMouseMoved.connect(mouseMoved)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    CP = CursorPlotSetup()
    x, y = CP.make_testdata()
    CP.plotData(x, y)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()