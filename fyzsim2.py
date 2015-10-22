#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, random
from PyQt4 import QtGui, QtCore
#from math import sin
#from numpy import *
from time import time

from scipy.integrate import ode
import scipy

from scipy import zeros, sin, cos, floor

WIDTH = 600
HEIGHT = 400

CENTER = (WIDTH/2, HEIGHT/2)


# u_tt = ν*u_xx + r*u*(1 - u/κ)

# τ=0.0001 ~ 24fps

τ = 0.0005

FPSGOAL = 24
FPSFRAK = 1 / FPSGOAL

SCALE = 40

# L = .5
# K = 1
# M = 20
# mu = 1

# # 
# r, rd = 4, 0
# phi, phid = 0, 5**(1/2) / 4
# psi, psid = scipy.pi / 4, 0

L = 1
K = 1
M = 10
mu = 1

# 
r, rd = 4, -1
phi, phid = 0, .1
psi, psid = scipy.pi / 4, 1

alpha, alphad = phi - psi, phid - psid

y0, t0 = [r, rd, phi, phid, alpha], 0

beta = r*r * phid*phid - L*L*(alphad - phid)

t1 = 4
dt = .01


class Example(QtGui.QWidget):

    def __init__(self):
        super(Example, self).__init__()
        self.initUI()
        self.initModel()

    def initUI(self):
        self.move(300, 300)
        self.setFixedSize(WIDTH, HEIGHT) # no resizing
        self.setWindowTitle('PLanety')

        PENWIDTH = 1
        self.bluePen = QtGui.QPen(QtCore.Qt.blue, PENWIDTH, QtCore.Qt.SolidLine)
        self.redPen = QtGui.QPen(QtCore.Qt.red, PENWIDTH, QtCore.Qt.SolidLine)
        self.blackPen = QtGui.QPen(QtCore.Qt.black, PENWIDTH, QtCore.Qt.SolidLine)

        self.honza = QtCore.QTimer()
        self.honza.timeout.connect( self.timerTick)
        self.honza.start(floor(1000*FPSFRAK))
        self.show()
    
    def mousePressEvent(self, event):
        pass #self.goat_sim.mousePressEvent(event)

    def keyPressEvent(self, event):
        pass #self.goat_sim.keyPressEvent(event)

    def initModel(self):
        self.lasttime = 0
        self.fps = 0
        self.timecounter = 0

        self.odes = ode(prava).set_integrator('dopri5')
        self.odes.set_initial_value(y0, t0)

        self.x1 = r * cos(phi) + L * cos(psi)
        self.y1 = r * sin(phi) + L * sin(psi)
        self.x2 = r * cos(phi) - L * cos(psi)
        self.y2 = r * sin(phi) - L * sin(psi)
        self.E = 0


    def equation(self, τ):
        odes = self.odes
        if odes.successful():
            odes.integrate(odes.t + τ)
            t = odes.t
            r = odes.y[0]
            phi = odes.y[2]
            psi = odes.y[2] - odes.y[4]
            
            self.x1 = r * cos(phi) + L * cos(psi)
            self.y1 = r * sin(phi) + L * sin(psi)
            self.x2 = r * cos(phi) - L * cos(psi)
            self.y2 = r * sin(phi) - L * sin(psi)
            
            ca = cos(odes.y[4])
            sa = sin(odes.y[4])
            mphid = (odes.y[0]*odes.y[0]*odes.y[3] - beta) / (L*L)

            men1p = (odes.y[0] * odes.y[0] + L * L + 2 * odes.y[0] * L * ca) ** (1 / 2)
            men2p = (odes.y[0] * odes.y[0] + L * L - 2 * odes.y[0] * L * ca) ** (1 / 2)
            E = odes.y[1] * odes.y[1] + odes.y[0] * odes.y[0] * odes.y[3] * odes.y[3] + L * L * mphid*mphid - K*M * (1/men1p + 1/men2p)
            self.E = E * mu

    def timerTick(self):
        now = time()

        while self.timecounter < FPSFRAK:
            self.equation(τ) # NEXT STEP EVALUATION
            self.timecounter += τ

        # FRAMES PER SECOND
        self.timecounter -= FPSFRAK
        self.fps = floor(1 / (now - self.lasttime))
        self.lasttime = now
        self.update()

    def drawFPS(self, qp):
        qp.setPen(QtGui.QColor(0,0,0))
        qp.setFont(QtGui.QFont('Decorative', 10))
        qp.drawText(QtCore.QPoint(10,20), 'fps: %d; τ=%f; x1=%.2f, y1=%.2f, x2=%.2f, y2=%.2f, E=%f' % (self.fps, τ, self.x1, self.y1, self.x2, self.y2, self.E))
    
        
    def drawPoints(self, qp):
        POINTSIZE = 10
        qp.setBrush(QtGui.QColor(200,0,0))
        qp.drawEllipse(SCALE*self.x1 + CENTER[0], SCALE*self.y1 + CENTER[1], POINTSIZE, POINTSIZE)

        qp.setBrush(QtGui.QColor(0,200,0))
        qp.drawEllipse(SCALE*self.x2 + CENTER[0], SCALE*self.y2 + CENTER[1], POINTSIZE, POINTSIZE)

        # qp.setPen(self.redPen)
        # qp.drawPoint(SCALE*self.x1 + CENTER[0], SCALE*self.y1 + CENTER[1])
        # qp.setPen(self.bluePen)
        # qp.drawPoint(SCALE*self.x2 + CENTER[0], SCALE*self.y2 + CENTER[1])

        qp.setBrush(QtGui.QColor(200,200,200))
        qp.drawEllipse( CENTER[0], CENTER[1], POINTSIZE*2, POINTSIZE*2)

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        qp.begin(self)
        self.drawFPS(qp)
        self.drawPoints(qp)
        qp.end()

# r         = y[0]
# rd        = y[1]
# phi       = y[2]
# phid      = y[3]
# alpha     = y[4]
E = 0
def prava(t, y):
    global E
    ca = cos(y[4])
    sa = sin(y[4])

    mphid = (y[0]*y[0]*y[3] - beta) / (L*L)

    men1p = (y[0] * y[0] + L * L + 2 * y[0] * L * ca) ** (1 / 2)
    men2p = (y[0] * y[0] + L * L - 2 * y[0] * L * ca) ** (1 / 2)

    men1 = men1p**3
    men2 = men2p**3

    # return [
    #     # y' = \dot{y}
    #     y[1],
    #     # y'' = ...
    #     y[0] * y[3] * y[3] - K*M/4 * ((2*y[0] + 2*L*ca) / men1 + (2*y[0] - 2*L*ca) / men2),
    #     # phi' = \dot{phi}
    #     y[3],
    #     # phi'' = ...
    #     (-1 / (y[0] * y[0])) * (2*y[0]*y[1]*y[3] + K*M/4 * ( -2*y[0]*L*sa/men1 + 2*y[0]*L*sa/men2 )),
    #     # alpha' = ...
    #     (1/(L*L)) * ( y[3] * (y[0] * y[0] + L*L) - beta),
    # ]

    return [
        # y' = \dot{y}
        y[1],
        # y'' = ...
        y[0] * y[3] * y[3] - K*M/2 * ((y[0] + L*ca) / men1 + (y[0] - L*ca) / men2),
        # phi' = \dot{phi}
        y[3],
        # phi'' = ...
        (1 / y[0]) * (-2*y[1]*y[3] + K*M*L*sa/2 * ( 1/men1 - 1/men2)),
        # alpha' = ...
        (1/(L*L)) * ( y[3] * (y[0] * y[0] + L*L) - beta),
    ]

def main():
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()


