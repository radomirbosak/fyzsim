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

τ = 0.001

FPSGOAL = 24
FPSFRAK = 1 / FPSGOAL


def paramset_old():
    global L, K, M, mu
    global r, rd, phi, phid, psi, psid, L, Ld
    global Lmin, Lmax
    global tuhost, trenie
    global TIMEFACTOR, SCALE
    global τ

    τ = 0.001

    TIMEFACTOR, SCALE = 10, 40
    L, K, M, mu = 1, 1, 1, 1

    r, rd = 5, 0
    phi, phid = 0, 1 / r**(3/2)
    psi, psid = scipy.pi / 4, 0
    L, Ld = 1, 0

    Lmin, Lmax = 0.9, 1.1

    tuhost = .025
    trenie = 1

def paramset_moon():
    global L, K, M, mu
    global r, rd, phi, phid, psi, psid, L, Ld
    global Lmin, Lmax
    global tuhost, trenie
    global TIMEFACTOR, SCALE
    global τ

    τ = 0.01

    TIMEFACTOR, SCALE = 1600, 2
    L, K, M, mu = 1, 1, 1, 1

    r, rd = 100, 0
    phi, phid = 0, 1 / r**(3/2)
    psi, psid = scipy.pi / 4, 0
    L, Ld = 1, 0

    Lmin, Lmax = 0.9, 1.1

    tuhost = .025
    trenie = 1.0

#paramset_moon()
paramset_old()

alpha, alphad = phi - psi, phid - psid

y0, t0 = [r, rd, phi, phid, alpha, L, Ld], 0

beta = r*r * phid - L*L*(alphad - phid)

class Example(QtGui.QWidget):

    def __init__(self):
        super(Example, self).__init__()
        self.initUI()
        self.initModel()

    def initUI(self):
        self.move(300, 300)
        self.setFixedSize(WIDTH, HEIGHT) # no resizing
        self.setWindowTitle('Planety s pružinou')

        PENWIDTH = 1
        self.bluePen = QtGui.QPen(QtCore.Qt.blue, PENWIDTH, QtCore.Qt.SolidLine)
        self.redPen = QtGui.QPen(QtCore.Qt.red, PENWIDTH, QtCore.Qt.SolidLine)
        self.blackPen = QtGui.QPen(QtCore.Qt.black, PENWIDTH, QtCore.Qt.SolidLine)

        self.honza = QtCore.QTimer()
        self.honza.timeout.connect( self.timerTick)
        self.honza.start(floor(1000*FPSFRAK))
        self.show()
    
    def mousePressEvent(self, event):
        pass

    def keyPressEvent(self, event):
        pass 

    def initModel(self):
        self.lasttime = 0
        self.fps = 0
        self.timecounter = 0
        self.model_total_time = 0

        self.odes = ode(prava).set_integrator('dopri5', nsteps=10000)
        self.odes.set_initial_value(y0, t0)

        self.x1 = r * cos(phi) + L * cos(psi)
        self.y1 = r * sin(phi) + L * sin(psi)
        self.x2 = r * cos(phi) - L * cos(psi)
        self.y2 = r * sin(phi) - L * sin(psi)
        self.E = 0
        self.L = 1
        self.phid, self.alphad = 1, .5

        self.HISTORY_MAX = 700
        self.HISTORY_SKIP = 100
        self.history_phase = 0
        self.history_tick = 0
        self.history = [(0,0,0)] * self.HISTORY_MAX

        self.measured_period = 0


    def equation(self, τ):
        odes = self.odes
        if odes.successful():
            odes.integrate(odes.t + τ)
            t = odes.t
            r = odes.y[0]
            phi = odes.y[2]
            psi = odes.y[2] - odes.y[4]
            L = odes.y[5]
            Ld = odes.y[6]
            self.phid = odes.y[3]
            self.alphad = 1 / (L*L) * (self.phid * (r*r + L*L) - beta)
            
            self.x1 = r * cos(phi) + L * cos(psi)
            self.y1 = r * sin(phi) + L * sin(psi)
            self.x2 = r * cos(phi) - L * cos(psi)
            self.y2 = r * sin(phi) - L * sin(psi)


            if self.history_tick >= self.HISTORY_SKIP:
                self.history[self.history_phase] = (r * cos(phi), r * sin(phi), phi)
                self.history_phase = (self.history_phase + 1) % self.HISTORY_MAX
                self.history_tick -= self.HISTORY_SKIP
                self.measured_period = self.find_period()
                # najdi periodu

            self.history_tick += TIMEFACTOR

            self.L = L
            
            ca = cos(odes.y[4])
            sa = sin(odes.y[4])
            mphid = (odes.y[0]*odes.y[0]*odes.y[3] - beta) / (L*L)

            men1p = (odes.y[0] * odes.y[0] + L * L + 2 * odes.y[0] * L * ca) ** (1 / 2)
            men2p = (odes.y[0] * odes.y[0] + L * L - 2 * odes.y[0] * L * ca) ** (1 / 2)
            E = odes.y[1] * odes.y[1] + odes.y[0] * odes.y[0] * odes.y[3] * odes.y[3] + L * L * mphid*mphid + Ld*Ld - K*M * (1/men1p + 1/men2p) + util(L)
            self.E = E * mu
        else:
            print("pokazilo sa dačo")
            exit(1)

    def get_psid(self):
        return self.phid - self.alphad

    def find_period(self):
        _, _, initial = self.history[(self.history_phase - 1) % self.HISTORY_MAX]
        for it in range(self.HISTORY_MAX):
            _, _, new = self.history[(self.history_phase - it - 1) % self.HISTORY_MAX]
            if abs(new - initial) > 2 * scipy.pi:
                break
        return (it + 1) * τ * self.HISTORY_SKIP

    def timerTick(self):
        now = time()

        while self.timecounter < FPSFRAK:
            self.equation(τ*TIMEFACTOR) # NEXT STEP EVALUATION
            self.timecounter += τ
            self.model_total_time += τ*TIMEFACTOR

        # FRAMES PER SECOND
        self.timecounter -= FPSFRAK
        self.fps = floor(1 / (now - self.lasttime))
        self.lasttime = now
        self.update()

    def drawFPS(self, qp):
        qp.setPen(QtGui.QColor(0,0,0))
        qp.setFont(QtGui.QFont('Decorative', 10))
        one = self.measured_period / (2 * scipy.pi / self.get_psid())
        txt = 'fps: {:3}; τ={:<6.4f}; L={:6.4f}, E={:9.6f}, fr={:6.4f}, T={:6.4f}'
        txt2 = 'T/(2pi/phid)={:6.4f}'
        qp.drawText(QtCore.QPoint(10,20), txt.format(self.fps, τ, self.L, self.E, trenie, self.measured_period))
        qp.drawText(QtCore.QPoint(10,40), txt2.format(one))
    
        
    def drawPoints(self, qp):
        POINTSIZE = 10
        POINTDIF = POINTSIZE / 2

        qp.setPen(self.blackPen)
        qp.drawLine(SCALE*self.x1 + CENTER[0] + 0, SCALE*self.y1 + CENTER[1] + 0, SCALE*self.x2 + CENTER[0] + 0, SCALE*self.y2 + CENTER[1] + 0)

        qp.setBrush(QtGui.QColor(200,0,0))
        qp.drawEllipse(SCALE*self.x1 + CENTER[0] - POINTDIF, SCALE*self.y1 + CENTER[1] - POINTDIF, POINTSIZE, POINTSIZE)

        qp.setBrush(QtGui.QColor(0,200,0))
        qp.drawEllipse(SCALE*self.x2 + CENTER[0] - POINTDIF, SCALE*self.y2 + CENTER[1] - POINTDIF, POINTSIZE, POINTSIZE)

        qp.setBrush(QtGui.QColor(200,200,200))
        qp.drawEllipse( CENTER[0]- POINTDIF*2, CENTER[1]- POINTDIF*2, POINTSIZE*2, POINTSIZE*2)

    def drawTail(self, qp):
        lx,ly,lphi = self.history[self.history_phase - 1]
        for hi in range(self.HISTORY_MAX - 1):
            x,y,phi = self.history[(self.history_phase - hi - 2) % self.HISTORY_MAX]

            hot = int(255 * hi / self.HISTORY_MAX)
            qp.setPen(QtGui.QColor(255, hot, hot))
            qp.drawLine(CENTER[0] + SCALE*lx, CENTER[1] + SCALE*ly, CENTER[0] + SCALE*x, CENTER[1] + SCALE*y)
            lx, ly, lphi = x, y, phi

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        qp.begin(self)
        self.drawFPS(qp)
        self.drawTail(qp)
        self.drawPoints(qp)
        qp.end()

def util(L):
    return tuhost * (1 / (L - Lmin) + 1 / (Lmax - L))

def util_der(L):
    return tuhost * (-1 / (L - Lmin)**2 + 1 / (Lmax - L)**2)

# r         = y[0]
# rd        = y[1]
# phi       = y[2]
# phid      = y[3]
# alpha     = y[4]
# l         = y[5]
# ld        = y[6]

E = 0
def prava(t, y):
    global E
    ca = cos(y[4])
    sa = sin(y[4])
    L = y[5]

    mphid = (y[0]*y[0]*y[3] - beta) / (L*L)

    men1p = (y[0] * y[0] + L * L + 2 * y[0] * L * ca) ** (1 / 2)
    men2p = (y[0] * y[0] + L * L - 2 * y[0] * L * ca) ** (1 / 2)

    men1 = men1p**3
    men2 = men2p**3

    clen = y[0] * y[0] * y[3] - beta

    return [
        # r' = \dot{r}
        y[1],
        # r'' = ...
        y[0] * y[3] * y[3] - K*M/2 * ((y[0] + L*ca) / men1 + (y[0] - L*ca) / men2),
        # phi' = \dot{phi}
        y[3],
        # phi'' = ...
        (1 / y[0]) * (-2*y[1]*y[3] + K*M*L*sa/2 * ( 1/men1 - 1/men2)),
        # alpha' = ...
        (1/(L*L)) * ( y[3] * (y[0] * y[0] + L*L) - beta),
        # L' = \dot{L}
        y[6],
        # L'' = ...
        clen*clen / (L**3) - K*M/2 * ((L + y[0]*ca)/men1 + (L - y[0]*ca)/men2) + 1 / (2*mu) * (-util_der(L)) - trenie * y[6],
    ]

def main():
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()


