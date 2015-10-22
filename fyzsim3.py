#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.integrate import ode
import scipy

from scipy import zeros, sin, cos

L = 1
K = 1
M = 10
mu = 1

# 
r, rd = 3, -1
phi, phid = 0, .25
psi, psid = scipy.pi / 4, -.1

alpha, alphad = phi - psi, phid - psid

y0, t0 = [r, rd, phi, phid, alpha], 0

beta = 2*r*r * phid*phid - 2*L*L*(alphad - phid)


t1 = 4
dt = .01


# r 	 = y[0]
# rd 	 = y[1]
# phi 	 = y[2]
# phid	 = y[3]
# alpha	 = y[4]

def prava(t, y):
	ca = cos(y[4])
	sa = sin(y[4])

	men1 = (y[0] * y[0] + L * L + 2 * y[0] * L * ca) ** (3 / 2)
	men2 = (y[0] * y[0] + L * L - 2 * y[0] * L * ca) ** (3 / 2)

	return [
		# y' = \dot{y}
		y[1],
		# y'' = ...
		y[0] * y[3] * y[3] - K*M/4 * ((2*y[0] + 2*L*ca) / men1 + (2*y[0] - 2*L*ca) / men2),
		# phi' = \dot{phi}
		y[3],
		# phi'' = ...
		(-1 / (y[0] * y[0])) * (2*y[0]*y[1]*y[3] + K*M/4 * ( -2*y[0]*L*sa/men1 + 2*y[0]*L*sa/men2 )),
		# alpha' = ...
		(1/(L*L)) * ( y[3] * (y[0] * y[0] + L*L) - beta / 2),
	]


odes = ode(prava).set_integrator('dopri5')
odes.set_initial_value(y0, t0)


pocet = t1 / dt

pp1 = pocet + 1

T = zeros(pp1)
R = zeros(pp1)
PHI = zeros(pp1)
PSI = zeros(pp1)

X1 = zeros(pp1)
Y1 = zeros(pp1)

X2 = zeros(pp1)
Y2 = zeros(pp1)

real_pocet = 0
while odes.successful() and odes.t < t1:
	T[real_pocet] = odes.t
	R[real_pocet] = odes.y[0]
	PHI[real_pocet] = odes.y[2]
	PSI[real_pocet] = PHI[real_pocet] - odes.y[4]

	odes.integrate(odes.t + dt)
	
	real_pocet += 1

T = T[:real_pocet]
R = R[:real_pocet]
PHI = PHI[:real_pocet]
PSI = PSI[:real_pocet]

X1 = R * cos(PHI) + L * cos(PSI)
Y1 = R * sin(PHI) + L * sin(PSI)
X2 = R * cos(PHI) - L * cos(PSI)
Y2 = R * sin(PHI) - L * sin(PSI)

print(real_pocet)
from matplotlib import pyplot
pyplot.plot(X1, Y1)
pyplot.plot(X2, Y2)
pyplot.plot([0], [0], 'ro')
pyplot.grid(True)
pyplot.show()