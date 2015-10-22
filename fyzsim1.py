#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
from numpy import *
from scipy.integrate import quad
from scipy.optimize import brentq as rootininterval

m = 1 # mass = 1kg
def calculate_period(a, b, fun, energy):
	"""Takes potential functions, energy level and x-bounds
	Returns 
	"""
	dt = lambda x: 1 / sqrt((energy - fun(x)) * 2 / m)
	return 2 * quad(dt, a, b)[0]

def total_energy(x, xd, pot):
	return pot(x) + 1/2 * m * xd*xd

def plot_problem_setting():
	"""Prints problem assignment as text
	`"""
	plt.subplot(2, 2, 4)

	plt.text(0.5, 0.5,
		"$m \\ddot{x} = - \\nabla U(x)$\n" + 
		"$U(x) = x^2(x-1)(x+1)^2$",
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=24)
	plt.title("Problem statement")
	

def plot_contour():
	def contour_section(x_range, energy_range, color='r'):
		x = linspace(x_range[0], x_range[1], 100)
		xd = linspace(-2, 2, 100)
		X, Y = meshgrid(x, xd)
		Z = [total_energy(xx, yy, ufun) for xx, yy in zip(X,Y)]
		levels = [l for l in energy_range]
		c = plt.contour(X,Y,Z, levels=levels, colors=color)
		plt.clabel(c, inline=1, fontsize=10)

	plt.subplot(2, 2, 3)
	plt.linestyles='solid'


	contour_section(
		(-1.5, -1),
		linspace(-1, -0.25, 4),
		'r')

	contour_section(
		(-1, 0),
		linspace(-0.1, -0.025, 4),
		'g')

	# green dot
	plt.plot([extrema[1]],[0], 'go')

	contour_section(
		(0, 1),
		linspace(-0.5, -0.1, 5),
		'b')

	# blue dot
	plt.plot([extrema[3]],[0], 'bo')

	contour_section(
		(-1.5, 1.5),
		linspace(0.2, 1, 5),
		'y')

	contour_section(
		(-1.5, 1.5),
		[0],
		'k')

	plt.xlim([-1.5, 1.2])
	plt.ylim([-2, 2])
	plt.grid(True)
	plt.xlabel('$x$')
	plt.ylabel('$\dot{x}$')
	plt.title('Phase portrait')

def plot_ux_more():
	def plot_ux_on_interval(a, b, color='r'):
		x = linspace(a, b, num=100)
		u = [ufun(z) for z in x]
		plt.plot(x, u, color)

	plt.subplot(2, 2, 1)

	plot_ux_on_interval(-1.5, -1, 'r')
	plot_ux_on_interval(-1, 0, 'g')
	plot_ux_on_interval(0, 1, 'b')
	plot_ux_on_interval(1, 1.2, 'y')

	plt.xlabel('position $x$')
	plt.ylabel('potential energy $U(x)$')
	plt.title('potential $U(x)$ depending on position $x$')

	plt.ylim([-1, 0.5])
	plt.xlim([-1.5, 1.2])

	plt.grid(True)
	#plt.show()

def plot_trajectories():
	def plot_section(energy_range, interval_range, plot_options='r'):
		eps = 10**-4
		points = 100

		e_axis = linspace(energy_range[0], energy_range[1], points)
		t_axis = []
		for energy in e_axis:
			ival = generate_interval(ufun, interval_range[0], interval_range[1], interval_range[2], energy)
			period = calculate_period(ival[0], ival[1], ufun, energy)
			t_axis.append(period)
		#plt.plot(e_axis, t_axis, plot_options)
		plt.plot(t_axis, e_axis, plot_options)

	def generate_interval(fun, fromx, middlex, tox, energy):
		leftbound = rootininterval(lambda x: fun(x) - energy, fromx, middlex) if fromx != -inf else -inf
		rightbound = rootininterval(lambda x: fun(x) - energy, middlex, tox) if tox != inf else inf
		return (leftbound, rightbound)

	plt.subplot(2, 2, 2)
	eps = 10**-4

	plot_section(
		(-0.8, -eps),
		(-inf, -100, extrema[0]),
		plot_options='r')

	plot_section(
		(extremaval[1], -eps),
		(extrema[0], extrema[1], extrema[2]),
		plot_options='g-')

	plot_section(
		(extremaval[3], -eps),
		(extrema[2], extrema[3], 1),
		plot_options='b')

	plot_section(
		(eps, 2.8),
		(-inf, 1, 100),
		plot_options='y')

	# plt.xlabel('Energy')
	# plt.ylabel('Period')
	plt.ylabel('energy')
	plt.xlabel('period length')

	# plt.xlim([-.6, .6])
	# plt.ylim([0, 15])
	plt.ylim([-1.0, .5])
	plt.xlim([0, 15])

	plt.title('Energy level on which the period length occurs')
	plt.grid(True)

# potential function U(x)
ufun = lambda x: x * x * (x - 1) * (x + 1) * (x + 1)
ucoef = [1, 1, -1, -1, 0, 0] # x^5 + x^4 - x^3 - x^2
udcoef = [5, 4, -3, -2, 0]
# ucoef[5] is the last (absolute) coefficient

extrema = sort([(c.real) for c in roots(udcoef) if not iscomplex(c)])
extremaval = [ufun(x) for x in extrema]
# extrema:
# u(-1) 	= 0
# u(-0.54) 	= -0.95
# u(0) 		= 0
# u(0.74) 	= -0.43
	

def main():
	plot_ux_more()
	plot_trajectories()
	plot_contour()
	plot_problem_setting()
	plt.show()

if __name__ == '__main__':
	main()