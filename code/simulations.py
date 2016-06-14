# -*- coding: utf-8 -*-
# plot_figures: Plots figures for the reference implementation of
# [1]
# [1]: "Multiple dynamical modes of thalamic relay neurons: rhythmic
# bursting and intermittent phase-locking", Wang, X-J, Neuroscience,
# 59(1), pg.21-31, 1994.
#
# Copyright (C) 2016 Georgios Is. Detorakis (gdetor@protonmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from neuron_model import *


########################################################################
# In the following code, if you want to use BDF or Adams methods for
# temporal integrations, please use: 
#            simulation.run(sim_time,
#                           fname_for_mem_pot,
#                           fname_for_time,
#                           'vode', 'Adams' [or 'BDF'] ) 
########################################################################
def simulate_figure1(diagram=False):
    print("Running simulation #1")

    if diagram is True:
        samples = 150
        freq = np.linspace(0.1, 15, samples)
        Iext = -1.0                                   

        name = "Fig1A"
        idx = 0
        for i, f in enumerate(freq):
            # period = np.round(1./f, 3) * 1000
            period = 1./f * 1000
            dur = np.linspace(0, period, samples)
            for j, d in enumerate(dur):
                print('Frequency: {}, Duration: {}'.format(f, d))
                nm = 'Fig1A_V'+str(idx)
                print('Fname: {}'.format(nm))
                simulation = WangModel(base+"params_figure1.cfg",
                                       frequency=f,
                                       duration=d,
                                       inp_type='pulse',
                                       amplitude=Iext,
                                       chunks=100)
                simulation.run(5*period, name+"_V"+str(idx),
                               name+"_T"+str(idx),
                               'vode',
                               'Adams')
                idx += 1

    freq = [10, 5, 0.5]           # Pulse frequency in Hz
    dur = [10, 120, 1200]         # Pulse duration in ms
    Iext = -1.0                   # Pulse current in μA/cm^2

    name = "Fig1B"
    for i in range(len(freq)):
        simulation = WangModel(base+"params_figure1.cfg",
                               frequency=freq[i],
                               duration=dur[i],
                               inp_type='pulse',
                               amplitude=Iext,
                               chunks=100)
        simulation.run(2000, name+"_V"+str(i), name+"_T"+str(i), 'vode', 'Adams')


def simulate_figure2():
    print("Running simulation #2")
    Iext = [3.0, 0.0, -0.47, -0.6, -0.8, -1.3, -1.4, -2.0]
    name = "Fig2"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure2.cfg",
                               inp_type='const',
                               amplitude=i)
        simulation.run(2000, name+"_V"+str(j), name+"_T"+str(j), 'vode', 'Adams')


def simulate_figure3():
    print("Running simulation #3")
    name = "Fig3"
    simulation = WangModel(base+"params_figure3.cfg",
                           inp_type="const",
                           amplitude=-0.91)
    simulation.run(8000, name+"_V", name+"_T", 'vode', 'Adams')


def simulate_figure4():
    print("Running simulation #4")
    samples = 100

    Iext = np.linspace(0, -2.0, samples)
    name = "Fig4A"
    for i, I in enumerate(Iext):
        nm = 'Fig4A_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure4.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(2000, name+"_V"+str(i), 
                       name+"_T"+str(i),
                       'vode',
                       'Adams')

    Iext = np.round(np.genfromtxt('../data/data.dat'), 2)
    name = "Fig4AReal"
    for i, I in enumerate(Iext):
        nm = 'Fig4AReal_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure4.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(2000, name+"_V"+str(i), 
                       name+"_T"+str(i),
                       'vode',
                       'Adams')

    samples = 100
    Iext = np.linspace(0, -2.0, samples)
    name = "Fig4B"
    for i, I in enumerate(Iext):
        nm = 'Fig4B_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure4b.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(2000, name+"_V"+str(i), 
                       name+"_T"+str(i),
                       'vode',
                       'Adams')

    Iext = np.round(np.genfromtxt('../data/data1.dat'), 2)

    name = "Fig4BReal"
    for i, I in enumerate(Iext):
        nm = 'Fig4BReal_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure4b.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(2000, name+"_V"+str(i), 
                       name+"_T"+str(i),
                       'vode',
                       'Adams')


def simulate_figure5():
    print("Running simulation #5")
    Iext = [-1.4, -1.51, -1.6]

    name = "Fig5A"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure5a.cfg",
                               frequency=10.0,
                               duration=80.0,
                               inp_type="pulse",
                               amplitude=i,
                               chunks=1)
        simulation.run(8000,
                       name+"_V"+str(j),
                       name+"_T"+str(j),
                       'vode',
                       'BDF')

    Iext = [-1.2, -1.51, -1.8]
    name = "Fig5B"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure5b.cfg",
                               frequency=10.0,
                               duration=80.0,
                               inp_type="pulse",
                               amplitude=i,
                               chunks=1)
        simulation.run(8000,
                       name+"_V"+str(j),
                       name+"_T"+str(j),
                       'vode',
                       'BDF')


if __name__ == '__main__':
    base = "params/"

    # Build data for Figure 1
    simulate_figure1(diagram=True)

    # Build data for Figure 2
    simulate_figure2()

    # # Build data for Figure 3
    simulate_figure3()

    # Build data for Figure 4
    simulate_figure4()

    # Build data for figure 5
    simulate_figure5()
