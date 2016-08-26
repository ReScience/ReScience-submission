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
def simulate_figure1():
    print("Running test for comparing integration methods")
    base = './params/'
    Iext = -0.55
    name = ['test_dopri', 'test_adams', 'test_bdf']

    simulation = WangModel(base+"params_figure3.cfg",
                           inp_type='const',
                           amplitude=Iext)
    simulation.run(3000, name[0]+"_V", 'dopri5')
    simulation.run(3000, name[1]+"_V", 'vode', 'Adams')
    simulation.run(3000, name[2]+"_V", 'vode', 'BDF')

    base = '../data/'
    def spk(x):
        tmp = np.zeros((x.shape[0], ))
        spikes = peak_local_max(x, min_distance=10)
        spikes = spikes[x[spikes] > 0]
        tmp[spikes] = 1
        return tmp

    dopri = np.load(base+'test_dopri_V.npy')
    adams = np.load(base+'test_adams_V.npy')
    bdf = np.load(base+'test_bdf_V.npy')
    spks1 = spk(dopri[:, 1]).astype('i')
    spks2 = spk(adams[:, 1]).astype('i')
    spks3 = spk(bdf[:, 1]).astype('i')
    s1 = np.where(spks1 > 0)[0]
    s2 = np.where(spks2 > 0)[0]
    s3 = np.where(spks3 > 0)[0]
    isi1 = np.diff(dopri[spks1==1, 0])
    isi2 = np.diff(adams[spks2==1, 0])
    isi3 = np.diff(bdf[spks3==1, 0])
    CV1 = np.mean(isi1) / np.std(isi1)
    CV2 = np.mean(isi2) / np.std(isi2)
    CV3 = np.mean(isi3) / np.std(isi3)
    print('sum(Spikes dopri5 - spikes Adams):  ', np.abs((s1 - s2).sum()))
    print('sum(Spikes dopri5 - spikes BDF):  ', np.abs((s1 - s3).sum()))
    print('CV dopri5:', CV1)
    print('CV Adams:', CV2)
    print('CV BDF:', CV3)


def simulate_figure2(diagram=False):
    print("Running simulation #2")

    if diagram is True:
        samples = 100
        freq = np.linspace(0.1, 15, samples)
        Iext = -1.0                                   

        name = "Fig2A"
        idx = 0
        for i, f in enumerate(freq):
            period = 1./f * 1000
            dur = np.linspace(0, period, samples)
            for j, d in enumerate(dur):
                print('Frequency: {}, Duration: {}'.format(f, d))
                nm = 'Fig2A_V'+str(idx)
                print('Fname: {}'.format(nm))
                simulation = WangModel(base+"params_figure2.cfg",
                                       frequency=f,
                                       duration=d,
                                       inp_type='periodic',
                                       amplitude=Iext,
                                       chunks=100,
                                       store_stim=False)
                simulation.run(15*period, name+"_V"+str(idx),
                               'vode',
                               'Adams')
                idx += 1

    freq = [10, 5, 0.5]           # Pulse frequency in Hz
    dur = [10, 120, 1200]         # Pulse duration in ms
    Iext = -1.0                   # Pulse current in μA/cm^2

    name = "Fig2B"
    for i in range(len(freq)):
        simulation = WangModel(base+"params_figure2.cfg",
                               frequency=freq[i],
                               duration=dur[i],
                               inp_type='periodic',
                               amplitude=Iext,
                               chunks=100)
        simulation.run(5000, name+"_V"+str(i), 'vode', 'Adams')


def simulate_figure3():
    print("Running simulation #3")
    Iext = [3.0, 0.0, -0.45, -0.455, -0.47, -0.55, -0.6, -0.8, -1.3, -1.4, -2.0]
    name = "Fig3"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure3.cfg",
                               inp_type='const',
                               amplitude=i)
        simulation.run(6000, name+"_V"+str(j), 'vode', 'Adams')


def simulate_figure4():
    print("Running simulation #4")
    Iext = np.linspace(-0.433, -0.55, 20)
    
    name = "Fig4A"
    for i, j in enumerate(Iext):
        simulation = WangModel(base+"params_figure4a.cfg",
                               inp_type="const",
                               amplitude=j,
                               chunks=1)
        simulation.run(1500, name+"_V"+str(i), 'vode', 'Adams')

    name = "Fig4B1"
    simulation = WangModel(base+"params_figure4b.cfg",
                           duration=100,
                           inp_type="pulse",
                           amplitude=-0.47,
                           pulse_ampl=-1.2,
                           chunks=1,
                           store_stim=True)
    simulation.run(1000, name+"_V", 'vode', 'Adams')

    name = "Fig4B2"
    simulation = WangModel(base+"params_figure4b.cfg",
                           duration=100,
                           inp_type="pulse",
                           amplitude=-0.47,
                           pulse_ampl=0.25,
                           chunks=1,
                           store_stim=True)
    simulation.run(1000, name+"_V", 'vode', 'Adams')


def simulate_figure5():
    print("Running simulation #5")
    name = "Fig5"
    simulation = WangModel(base+"params_figure5.cfg",
                           inp_type="const",
                           amplitude=-0.95)
    simulation.run(5000, name+"_V", 'vode', 'Adams')


def simulate_figure6():
    print("Running simulation #6")

    Iext = np.genfromtxt('../data/data1or.dat')[:, 0]
    name = "Fig6AReal"
    for i, I in enumerate(Iext):
        nm = 'Fig6AReal_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure6a.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(5000, name+"_V"+str(i), 'vode', 'Adams')

    Iext = np.genfromtxt('../data/data2or.dat')[:, 0]
    name = "Fig6BReal"
    for i, I in enumerate(Iext):
        nm = 'Fig6BReal_V'+str(i)
        print('Fname: {}, Iext: {}'.format(nm, I))
        simulation = WangModel(base+"params_figure6b.cfg",
                               inp_type='const',
                               amplitude=I)
        simulation.run(5000, name+"_V"+str(i), 'vode', 'Adams')


def simulate_figure7():
    print("Running simulation #7")
    Iext = [-1.4, -1.5, -1.6]

    name = "Fig7A"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure7a.cfg",
                               frequency=10.0,
                               duration=80.0,
                               inp_type="periodic",
                               amplitude=i,
                               chunks=1,
                               store_stim=True)
        simulation.run(40000, name+"_V"+str(j), 'vode', 'Adams')

    Iext = [-1.2, -1.5, -1.8]
    name = "Fig7B"
    for j, i in enumerate(Iext):
        simulation = WangModel(base+"params_figure7b.cfg",
                               frequency=10.0,
                               duration=80.0,
                               inp_type="periodic",
                               amplitude=i,
                               chunks=1,
                               store_stim=True)
        simulation.run(20000, name+"_V"+str(j), 'vode', 'Adams')


if __name__ == '__main__':
    base = "./params/"

    # Build data for Figure 1
    simulate_figure1()

    # Build data for Figure 2
    simulate_figure2(diagram=True)

    # Build data for Figure 3
    simulate_figure3()

    # Build data for Figure 4
    simulate_figure4()

    # Build data for Figure 5
    simulate_figure5()

    # Build data for figure 6
    simulate_figure6()

    # Build data for figure 7
    simulate_figure7()
