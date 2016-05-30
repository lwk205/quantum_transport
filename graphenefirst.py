# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-17
Contact: gengyusheng@gmail.com
内容： 单层graphene体系物理性质计算及图示
"""

from __future__ import division
import random
import kwant
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from common_functions import calculate, draw
from make_transport_system import monolayer_graphene

class SimplenameSpace:
    def __init__(self):
        self.EL = 0.0
        self.ER = 0.0
        self.phi = 0.0
        self.W = 0.0
        self.EF = 0.0


def attr_vers_cond(sys,pars,values,attr):
    Ts = [
        calculate.calculate_conductance(sys, 0.0, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    draw.simple_plot_data_2d_without_ax([values, Ts],  xlabel=attr,ylabel="Conductance",ylim=[0.0,1.2])
    # return np.array([values, Ts])

def landau_levels(sys,values,pars, sparse):
     egs = [np.sort(calculate.system_eigvalues(sys, pars, sparse)) for pars.phi in values]
     draw.simple_plot_data_2d_without_ax([values, egs],xlabel="$\phi$" , ylabel="Energy")

def landau_dos(sys,pars,sparse):
    egs = calculate.system_eigvalues(sys, pars, sparse )
    delta = 0.00001
    N = 0
    egg = np.arange(-4.1,4.1,delta)
    ns = []
    for eg in egg:
        n = 0

        if N < len(egs) :
            while eg <= egs[N] < eg + delta:
                N = N+1
                n = n + 1
                if N >= len(egs):
                    break
            ns.append(n)
        else:
            ns.append(0)


    draw.simple_plot_data_2d_without_ax([egg, ns],ylim=[0,3.5],xlim=[-3.0,3.0])
    return None


def plot_bands(lead,pars):
    data = calculate.caculate_energy_band(lead,pars)
    draw.simple_plot_data_2d_without_ax(data, xlabel="k", ylabel="Energy")

def random_map(sys):
    pars = SimplenameSpace()
    pars.EL = -0.1
    pars.ER = 0.1
    pars.phi = 0.007*2.3
    pars.EF = 0.0
    pars.W = 0.6
    wd = calculate.wave_density(sys,pars,lead_nr=0)
    T = calculate.calculate_conductance(sys,pars.EF,1,0,pars)
    return  T, wd

def eig_solves(sys):
    pars = SimplenameSpace()
    pars.EL = -0.1
    pars.ER = -0.1
    pars.phi = 0.007*2.3
    hams = calculate.get_system_hamiltonian(sys,pars,sparse=True)
    eigs, evs = eigsh(hams, which="SM",k=20)
    # print np.sort(eigs)
    return eigs, evs

def many_times_random_map(sys,N=200):
    wd0, T0 = 0, 0
    N = 200
    random.seed(10000)
    for i in range(N):
        T, wd = random_map(sys)
        wd0 = wd0 + wd
        # print i, T
        T0 = T0 + T
    kwant.plotter.map(sys, wd0 / N)
    # print T0 / N

def get_current(sys_vector,sysvx,sysvy):
    vx_coomatrix = calculat

def main():
    Nx = 60
    Ny = 40

    sys, sysvx, sysvy, lead_left, lead_right = monolayer_graphene.make_monolayer_graphene_system_with_v(Nx, Ny)
    pars = SimplenameSpace()
    pars.EL = -0.1
    pars.ER = -0.1
    pars.phi = 0.007*2.3
    hams = calculate.get_system_hamiltonian(sys,pars,sparse=True)
    ham = calculate.get_system_hamiltonian(sys,pars,sparse = False)
    eigs, evs = eigsh(hams, which="SM",k=20)
    # print evs.shape
    return evs, eigs, sys, ham



if __name__ == "__main__" :
    evs, eigs, sys, ham = main()
    wf = np.abs(evs)**2
    kwant.plotter.map(sys, wf[:,0])

    
