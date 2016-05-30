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


def plot_bands(lead,pars):
    data = calculate.caculate_energy_band(lead,pars)
    draw.simple_plot_data_2d_without_ax(data, xlabel="k", ylabel="Energy")


def eig_solves(sys):
    pars = SimplenameSpace()
    pars.EL = -0.1
    pars.ER = -0.1
    pars.phi = 0.007*2.3
    hams = calculate.get_system_hamiltonian(sys,pars,sparse=True)
    eigs, evs = eigsh(hams, which="SM",k=20)
    # print np.sort(eigs)
    return eigs, evs


def get_current(sys_vector, sysvx, sysvy, pars):
    vx_coomatrix = calculate.get_system_hamiltonian(sysvx,pars,sparse=True)
    vy_coomatrix = calculate.get_system_hamiltonian(sysvy, pars, sparse=True)

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

    
