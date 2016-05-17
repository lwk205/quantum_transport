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
        calculate.calculate_conductance(sys, Ef=pars.EF, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    draw.simple_plot_data_2d_without_ax([values, Ts],  xlabel=attr,ylabel="Conductance")
    # return np.array([values, Ts])

def landau_levels(sys,values,pars,ylim):
     egs = [calculate.eig_system_hamiltonian(sys.finalized(), pars) for pars.phi in values]
     draw.simple_plot_data_2d_without_ax([values, egs], ylim = ylim,xlabel="$\phi$" , ylabel="Energy")

def landau_dos(sys,values,pars):
    egs = calculate.eig_system_hamiltonian(sys.finalized(),pars)
    delta = 0.0001
    N = 0
    while N<len(egs):
        n = 0
        while egs[N+1] < egs[N] + delta:
            pass


def plot_bands(lead,pars):
    data = calculate.caculate_energy_band(lead,pars)
    draw.simple_plot_data_2d_without_ax(data, xlabel="k", ylabel="Energy")

def main():
    Nx = 10
    Ny = 10

    sys, lead_left, lead_right = monolayer_graphene.make_monolayer_graphene_system(Nx,Ny)


    pars = SimplenameSpace()
    phis = np.linspace(0.0,7.3,100)
    landau_levels(sys,phis,pars,ylim=None)

    pars = SimplenameSpace()
    plot_bands(lead_left,pars)

if __name__ == "__main__":
    main()
