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
        calculate.calculate_conductance(sys, pars.EF, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    draw.simple_plot_data_2d_without_ax([values, Ts],  xlabel=attr,ylabel="Conductance")
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

def main():
    Nx = 32
    Ny = 32

    sys, lead_left, lead_right = monolayer_graphene.make_monolayer_graphene_system(Nx,Ny)
    pars = SimplenameSpace()
    pars.EL = 0.0
    pars.phi = 0.0
    EFs = np.linspace(0.000001,0.2)
    attr_vers_cond(sys,pars,EFs,"EF")

    #
    # pars = SimplenameSpace()
    # pars.EF = 0.0
    # pars.EL = 0.0
    # pars.ER = 0.0
    # phis = np.linspace(0.0,7,2)
    # landau_levels(sys,phis,pars,sparse=False)

    # pars = SimplenameSpace()
    # plot_bands(lead_left,pars)
    # pars = SimplenameSpace()
    # pars.phi = 0.0
    # pars.EL = 1.0
    # pars.ER = 0.0
    # ham = calculate.get_system_hamiltonian(sys.finalized(),pars, sparse=True)
    # landau_dos(sys,pars)
    # egs = calculate.eig_system_hamiltonian(sys.finalized(),pars)
    # return ham

if __name__ == "__main__":
    main()
