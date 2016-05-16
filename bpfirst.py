# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-15
Contact: gengyusheng@gmail.com
内容： 单层黑磷输运性质计算
"""

from __future__ import division
import random
import kwant
import matplotlib.pyplot as plt
import numpy as np
from common_functions import calculate, draw
from make_transport_system import monolayer_black_phosphorus


class SimplenameSpace:
    def __init__(self):
        self.EL = 0.0
        self.ER = 0.0
        self.phi = 0.0
        self.W = 0.0
        self.EU = 0.0
        self.ED = 0.0
        self.vtop = 0.5
        self.vbottom = -0.5
        self.vi = 0.01
        self.EF = 0.25


def attr_vers_cond(sys,pars,values,ax,attr):
    Ts = [
        calculate.calculate_conductance(sys, Ef=pars.EF, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    ax = draw.simple_plot_data_2d([values, Ts], ax=ax, xlabel=attr)
    return ax



def main():
    Nx = 20; Ny = 40;
    
    sys, lead = monolayer_black_phosphorus.make_system_two_terminal(Nx, Ny)

#    pars.vtop = 0.5
#    pars.vbottom = -pars.vtop
#    pars.vi = 0.01
#    pars.EF = 0.35

    ##plot_bands
    # bands = calculate.caculate_energy_band(lead,pars)
    # draw.simple_plot_data_2d(bands,axes[0],xlabel="k")

    fig, axes = plt.subplots(1,3)
    ##磁场EF变化
    pars = SimplenameSpace()
    pars.EL, pars.ER, pars.W = [-0.3, -0.6, 0.0]
    EFs = np.linspace(0.0, 0.6 , 100)
    attr_vers_cond(sys,pars,EFs,axes[0],'EF')

#    ##磁场phi变化
#    pars = SimplenameSpace()
#    pars.EL, pars.ER, pars.W = [0.0, 0.0, 0.0]
#    phis = np.linspace(-3, 3, 50)
#    attr_vers_cond(sys,pars,phis,axes[1],'phi')

    ##右边势能ER变化    
    pars = SimplenameSpace()
    pars.EF = 0.0
    pars.EL, pars.W ,pars.phi = [-0.3, 0.0, 0.00]
    ERs = np.linspace(-3.0, 0.3, 100)
    attr_vers_cond(sys, pars, ERs, axes[1], 'ER')

    ##无序强度W变化
    pars = SimplenameSpace()
    pars.EF = 0.0
    pars.EL,pars.ER,pars.phi=[-0.3, -0.6, 0.00]
    Ws = np.linspace(0, 5.0, 100)
    attr_vers_cond(sys, pars, Ws, axes[2], 'W')
    
    fig.set_figwidth(12)
    fig.set_figheight(3)
#
#    pars = SimplenameSpace()
#    wave_density = calculate.wave_density(sys,pars.EF,pars,lead_nr=0)
#    kwant.plotter.map(sys,wave_density)
#    #plt.show()
    



if __name__ == "__main__":
    main()

