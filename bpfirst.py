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
        self.vtop = 0.0
        self.vbottom = 0.0
        self.vi = 0.0
        self.EF = 0.0


def attr_vers_cond(sys,pars,values,ax,attr):
    Ts = [
        calculate.calculate_conductance(sys, Ef=pars.EF, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    ax = draw.simple_plot_data_2d([values, Ts], ax=ax, xlabel=attr, ylim=[-0.1, 6.1])
    return ax



def main():
    Nx = 10; Ny = 30;
    
    sys, lead = monolayer_black_phosphorus.make_system_two_terminal(Nx, Ny)
    #kwant.plot(sys)
    pars = SimplenameSpace()
    pars.vtop = 0.5
    pars.vbottom = -pars.vtop
    pars.vi = 0.01
    pars.EF = 0.3
    ##plot_bands
    data = calculate.caculate_energy_band(lead,pars)


    fig, axes = plt.subplots(1,4)
    
    

    ##bands
    draw.simple_plot_data_2d(data,axes[0],xlabel="k")

    ##磁场phi变化
    pars.EL, pars.ER, pars.W = [-0.01, 0.0, 0.0]
    phis = np.linspace(-3, 3, 10)
    attr_vers_cond(sys,pars,phis,axes[1],'phi')

    ##右边势能变化
    pars.EL, pars.W ,pars.phi = [0,0,0.02]
    ERs = np.linspace(-2.1, 2.1, 50)
    attr_vers_cond(sys, pars, ERs, axes[2], 'ER')

    ##无序强度变化
    pars.EL,pars.ER,pars.phi=[0,0,0]
    Ws = np.linspace(0, 1.0, 10)
    attr_vers_cond(sys, pars, Ws, axes[3], 'W')
    
    fig.set_figwidth(12)
    fig.set_figheight(3)
    
    #plt.show()
    
    return fig


if __name__ == "__main__":
    fig = main()
    # test_four()
