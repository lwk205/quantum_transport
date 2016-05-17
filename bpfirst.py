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
        self.vtop = -0.5
        self.vbottom = 0.5
        self.vi = 0.0
        self.EF = 0.0


def attr_vers_cond(sys,pars,values,attr):
    Ts = [
        calculate.calculate_conductance(sys, Ef=pars.EF, out_leads=1, in_leads=0, pars=pars)
        for vars(pars)[attr]  in values
        ]
    #draw.simple_plot_data_2d_without_ax([values, Ts],  xlabel=attr,ylabel="Conductance")
    return np.array([values, Ts])



def main():
    Nx = 10; Ny = 20;
    
    sys, lead = monolayer_black_phosphorus.make_system_two_terminal(Nx, Ny)

#    pars.vtop = 0.5
#    pars.vbottom = -pars.vtop
#    pars.vi = 0.01
#    pars.EF = 0.35
    El =-0.4
    vtop = -0.5
    vbot = 0.5#-vtop
    #####plot_bands
#    pars = SimplenameSpace()
#    pars.vbottom = vbot
#    pars.vtop = vtop
#    pars.EL = El
#    bands = calculate.caculate_energy_band(lead,pars)
#    print np.min(np.abs(bands[1]))
#    draw.simple_plot_data_2d_without_ax(bands,xlabel="k",ylabel="energy")


    ####磁场EF变化
#    pars = SimplenameSpace()
#    pars.vbottom = vbot
#    pars.vtop = vtop    
#    pars.EL, pars.ER, pars.W = [El, 0.0, 0.0]
#    EFs = np.linspace(-0.5, 0.0 , 50)
#    attr_vers_cond(sys,pars,EFs,'EF')

    ##磁场phi变化
#    pars = SimplenameSpace()
#    pars.EL, pars.ER, pars.W = [-0.4, 0.0, 0.0]
#    pars.EF = 0.0
#    phis = np.linspace(-3, 3, 50)
#    attr_vers_cond(sys,pars,phis,'phi')

    ##右边势能ER变化    
#    pars = SimplenameSpace()
#    pars.vbottom = vbot
#    pars.vtop = vtop
#    pars.EF = -0.0
#    pars.ER, pars.W ,pars.phi = [-0.4, 0.0, 0.04]
#    ERs = np.linspace(-0.55, -0.2, 50)
#    attr_vers_cond(sys, pars, ERs, 'EL')

#    ##无序强度W变化

    def random_cal():
        pars = SimplenameSpace()
        pars.vbottom = vbot
        pars.vtop = vtop
        pars.EF = -0.0
        pars.EL,pars.ER,pars.phi=[-0.20, -0.4, 0.04]
        Ws = np.linspace(0, 10.0, 100)
        data = attr_vers_cond(sys, pars, Ws, 'W')
        return data
        
    N = 200
    All = 0.0
    random.seed(10000)
    for i in range(N):
        All = All +random_cal()
        print i
    
    
    
        
#
#
#    pars = SimplenameSpace()
#    pars.vbottom = vbot
#    pars.vtop = vtop
#    pars.EF = -0.1
#    wave_density = calculate.wave_density(sys,pars.EF,pars,lead_nr=0)
#    kwant.plotter.map(sys,wave_density)
    #plt.show()
    return All/N



if __name__ == "__main__":
    All = main()
    plt.figure()
    plt.plot(All[0,:],All[1,:])
    plt.show()

