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

def main():
    Nx = 60; Ny = 60;

    sys, lead = monolayer_black_phosphorus.make_system_two_terminal(Nx, Ny)
    #kwant.plot(sys)
    pars = SimplenameSpace()
    pars.vtop = 0.5
    pars.vbottom = -pars.vtop
    pars.vi = 0.01
    pars.EF = 0.3

    bands = calculate.caculate_energy_band(lead, pars)

    fig, axes = plt.subplots(1, 2)
    ##bands
    draw.simple_plot_data_2d(bands, axes[0], xlabel="k")

    wave_density = calculate.wave_density(sys,0.3,pars)
    kwant.plotter.map(sys,wave_density,ax=axes[1])
    fig.set_figwidth(12)
    fig.set_figheight(6)
    fig.show()



if __name__ == "__main__":
    main()
    # test_four()
