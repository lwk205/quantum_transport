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
from matplotlib.mlab import griddata
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from common_functions import calculate, draw
from make_transport_system import monolayer_graphene
from scipy.sparse import coo_matrix 


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



def get_core_matrix(sysvx, sysvy, hop_sites, pars):
    sites = sysvy.sites
    sites_index = [[sites.index(site[0]), sites.index(site[1])] for site in hop_sites]
    site_pos = np.array([(site[0].pos + site[1].pos) / 2.0 for site in hop_sites])
    vx_coomatrix = calculate.get_system_hamiltonian(sysvx, pars, sparse=True)
    vy_coomatrix = calculate.get_system_hamiltonian(sysvy, pars, sparse=True)
    return site_pos[:,0],site_pos[:,1],sites_index, vx_coomatrix,vy_coomatrix


def get_current(ev,sites_index, vx_coomatrix, vy_coomatrix ):


    ######构造ev coo_matrix(一种稀疏矩阵的存储方式) 这个地方是计算流密度的关键之处
    row = []
    col = []
    ev_data = []
    cols = len(sites_index)
    rows = len(ev)
    for i, index in enumerate(sites_index):
        j1, j2 = index
        col.append(i)
        row.append(j1)
        ev_data.append(ev[j1])
        col.append(i)
        row.append(j2)
        ev_data.append(ev[j2])


    ev_coomatrix = coo_matrix((ev_data, (row, col)), shape=(rows, cols))

    vx = (ev_coomatrix.T.conjugate()).dot(vx_coomatrix.dot(ev_coomatrix)).diagonal()
    vy = (ev_coomatrix.T.conjugate()).dot(vy_coomatrix.dot(ev_coomatrix)).diagonal()

    ###datas
    vx = np.real(vx)
    vy = np.real(vy)



    return  vx, vy


def plot_current_field(xs,ys,vx,vy):
    ###define grid
    x_grid = np.linspace(xs.min(),xs.max(),200)
    y_grid = np.linspace(ys.min(),ys.max(),200)

    vxs = griddata(xs, ys, vx, x_grid, y_grid, interp='linear')
    vys = griddata(xs, ys, vy, x_grid, y_grid, interp='linear')


    speed = np.sqrt(vxs**2+vys**2)
    lw = 5*speed / speed.max()
    plt.figure()
    plt.streamplot(x_grid, y_grid, vxs,vys, color=lw, linewidth=2, cmap=plt.cm.autumn)
    plt.colorbar()
    plt.show()

def get_current_with_lead(sys,sysvx,sysvy,hop_sites,pars):
    wf = kwant.wave_function(sys, pars.EF, args=[pars])
    eva = wf(0)
    T = calculate.calculate_conductance(sys, pars.EF, 1, 0, pars)
    print T, eva.shape[0]

    xs, ys, sites_index, vx_coomatrix, vy_coomatrix = get_core_matrix(sysvx, sysvy, hop_sites, pars)

    vx0 = 0;
    ev0 = 0;
    vy0 = 0;
    for i in range(eva.shape[0]):
        ev = eva[i, :]
        vx, vy = get_current(ev, sites_index, vx_coomatrix, vy_coomatrix)

        vx0 = vx + vx0
        vy0 = vy + vy0
        ev0 = ev0 + np.abs(ev) ** 2
    return xs, ys, ev0, vx0, vy0

def random_many_times(N,sys,sysvx,sysvy,hop_sites,pars):
    ##############calculate randoms###########
    ev0,vx0,vy0 = 0, 0, 0
    for i in range(N):
        xs, ys, ev, vx, vy = get_current_with_lead(sys,sysvx,sysvy,hop_sites,pars)
        ev0 += ev
        vx0 += vx
        vy0 += vy

    return xs, ys, ev0/N, vx0/N, vy0/N

def main():

    Nx = 60
    Ny = 40
    hop_sites, sys, sysvx, sysvy, lead_left, lead_right = \
           monolayer_graphene.make_monolayer_graphene_system_with_v(Nx,Ny)
    # kwant.plot(sys)
    ####elecron wave####
    pars = SimplenameSpace()
    pars.EL = -0.1
    pars.ER = 0.1
    pars.phi = 0.007 * 2.3
    pars.EF = 0.0



    ########non_random
    # pars.W = 0.0
    # xs, ys, ev, vx, vy = get_current_with_lead(sys, sysvx, sysvy, hop_sites, pars)


    ######randoms###
    pars.W = 1.0
    N = 10
    random.seed(1000)
    xs, ys, ev, vx, vy = random_many_times(N, sys, sysvx, sysvy, hop_sites, pars)

    plot_current_field(xs, ys, vx, vy)
    kwant.plotter.map(sys,ev)


if __name__ == "__main__" :
    main()

    
