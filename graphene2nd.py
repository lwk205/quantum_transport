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


def eig_solves(sys):
    pars = SimplenameSpace()
    pars.EL = -0.
    pars.ER = -0.
    pars.phi = 0.000
    hams = calculate.get_system_hamiltonian(sys,pars,sparse=True)
    eigs, evs = eigsh(hams, which="SM",k=20)
    # print np.sort(eigs)
    return eigs, evs


#def get_current(sys_vector, sysvx, sysvy, pars):
   

#def main():
Nx = 60
Ny = 40

hop_sites, sys, sysvx, sysvy, lead_left, lead_right = monolayer_graphene.make_monolayer_graphene_system_with_v(Nx, Ny)

sites = sys.sites
hop_sites = list(hop_sites)
sites_index = [[sites.index(site[0]),sites.index(site[1])] for site in hop_sites]
site_pos = np.array([(site[0].pos + site[1].pos)/2.0 for site in hop_sites])

pars = SimplenameSpace()
pars.EL = -0.0
pars.ER = -0.0
pars.phi = 0.05#0.007*2.3
hams = calculate.get_system_hamiltonian(sys,pars,sparse=True)
eigs, evs = eigsh(hams, which="SM",k=20)
vx_coomatrix = calculate.get_system_hamiltonian(sysvx,pars,sparse=True)
vy_coomatrix = calculate.get_system_hamiltonian(sysvy, pars, sparse=True)

####get current####
ev = evs[:,0]
##create the matrix to compute conjT(ev).v.ev##
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

ev_coomatrix = coo_matrix((ev_data, (row,col) ), shape=(rows,cols))

vx = (ev_coomatrix.T.conjugate()).dot(vx_coomatrix.dot(ev_coomatrix)).diagonal()
vy = (ev_coomatrix.T.conjugate()).dot(vy_coomatrix.dot(ev_coomatrix)).diagonal()

vx = np.real(vx)
vy = np.real(vy)

X  = np.ones((cols,1))
X = np.kron(X,site_pos[:,0])
Y  = np.ones((1,cols))
Y = np.kron(Y,site_pos[:,1].reshape(cols,1))
Ux = np.diag(vx)
Uy = np.diag(vy)
speed = np.sqrt(Ux*Ux + Uy*Uy)
lw = 5*speed 

plt.figure()
plt.streamplot(X, Y, Ux, Uy, density=0.6, color='k', linewidth=lw)
plt.show()

#data_save = np.array([site_pos[:,0],site_pos[:,1],vx,vy]).T
#np.savetxt("current_dens.txt",data_save)


    # print evs.shape
#    return evs, eigs, sys, ham



#if __name__ == "__main__" :
#    evs, eigs, sys, ham = main()
#    wf = np.abs(evs)**2
#    kwant.plotter.map(sys, wf[:,0])

    
