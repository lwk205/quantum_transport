# -*- coding: utf-8 -*-
"""
Created on Fri May 27 06:16:05 2016

@author: rainson
"""

#from make_transport_system import monolayer_graphene
from make_transport_system import lattices
import kwant
from cmath import exp
from scipy.sparse.linalg import eigsh
import numpy as np
from scipy.linalg import eigh
import random

Nx = 40
Ny = 30
t = -1.0
phi = 0.00*2.3
EL = 0.0
ER = 0.0
W = 0.0
sys = kwant.Builder()
sysvx = kwant.Builder()
sysvy = kwant.Builder()

graphene, subs = lattices._make_graphene_lattice()
a, b, c, d = subs 




def make_system(Nx,Ny):
    
    ####set onsite#########
    sites = [sub(x,y) 
             for sub in subs 
             for x in range(Nx)
             for y in range(Ny)]
                 
                 
    for site in sites:
        pos_x, pos_y = site.pos
        pot = EL+(ER-EL)*pos_x/Nx + W*(2.0*random.random()-1.0)
        sys[site] =  pot
        sysvx[site] = 0.0
        sysvy[site] = 0.0
        
    
    ##############set hoppings##########        
    hoppings = [(b(x,y),a(x,y)) for x in range(Nx) for y in range(Ny)]
    hoppings.extend([(b(x+1,y),a(x,y)) for x in range(Nx-1) for y in range(Ny)])
    hoppings.extend([(a(x,y+1),d(x,y)) for x in range(Nx) for y in range(Ny-1)])
    hoppings.extend([(c(x,y),b(x,y)) for x in range(Nx) for y in range(Ny)])
    hoppings.extend([(d(x,y),c(x,y)) for x in range(Nx) for y in range(Ny)])
    hoppings.extend([(d(x,y),c(x+1,y)) for x in range(Nx-1) for y in range(Ny)])
    
    for hopping in hoppings:
        site_to, site_from = hopping
        x_to, y_to = site_to.pos
        x_from, y_from = site_from.pos
        hop_value = t * exp(-0.5j*phi*(x_to-x_from)*(y_to+y_from))
        sys[hopping] = hop_value
        sysvx[hopping] = -1.0j*hop_value*(x_to-x_from)
        sysvy[hopping] = -1.0j*hop_value*(x_to-x_from)
    
    return sys, sysvx, sysvy
    
    
def make_leads(Ny):
    a, b, c, d = subs 
    hoppings = [
            ((0, 0), b, a),
            ((1, 0), b, a),
            ((0, -1), d, a),
            ((0, 0), c, b),
            ((0, 0), d, c),
            ((-1, 0), d, c)
            ]

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[sub(0,ny) for sub in subs  for ny in range(Ny)] ]= EL
    lead_left[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = t
    lead_right[[sub(Nx-1,ny) for sub in subs  for ny in range(Ny)] ]= ER
    lead_right[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = t
    return lead_left, lead_right
    
    
 
sys, sysvx, sysvy = make_system(Nx,Ny)
lead_left, lead_ringht = make_leads(Ny)
sys = sys.finalized()
sysvx = sysvx.finalized()
sysvy = sysvy.finalized()
ham = sys.hamiltonian_submatrix()
eig_values, eig_vectors = eigsh(ham, k=15, which = "SM")
#N = len(eig_values)
kwant.plotter.map(sys,np.abs(eig_vectors[:,4])**2)

#kwant.plot(sys)

