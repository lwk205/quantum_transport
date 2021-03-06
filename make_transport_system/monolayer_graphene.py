# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-17
Contact: gengyusheng@gmail.com
内容： 单层graphene的体系构建
"""

from __future__ import division
import kwant
import random
from cmath import exp
import math
from lattices import  _make_graphene_lattice
import matplotlib.pylab



def make_monolayer_graphene_system(Nx,Ny):
    """
    :param 长度Nx:
    :param 宽度Ny:
    :return sys without lead , lead_left, lead_right:
    """
    t = 1.0


    def potential(site,pars):
        x, y = site.pos
        EL, ER = pars.EL, pars.ER
        return EL+(ER-EL)*x/(Nx-1) + pars.W*(2.0*random.random()-1.0)


    def hopping(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))

    def right_pot(site,pars):
        return pars.ER

    def left_pot(site,pars):
        return pars.EL



    graphene, Subs = _make_graphene_lattice()

    sys = kwant.Builder()
    sys[[a(nx,ny) for a in Subs for nx in range(Nx) for ny in range(Ny)]] = potential
    sys[graphene.neighbors(eps=0.0001)] = hopping

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[a(0,ny) for a in Subs  for ny in range(Ny)] ]= left_pot
    lead_left[graphene.neighbors(eps=0.0001)] = -t
    lead_right[[a(Nx,ny) for a in Subs  for ny in range(Ny)] ]= right_pot
    lead_right[graphene.neighbors(eps=0.0001)] = -t

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)

    return sys.finalized(),lead_left,lead_right



def make_monolayer_graphene_system_1(Nx,Ny):
    """
    :param 长度Nx:
    :param 宽度Ny:
    :return sys without lead , lead_left, lead_right:
    """
    t = 1.0


    def potential(site,pars):
        x, y = site.pos
        EL, ER = pars.EL, pars.ER
        return EL+(ER-EL)*x/Nx + pars.W*(2.0*random.random()-1.0)


    def hopping(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))



    graphene = kwant.lattice.honeycomb()
    Subs = graphene.sublattices


    sys = kwant.Builder()
    sys[[a(nx,ny) for a in Subs for nx in range(Nx) for ny in range(Ny)]] = potential
    sys[graphene.neighbors(eps=0.0001)] = hopping

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[a(nx,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_left[graphene.neighbors(eps=0.0001)] = -t
    lead_right[[a(nx,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_right[graphene.neighbors(eps=0.0001)] = -t

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)

    return sys.finalized(),lead_left,lead_right


def make_monolayer_graphene_system_2(Nx,Ny):
    """
    :param 长度Nx:
    :param 宽度Ny:
    :return sys without lead , lead_left, lead_right:
    """
    t = 1.0


    def potential(site,pars):
        x, y = site.pos
        EL, ER = pars.EL, pars.ER
        return EL+(ER-EL)*x/Nx + pars.W*(2.0*random.random()-1.0)


    def hopping(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))



    graphene, Subs = _make_graphene_lattice()
    a, b, c, d = Subs
    hoppings = [
                ((0, 0), b, a),
                ((1, 0), b, a),
                ((0, -1), d, a),
                ((0, 0), c, b),
                ((0, 0), d, c),
                ((-1, 0), d, c)
                ]

    sys = kwant.Builder()

    sys[[a(nx,ny) for a in Subs for nx in range(Nx) for ny in range(Ny)]] = potential
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = hopping

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[a(0,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_left[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -t
    lead_right[[a(Nx-1,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_right[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -t

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)

    return sys.finalized(),lead_left,lead_right

def make_monolayer_graphene_system_with_v(Nx,Ny):
    """
    :param 长度Nx:
    :param 宽度Ny:
    :return sys without lead , lead_left, lead_right:
    """
    t = 1.0

    def potential(site,pars):
        x, y = site.pos
        EL, ER = pars.EL, pars.ER
        return EL+(ER-EL)*x/Nx + pars.W*(2.0*random.random()-1.0)

    def left_pot(site,pars):
        return pars.EL

    def right_pot(site,pars):
        return pars.ER


    def hopping(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))

    def hopping_vx(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))*(x2-x1)*1.0j

    def hopping_vy(site1,site2,pars) :
        x1, y1 = site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))*(y2-y1)*1.0j

    graphene, subs = _make_graphene_lattice()
    a, b, c, d = subs
    hoppings = [
                ((0, 0), b, a),
                ((1, 0), b, a),
                ((0, -1), d, a),
                ((0, 0), c, b),
                ((0, 0), d, c),
                ((-1, 0), d, c)
                ]

    sys = kwant.Builder()
    sysvx = kwant.Builder()
    sysvy = kwant.Builder()

    sys[[sub(nx,ny) for sub in subs for nx in range(Nx) for ny in range(Ny)]] = potential
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = hopping

    sysvx[[sub(nx, ny) for sub in subs for nx in range(Nx) for ny in range(Ny)]] = 0.0
    sysvx[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = hopping_vx

    sysvy[[sub(nx, ny) for sub in subs for nx in range(Nx) for ny in range(Ny)]] = 0.0
    sysvy[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = hopping_vy

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[sub(0,ny) for sub in subs  for ny in range(Ny)] ]= left_pot
    lead_left[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -t
    lead_right[[sub(Nx-1,ny) for sub in subs  for ny in range(Ny)] ]= right_pot
    lead_right[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -t

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)

    return list(sys.hoppings()),sys.finalized(), sysvx.finalized(), sysvy.finalized(), lead_left, lead_right


def test():
    Nx = 3; Ny = 3;
    hopping_sites, sys, sysvx, sysvy, lead_left, lead_right = make_monolayer_graphene_system_with_v(Nx,Ny)
    # kwant.plot(sys)
    # kwant.plot(sysvx)
    for site in hopping_sites:
        print site





if __name__ == "__main__":
    test()