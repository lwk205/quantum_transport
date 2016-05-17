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
        x, y = site.tag
        EL, ER = pars.EL, pars.ER
        return EL+(ER-EL)*x/(Nx-1) + pars.W*(2.0*random.random()-1.0)

    def hopping(site1, site2, pars):
        x1, y1 =site1.pos
        x2, y2 = site2.pos
        phi = pars.phi
        return -t*exp(-0.5j*phi*(x1-x2)*(y1+y2))


    graphene, Subs = _make_graphene_lattice()

    sys = kwant.Builder()
    sys[[a(nx,ny) for a in Subs for nx in range(Nx) for ny in range(Ny)]] = potential
    sys[graphene.neighbors(eps=0.0001)] = hopping

    sym_left = kwant.lattice.TranslationalSymmetry((-1,0))
    sym_right = kwant.lattice.TranslationalSymmetry((1,0))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[a(nx,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_left[graphene.neighbors(eps=0.0001)] = hopping
    lead_right[[a(nx,ny) for a in Subs  for ny in range(Ny)] ]= potential
    lead_right[graphene.neighbors(eps=0.0001)] = hopping

    return sys,lead_left,lead_right

def test():
    Nx = 10; Ny = 10;
    sys,lead0,lead1 = make_monolayer_graphene_system(Nx,Ny)
    kwant.plot(sys)
    sys.attach_lead(lead0)
    sys.attach_lead(lead1)
    kwant.plot(sys)

if __name__ == "__main__":
    test()