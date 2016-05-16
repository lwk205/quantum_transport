# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-12
Contact: gengyusheng@gmail.com
内容： 单层黑磷的体系构建
"""

from __future__ import division
import kwant
import random
from cmath import exp
from lattices import  _make_lattice_bp_single_layer
import matplotlib.pylab


def make_system_two_terminal( Nx=10,Ny=10):
    """
    make_sys
    pars 为 pars.EL, pars.ER, pars.phi, pars.W
    """
    lattice, sublattices = _make_lattice_bp_single_layer()
    A = sublattices

    def potential( site, pars):
        EL, ER, W,vtop,vbottom,vi = pars.EL, pars.ER, pars.W,pars.vtop,\
                                    pars.vbottom, pars.vi
        x0,y0=site.pos
        PNU = EL+(ER-EL)*x0/Nx+W*(2.0*random.random()-1.0)
        
        x,y=site.tag
        if site.family==A[0]:
            if y==Ny-1:
                return vtop+PNU
            else:
                return vi+PNU
        if site.family==A[1]:
            if y==Ny-1:
                return vtop+PNU
            else:
                return vi+PNU
            
        if site.family==A[2]:
            if y==0:
                return vbottom+PNU
            else:
                return -vi+PNU
        if  site.family==A[3]:
            if y==0:
                return vbottom+PNU
            else:
                return -vi+PNU
            


    def left_pot(site, pars):
        return pars.EL


    def ringht_pot( site, pars):
        return pars.ER


    def set_hopping(sys):
        t1, t2, t3, t4, t5 = -1.0 ,3.0041, -0.1680, -0.0861, -0.0451
        

        #t1: nearest hopping
        def hop_t1(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t1*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[0],A[2])]=hop_t1
        sys[kwant.builder.HoppingKind((0,0),A[3],A[1])]=hop_t1
        sys[kwant.builder.HoppingKind((1,0),A[0],A[2])]=hop_t1
        sys[kwant.builder.HoppingKind((1,0),A[3],A[1])]=hop_t1

        #t2: 2nd nearest hopping
        def hop_t2(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t2*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[1],A[2])]=hop_t2
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[3])]=hop_t2

        #t3: 3st nearest hopping
        def hop_t3(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t3*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,1),A[3],A[1])]=hop_t3
        sys[kwant.builder.HoppingKind((1,1),A[3],A[1])]=hop_t3
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[2])]=hop_t3
        sys[kwant.builder.HoppingKind((1,-1),A[0],A[2])]=hop_t3

        #t4: 4th nearest hopping
        def hop_t4(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t4*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,1),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((1,1),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((0,0),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((1,0),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((1,-1),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((0,0),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((1,0),A[0],A[1])]=hop_t4

        #t5: 5th nearest hopping
        def hop_t5(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t5*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[0],A[3])]=hop_t5
        sys[kwant.builder.HoppingKind((0,-1),A[2],A[1])]=hop_t5

        #return sys


    sys = kwant.Builder()
    sys[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= potential
    set_hopping(sys)

    sym_left = kwant.lattice.TranslationalSymmetry(lattice.vec((-1,0)))
    sym_right = kwant.lattice.TranslationalSymmetry(lattice.vec((1,0)))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_left[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= potential
    set_hopping(lead_left)

    lead_right[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= potential
    set_hopping(lead_right)

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)

    return sys.finalized(), lead_left

def make_system_four_terminal( Nx=10,Ny=10):
    """
    make_sys
    pars 为 pars.EL, pars.ER, pars.phi, pars.W
    """
    lattice, sublattices = _make_lattice_bp_single_layer()


    def potential( site, pars):
        EL, ER, W = pars.EL, pars.ER, pars.W
        x,y=site.pos
        return EL+(ER-EL)*x/Nx+W*(2.0*random.random()-1.0)


    def left_pot(site, pars):
        return pars.EL


    def ringht_pot( site, pars):
        return pars.ER

    def up_pot(site,pars):
        return pars.EU

    def down_pot(site,pars):
        return pars.ED



    def set_hopping(sys):
        t1, t2, t3, t4, t5 = -1.0 ,3.0041, -0.1680, -0.0861, -0.0451
        A = sublattices

        #t1: nearest hopping
        def hop_t1(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t1*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[0],A[2])]=hop_t1
        sys[kwant.builder.HoppingKind((0,0),A[3],A[1])]=hop_t1
        sys[kwant.builder.HoppingKind((1,0),A[0],A[2])]=hop_t1
        sys[kwant.builder.HoppingKind((1,0),A[3],A[1])]=hop_t1

        #t2: 2nd nearest hopping
        def hop_t2(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t2*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[1],A[2])]=hop_t2
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[3])]=hop_t2

        #t3: 3st nearest hopping
        def hop_t3(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t3*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,1),A[3],A[1])]=hop_t3
        sys[kwant.builder.HoppingKind((1,1),A[3],A[1])]=hop_t3
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[2])]=hop_t3
        sys[kwant.builder.HoppingKind((1,-1),A[0],A[2])]=hop_t3

        #t4: 4th nearest hopping
        def hop_t4(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t4*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,1),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((1,1),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((0,0),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((1,0),A[3],A[2])]=hop_t4
        sys[kwant.builder.HoppingKind((0,-1),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((1,-1),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((0,0),A[0],A[1])]=hop_t4
        sys[kwant.builder.HoppingKind((1,0),A[0],A[1])]=hop_t4

        #t5: 5th nearest hopping
        def hop_t5(site1,site2,pars):
            phi = pars.phi
            x1, y1 = site1.pos
            x2, y2 = site2.pos
            return t5*exp(-0.5j*phi*(x1-x2)*(y1+y2))

        sys[kwant.builder.HoppingKind((0,0),A[0],A[3])]=hop_t5
        sys[kwant.builder.HoppingKind((0,-1),A[2],A[1])]=hop_t5

        #return sys


    sys = kwant.Builder()
    sys[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= potential
    set_hopping(sys)

    sym_left = kwant.lattice.TranslationalSymmetry(lattice.vec((-1,0)))
    sym_right = kwant.lattice.TranslationalSymmetry(lattice.vec((1,0)))
    sym_up = kwant.lattice.TranslationalSymmetry(lattice.vec((0, 1)))
    sym_down = kwant.lattice.TranslationalSymmetry(lattice.vec((0, -1)))
    lead_left = kwant.Builder(sym_left)
    lead_right = kwant.builder.Builder(sym_right)
    lead_up = kwant.Builder(sym_up)
    lead_down = kwant.Builder(sym_down)

    lead_left[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= left_pot
    set_hopping(lead_left)

    lead_right[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= ringht_pot
    set_hopping(lead_right)

    lead_up[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= up_pot
    set_hopping(lead_up)

    lead_down[[a(nx,ny) for a in sublattices for nx in range(Nx) for ny in range(Ny)] ]= down_pot
    set_hopping(lead_down)

    sys.attach_lead(lead_left)
    sys.attach_lead(lead_right)
    sys.attach_lead(lead_down)
    sys.attach_lead(lead_up)

    return sys.finalized(), lead_up

def test():
    sys, lead = make_system_two_terminal(Nx =10, Ny = 10)
    sys1, lead1 = make_system_four_terminal(Nx = 5, Ny = 5)
    kwant.plot(sys1)
    kwant.plot(sys)


if __name__ == '__main__':
    test()