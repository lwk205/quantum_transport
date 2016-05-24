# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-12
Contact: gengyusheng@gmail.com
内容： 各种计算的函数
"""

from __future__ import division
import kwant
import numpy as np
import scipy.linalg as sl
import  scipy.sparse.linalg as sla



def calculate_conductance(sys,Ef, out_leads, in_leads, pars):
    """
    :param
            sys : 利用kwant产生的体系
            EF : 费米能的大小
            out_leads : 电子出射的导线 序列或者是一个整数
            in_leads : 入射的导线 序列或者是一个整数
    :return
            T: kwant.solvers.common.SMatrix ，如果只有单个入射出射导线的话就是一个实数
    """
    smatrix = kwant.smatrix(sys, Ef, args=[pars] )
    T = smatrix.transmission(out_leads, in_leads)
    return T


def caculate_energy_band(lead, pars):
    """
    :param lead:
    :param pars:
    :return data:
    """
    from math import pi
    bands=kwant.physics.Bands(lead.finalized(),args = [pars])
    N = 100
    ks=np.linspace(-pi,pi,N)
    enbands=[list(bands(k)) for k in ks]
    return [ks, enbands]

def wave_density(sys,  pars, lead_nr=0):
    wf = kwant.wave_function(sys, pars.EF, args=[pars])
    return (abs(wf(lead_nr))**2).sum(axis=0)


def get_system_hamiltonian(sys, pars, sparse = False):
    ham_mat = sys.hamiltonian_submatrix(args=[pars], sparse=sparse)
    return ham_mat

def system_eigvalues(sys, pars,sparse=False):
    ham=get_system_hamiltonian(sys,pars,sparse=sparse)
    if sparse:
        return sla.eigsh(ham,k=15,which='SM',return_eigenvectors=False)
    else:
        return sl.eigh(ham,eigvals_only=True)


def test():
    pass

if __name__ == "__main__":
    test()