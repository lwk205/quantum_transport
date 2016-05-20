# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-12
Contact: gengyusheng@gmail.com
内容： 单层黑磷的体系构建
"""

from __future__ import division
import kwant
from math import sqrt


def _make_lattice_bp_single_layer():
    """
    :return:
    """
    a=3.297399977
#    b=11.5564002991
    c=4.5771999359
    BP2D = kwant.lattice.general([(a, 0.0),(0.0,c)],\
                 [(0.0*a,0.912069976*c),\
                 (0.5*a,0.412070006*c), \
                 (0.5*a,0.587930024*c),\
                 (0.0*a,0.087930001*c)])
    subs = BP2D.sublattices
    return BP2D, subs


def _make_graphene_lattice() :

    graphene = kwant.lattice.general(
                                    [(1, 0), (0.0, sqrt(3))],
                                    [
                                        (0.5, -0.5 / sqrt(3)), (0, 0),
                                        (0.0, 1 / sqrt(3)), (0.5, 0.5 * sqrt(3))
                                    ]
                                      )

    Subs = graphene.sublattices
    return graphene, Subs