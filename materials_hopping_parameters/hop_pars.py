# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-12
Contact: gengyusheng@gmail.com
内容： 不同材料的hopping系数
"""

def monolayer_black_phosphorus():
    """
    单层黑磷的hopping系数
    cite from: prb,89,201408(R)(2014)
    """
    t1, t2, t3, t4, t5 = -1.0, 3.0041, -0.1680, -0.0861, -0.0451
    return [t1, t2, t3, t4, t5]

def bilayer_black_phosphorus():
    """
    单层黑磷的hopping系数
    cite from: prb,89,201408(R)(2014)
    """
    pass