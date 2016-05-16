# -*- coding: utf-8 -*-
"""
Created by Rainson on 2016-05-12
Contact: gengyusheng@gmail.com
内容： 作图的函数
"""

#import matplotlib as mpl; mpl.use("Agg")
from __future__ import division
from matplotlib import pylab as plt
import kwant

def simple_plot_data_2d(data,ax=None, xlabel = None, ylabel=None, xlim= None, ylim= None):
    """
    简单的作图程序，主要用于直观显示数据，具体的参数调节需要编写复杂的作图程序
    :param data:
    :param xlim:
    :param ylim:
    :return:
    """
    if ax is None:
        ax = plt.subplot(111)

    ax.plot(data[0],data[1])

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    return  ax

def simple_plot_data_2d_without_ax(data, xlabel = None, ylabel=None, xlim= None, ylim= None):
    """
    简单的作图程序，主要用于直观显示数据，具体的参数调节需要编写复杂的作图程序
    :param data:
    :param xlim:
    :param ylim:
    :return:
    """
    plt.figure()
    plt.plot(data[0],data[1])

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    plt.show()

def plot_3d_map(data):
    """
    二维map函数，可用于画相图之类的东西
    :param data:
    :return:
    """


def test():
    data = [range(10), range(10, 20)]
    data1 = [range(10), range(10)]
    ax = simple_plot_data_2d(data,ax=None,xlabel='x',ylabel='y',xlim=[0,10],ylim=None)
    ax = simple_plot_data_2d(data1, ax, 'x', 'y', [0, 10], ylim=None)

    ax.figure.show()
    return ax


if __name__ == "__main__":
    ax = test()