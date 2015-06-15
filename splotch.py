import numpy as np
import matplotlib.pyplot as plt
from pylab import *

def splotch(ax, x, y, r, basesize, buffer, mc, key, ledge, bedge, nprot, color, m):
        (prot1,prot2)=key.split('-')
        ax.scatter(x,y,s=min(20,int(basesize/nprot))*r,alpha=0.3,c=color,edgecolor=color,marker=m,picker=1)
        ax.set_xlim(0, mc[0]+buffer)
        ax.set_ylim(0, mc[1]+buffer)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        [(ax.xaxis.get_major_ticks()[t].label).set_rotation('vertical') for t in range(len(ax.xaxis.get_major_ticks()))]
        tickfontsize = max(3,10-int(nprot/2))
        [(ax.xaxis.get_major_ticks()[t].label).set_fontsize(tickfontsize) for t in range(len(ax.xaxis.get_major_ticks()))]
        [(ax.yaxis.get_major_ticks()[t].label).set_fontsize(tickfontsize) for t in range(len(ax.yaxis.get_major_ticks()))]
        if nprot>1:
            ax.get_yaxis().set_label_coords(-0.4-.01*nprot,0.5)
            ax.get_xaxis().set_label_coords(0.5,-0.45-.03*nprot)
        labelfontsize = max(18-nprot,4)
        if ledge:
            ax.set_ylabel(prot2,fontsize=labelfontsize)
        else:
            ax.get_yaxis().set_ticks([])
        if bedge:
            ax.set_xlabel(prot1,fontsize=labelfontsize)
        else:
            ax.get_xaxis().set_ticks([])
