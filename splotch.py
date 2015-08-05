import numpy as np
from bokeh.plotting import figure
from bokeh.models import HoverTool 

TOOLS="hover,pan,wheel_zoom,box_zoom,reset,save"
minsize=200

def splotch(f, x, y, r, basesize, buffer, mc, key, bedge, ledge, nprot, c, m, l, uselegend):
    (prot1,prot2)=key.split('-')
    if f==None:
        f = figure(x_range=(0,mc[0]+buffer),y_range=(0,mc[1]+buffer),tools=TOOLS, plot_width=max(800/nprot,minsize),plot_height=max(800/nprot,minsize),title=None,min_border=10)
        f.xaxis.major_label_orientation = 3.14/4
    if not uselegend:
        l = None
    else:
        l = l.rstrip('/').split('/')[-1]
    s1 = f.scatter(x,y,size=np.sqrt(basesize/nprot*r),alpha=0.25,color=c,marker=m,legend=l) 
    s1.select(dict(type=HoverTool)).tooltips = {"("+prot1+","+prot2+")":"(@x,@y)"}
    if bedge:
        f.xaxis.axis_label = prot1
    if ledge:
        f.yaxis.axis_label = prot2
    f.legend.orientation = "top_left"
    return f

def makelegend(f,nprot,c,m,l):
#makes empty plot and sticks a legend in it
    if f==None:
        f = figure(tools=TOOLS, plot_width=150,plot_height=max(800/nprot,minsize),min_border=10)
    l = l.rstrip('/').split('/')[-1]
    f.scatter([],[],color=c,marker=m,alpha=0.4,legend=l)
    f.outline_line_color = None
    return f
