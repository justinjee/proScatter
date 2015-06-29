import numpy as np
from bokeh.plotting import figure
from bokeh.models import HoverTool 

def splotch(f, x, y, r, basesize, buffer, mc, key, bedge, ledge, nprot, c, m, l, uselegend):
    (prot1,prot2)=key.split('-')
    TOOLS="hover,pan,wheel_zoom,box_zoom,reset,save"
    if f==None:
        f = figure(x_range=(0,mc[0]+buffer),y_range=(0,mc[1]+buffer),tools=TOOLS, plot_width=1000/nprot,plot_height=1000/nprot,title=None,min_border=10)
        f.xaxis.major_label_orientation = 3.14/4
    if not uselegend:
        l = None
    else:
        l = l.rstrip('/').split('/')[-1]
    s1 = f.scatter(x,y,size=np.sqrt(basesize/nprot*r),alpha=0.4,color=c,marker=m,legend=l) 
    s1.select(dict(type=HoverTool)).tooltips = {"("+prot1+","+prot2+")":"(@x,@y)"}
    if bedge:
        f.xaxis.axis_label = prot1
    if ledge:
        f.yaxis.axis_label = prot2
    f.legend.orientation = "top_left"
    return f
