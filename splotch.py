import numpy as np
from bokeh.plotting import figure

def splotch(f, x, y, r, basesize, buffer, mc, key, bedge, ledge, nprot, c, m):
    (prot1,prot2)=key.split('-')
    TOOLS="pan,wheel_zoom,box_zoom,reset,save"
    if f==None:
        f = figure(x_range=(0,mc[0]+buffer),y_range=(0,mc[1]+buffer),tools=TOOLS, plot_width=1000/nprot,plot_height=1000/nprot,title=None,min_border=10)
        f.xaxis.major_label_orientation = 3.14/4
        #f.xaxis.minor_tick_line_color = None
        #f.yaxis.minor_tick_line_color = None
    f.scatter(x,y,size=np.sqrt(basesize/nprot*r),alpha=0.4,color=c,marker=m) 
    if bedge:
        f.xaxis.axis_label = prot1
    if ledge:
        f.yaxis.axis_label = prot2
    return f
