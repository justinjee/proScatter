import numpy as np
from bokeh.plotting import * 
from bokeh.models import HoverTool, BoxZoomTool, ResetTool, PanTool, WheelZoomTool 

#TOOLS="hover,pan,wheel_zoom,box_zoom,reset,save"
minsize=200
maxsize=800

def splotch(f, fkey, x, y, r, basesize, buffer, mc, key, bedge, ledge, nprot, c, m, l, uselegend, unjoin=False):
    (prot1, prot2) = key.split('-')
    hover = HoverTool(
        tooltips="""
        <div>
            <span style="font-size: 12px">("""+prot1+""" @x, """+prot2+""" @y)</span>
        </div>
        """
    )
    TOOLS = [hover, BoxZoomTool(), ResetTool(), PanTool(), WheelZoomTool()]
    if not f:
        if fkey and not unjoin:
            f = figure(x_range=fkey.x_range, y_range=fkey.y_range, tools=TOOLS,
                    plot_width=max(maxsize/nprot, minsize),
                    plot_height=max(maxsize/nprot, minsize),
                    title=None, min_border=10)
            f.xgrid.bounds = (0, mc[0])
            f.ygrid.bounds = (0, mc[1])
        else:
            f = figure(x_range=(0, mc[0]+buffer), y_range=(0, mc[1]+buffer), tools=TOOLS, \
                    plot_width=max(maxsize/nprot, minsize),
                    plot_height=max(maxsize/nprot, minsize),
                    title=None,min_border=10)
        f.xaxis.major_label_orientation = 3.14/4
    if not uselegend:
        l = None
    else:
        l = l.rstrip('/').split('/')[-1]
    s1 = f.scatter(x, y, size=np.sqrt(basesize/nprot*r), alpha=0.25, color=c, marker=m, legend=l) 
    #s1.select(dict(type=HoverTool)).tooltips = {"("+prot1+","+prot2+")":"(@x,@y)"}
    if bedge:
        f.xaxis.axis_label = prot1
        f.plot_height += 11
    else:
        f.xaxis.visible = None
    if ledge:
        f.yaxis.axis_label = prot2
        f.plot_width += 11
    else:
        f.yaxis.visible = None 
    f.legend.orientation = "top_left"
    return f

def stripfolder(s):
    return s.rstrip('/').split('/')[-1]

def makelegend(f, nprot, c, m, l):
    '''
    makes empty plot and sticks a legend in it
    '''
    if f is None:
        f = figure(plot_width=150, plot_height=max(maxsize/nprot, minsize), min_border=10)
    l = stripfolder(l)
    f.scatter([], [], color=c, marker=m, alpha=0.4, legend=l)
    f.outline_line_color = None
    return f
