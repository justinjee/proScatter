import numpy as np
import math
from bokeh.plotting import * 
from bokeh.models import HoverTool, BoxZoomTool, ResetTool, PanTool, WheelZoomTool 

MINSIZE=200
MAXSIZE=800
POINTSIZE=10
def splotch_df(df, prot1, prot2, xr, yr, x_on=True, y_on=True):

    hover = HoverTool(
        tooltips="""
        <div>
            <span style="font-size: 12px">("""+prot1+""" @x, """+prot2+""" @y)</span>
        </div>
        """
    )

    TOOLS=[hover, BoxZoomTool(), ResetTool(), PanTool(), WheelZoomTool()]

    fig = figure(tools=TOOLS, width=250, height=250, x_range=xr, y_range=yr, title=None, min_border=10)
    scatter = fig.scatter(df['res1'], df['res2'],
                size=np.sqrt(POINTSIZE * df['SpecNum']),
                alpha=0.25)
    fig.xaxis.axis_label = prot1
    fig.yaxis.axis_label = prot2
    fig.xaxis.visible = x_on 
    fig.yaxis.visible = y_on
    fig.xaxis.major_label_orientation = math.pi / 4

    return fig

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
