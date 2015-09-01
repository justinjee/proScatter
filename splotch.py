import numpy as np
import math
from bokeh.plotting import * 
from bokeh.models import HoverTool, BoxZoomTool, ResetTool, PanTool, WheelZoomTool 

MINSIZE=150
MAXSIZE=600
POINTSIZE=10
AXISSIZE=50
colors = ['blue','red','green']
markers = ['o','o','+']

def splotch_df(df, prot1, prot2, xr, yr, gridsize=1, x_on=True, y_on=True):

    hover = HoverTool(
        tooltips="""
        <div>
            <span style="font-size: 12px">("""+prot1+""" @x, """+prot2+""" @y)</span>
        </div>
        """
    )

    TOOLS=[hover, BoxZoomTool(), ResetTool(), PanTool(), WheelZoomTool()]

    fig = figure(tools=TOOLS, 
                 width=max(MINSIZE,MAXSIZE/gridsize)+int(y_on)*AXISSIZE, 
                 height=max(MINSIZE,MAXSIZE/gridsize)+int(x_on)*AXISSIZE, 
                 x_range=xr, y_range=yr, title=None, min_border=10)
    for i, exp in enumerate(set(df['label'])):
        exp_indices = df['label']==exp
        scatter = fig.scatter(df[exp_indices]['res1'], df[exp_indices]['res2'],
                size=np.sqrt(POINTSIZE * df[exp_indices]['SpecNum']),
                alpha=0.2, color=colors[i], marker=markers[i])
    fig.xaxis.axis_label = prot1
    fig.yaxis.axis_label = prot2
    fig.xaxis.visible = x_on 
    fig.yaxis.visible = y_on
    fig.xaxis.major_label_orientation = math.pi / 4

    return fig

def stripfolder(s):
    return s.rstrip('/').split('/')[-1]

def makelegend(df):
    '''
    makes empty plot and sticks a legend in it
    '''
    f = figure(plot_width=MINSIZE, plot_height=MINSIZE, min_border=10)
    for i, exp in enumerate(set(df['label'])):
        l = stripfolder(exp)
        f.scatter([], [], color=colors[i], marker=markers[i], alpha=0.4, legend=l)
        f.outline_line_color = None
    return f
