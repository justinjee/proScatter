import numpy as np
import math
from bokeh.plotting import * 
from bokeh.models import HoverTool, BoxZoomTool, ResetTool, PanTool, WheelZoomTool 

#TOOLS="hover,pan,wheel_zoom,box_zoom,reset,save"
MINSIZE=200
MAXSIZE=800

def splotch_df(df, basesize=50):
    TOOLS = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,previewsave"
    rows = []
    row = []
    num_prot = len(set(df['prot1']))
    grouped = df.groupby(['prot1', 'prot2'])

    for i, (key,group) in enumerate(grouped):
        fig = figure(tools=TOOLS, width=250, height=250, title=None, min_border=10)
        scatter = fig.scatter(group['res2'], group['res1'],
                    size=np.sqrt(basesize / num_prot * group['SpecNum']),
                    alpha=0.25)
        fig.xaxis.axis_label = key[1]
        fig.yaxis.axis_label = key[0]
        fig.xaxis.visible = False
        fig.yaxis.visible = False
        fig.xaxis.major_label_orientation = math.pi / 4
        row.append(fig)
        if (i + 1) % num_prot == 0:
            fig.yaxis.visible = True
            rows.append(row[::-1])
            row = []

    for fig in rows[-1]:
        fig.xaxis.visible = True

    return gridplot(rows)

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
