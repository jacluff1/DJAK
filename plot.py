__all__ = ['data_set','note_set','Myax','MyPlot']

import matplotlib.pyplot as plt
import numpy as np

class data:
    def __init__(self,X,Y,style,color,label):
        self.X = X
        self.Y = Y
        self.style = str(style)
        self.color = str(color)
        self.label = str(label)

class note:
    def __init__(self,note,x,y,color):
        self.note = note
        self.x = x
        self.y = y
        self.color = str(color)
        self.fontsize = 20 

class ax:
    def __init__(self,datasets,subplot,xlabel,ylabel,title):
        self.datasets = datasets
        self.subplot = subplot
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        self.notes = []

def plot(axies,name=None,size=(15,7.5),fs=20):
    """ plot multiple axis with multiple plots/axis

    args
    ----
    axies: a list or array of 'ax' class objects
    name: ** if name is defined, figure will be saved and closed
    size: ** figure size

    returns
    -------
    either shows or saves figure

    axis class objects
    ------------------
    datasets: list or array of 'data' class objects
    subplot: plot number (examples: 121, 122, etc..) 
    xlabel: string to lable x-axis
    ylabel: string to lable y-axis
    title: string to title axis
    notes: optional 'note' class to annotate subplot

    note class object
    -----------------
    note: annotated string
    x: x coordinate
    y: y coordinate
    color: color string
    fontsize: ** size of annotation

    data
    ----
    X: list or array od x-axis data
    Y: list or array of y-axis data
    style: string, linestyle
    color: string, color
    label: string, legend label
    """
    
    fig = plt.figure(figsize=size)
    if len(axies) < 0:
        print("need list of 'axies'")
        return ValueError
    
    for a in axies:
        ax = plt.subplot(a.subplot)
        ax.set_xlabel(str(a.xlabel), fontsize=fs)
        ax.set_ylabel(str(a.ylabel), fontsize=fs)
        ax.set_title(str(a.title), fontsize=fs+2)
        
        xmin,xmax = 0,0
        for d in a.datasets:
            ax.plot(d.X,d.Y,d.style,color=d.color,label=d.label)
            if len(d.X) >= 2:
                if min(d.X) < xmin:
                    xmin = min(d.X)
                if max(d.X) > xmax:
                    xmax = max(d.X)
        ax.set_xlim(xmin,xmax)
        lgd = ax.legend(loc='upper center', bbox_to_anchor=(.5,-.11), numpoints=1, ncol=2, fancybox=True, shadow=True, fontsize=fs)
        
        if len(a.notes) > 0:
            for n in a.notes:
                ax.annotate(n.note, xy=(n.x,n.y), color=n.color, fontsize=n.fontsize)
        
        plt.tight_layout()
        
    if name == None:
        plt.show()
    else:
        fig.savefig(str(name)+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    
    return
