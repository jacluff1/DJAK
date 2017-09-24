import numpy as np
import pandas as pd
import csv

__all__ = ['CSV','sort','nearest','maxima','minima','extrema',
           'var']

def CSV(title,skip_row=0,path_in='../csv/',path_out='../npy/',numpyA=float,saveA=True):
    """ imports csv, converts to np array, and saves it as .npy file"""
    #path_in = str(path_in)
    #path_out = str(path_out)
    #title = str(title)

    DATA = []
    with open(path_in+title, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)

    DATA = np.array(DATA[skip_row:]).astype(numpyA)
    if saveA == True:
        np.save(path_out+title+".npy",DATA)

    return DATA

def sort(Array,column): # sorts array by given column number
    return Array[Array[:,column].argsort()]

def nearest(Array1,value): # finds the index of the nearest element in an array to a given value
    index = (np.abs(Array1-value)).argmin()
    return index

def maxima(X,Y,x_approx,dx,ret='i_max'):# finds local maximum of array
    dx = int(dx)
    i_approx = nearest(X,x_approx) # good
    i_add = len(X[:i_approx-dx])
    X_n = np.array(X[i_approx-dx:i_approx+dx+1])
    Y_n = np.array(Y[i_approx-dx:i_approx+dx+1])
    y_max = np.max(Y_n) # good
    i_max = nearest(Y_n,y_max) # good
    x_max = X_n[i_max] # good
    if ret == 'i_max':
        return i_max + i_add
    if ret == 'y_max':
        return y_max
    if ret == 'x_max':
        return x_max
    if ret == 'all':
        return (i_max+i_add),x_max,y_max

def minima(X,Y,x_approx,dx,ret='i_min'): # finds local minimum of an array
    dx = int(dx)
    i_approx = nearest(X,x_approx)
    i_add = len(X[:i_approx-dx])
    X_n = np.array(X[i_approx-dx:i_approx+dx+1])
    Y_n = np.array(Y[i_approx-dx:i_approx+dx+1])
    y_min = np.min(Y_n)
    i_min = nearest(Y_n,y_min)
    x_min = X_n[i_min]
    if ret == 'i_min':
        return i_min + i_add
    if ret == 'y_min':
        return y_min
    if ret == 'x_min':
        return x_min
    if ret == 'all':
        return (i_min+i_add),x_min,y_min

def extrema(X,Y,x_range,xmin_approx,xmax_approx,ret='sep',sortcol=0): # finds the exact extrema given a list of approximate minima and maxima
    dx = int(x_range)
    i_min,i_max = [],[]
    for x in xmin_approx:
        i_min.append(nearest(X,x))
    for x in xmax_approx:
        i_max.append(nearest(X,x))

    y_min,y_max = [],[]
    for i in i_min:
        y_min.append(np.min(Y[i-dx:i+dx+1]))
    for i in i_max:
        y_max.append(np.max(Y[i-dx:i+dx+1]))

    x_min,x_max = [],[]
    for y in y_min:
        x_min.append(X[nearest(Y,y)])
    for y in y_max:
        x_max.append(X[nearest(Y,y)])

    minima,maxima = [],[]
    for i in range(len(x_min)):
        minima.append([x_min[i],y_min[i]])
    for i in range(len(x_max)):
        maxima.append([x_max[i],y_max[i]])

    if ret == 'sep':
        return minima,maxima
    else:
        minima = np.array(minima)
        maxima= np.array(maxima)
        ext = sort(np.vstack((minima,maxima)),sortcol)
        return ext

class var:
    def __init__(self,val,err):
        self.val    =   float(val)
        self.err    =   float(err)
        ylen        =   len(str(int(self.err)))
        n           =   0

        if self.err > 1:
            self.err_scale  =   10**(ylen-1)
            y               =   int(round(self.err / self.err_scale))
            if self.err > 1000:
                self.perr   =   '%s $\\times$ 10$^{%s}$' % (y,(ylen-1))
            else:
                self.perr   =   str(y*self.err_scale)

        elif self.err < 1:
            y               =   self.err
            while y < 1:
                n           +=  1
                y           *=  10
            self.err_scale  =   10**(-n)
            y               =   int(round(self.err / self.err_scale))
            if self.err < .001:
                self.perr   =   '%s $\\times$ 10$^{-%s}$' % (y,n)
            else:
                self.perr   =   str(y * self.err_scale)

        if self.val > 1:
            x               =   int(round(self.val / self.err_scale))
            if self.err > 1000:
                self.pval   =   '%s $\\times$ 10$^{%s}$' % (x,(ylen-1))
            else:
                self.pval   =   str(x*self.err_scale)

        elif self.val < 1:
            x               =   int(round(self.val / self.err_scale))
            if self.err < .001:
                self.pval   =   '%s $\\times$ 10$^{-%s}$' % (x,n)
            else:
                self.pval   =   str(x * self.err_scale)

def directory_checker(dirname):
    """ check if a directory exists, makes it if it doesn't """
    dirname =   str(dirname)
    if not os.path.exists(dirname):
        try:
            os.mkdir(dirname)
        except:
            os.stat(dirname)
