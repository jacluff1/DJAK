import numpy as np
import csv

__all__ = ['CSV','sort','nearest','maxima','minima','extrema',
           'var']

def CSV(path_in,path_out,title,numpyA=float,saveA=True): # imports csv, converts to np array, and saves it as .npy file
    #path_in = str(path_in)
    #path_out = str(path_out)
    #title = str(title)
    DATA = []
    with open(path_in+title+".csv", newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    DATA = np.array(DATA).astype(numpyA)
    if saveA == True:
        np.save(path_out+title+".npy",DATA)
    else:
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
    def __init__(self,name,val,err,units):
        self.val = val
        self.err = err
        self.units = str(units)
        self.name = str(name)
        
        def prec(self):
            p = abs(int(np.log10(abs(self.err)))) + 1
            return p
        
        if type(self.val) == float or type(self.val) == np.float64:
            p = prec(self)
            self.av = self.val
            self.pval = round(self.val,p)
            self.perr = round(self.err,p)
            self.pav = self.pval
        
        elif type(self.err) == float:
            self.val = np.array(self.val)
            self.av = np.average(self.val)
            p = prec(self)
            self.perr = round(self.err,p)
            self.pval = np.round(self.val,p)
            self.pav = np.round(self.av,p)
        
        elif self.err == 'std':
            self.val = np.array(self.val)
            self.err = np.std(self.val)
            self.av = np.average(self.val)
            p = prec(self)
            self.perr = round(self.err,p)
            self.pval = np.round(self.val,p)
            self.pav = round(self.av,p)
        else:
            print("need additional boolean conditions")
            print(self.name,self.val,self.err)
            self.pval = val
            self.perr = err
            self.pav = self.pval
        
def printvar(var):
    print(var.name,var.av,var.err,var.units)

def printerr(var):
    print(var.name,var.err,var.err/var.val)

