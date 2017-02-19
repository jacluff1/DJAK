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
        self.av = np.average(self.val)
        self.std = np.std(self.val)
        
        def prec(self):
            p = abs(int(np.log10(abs(self.err)))) + 1
            return p

        def int_prec(x,p):
            x = str(x)
            x_keep = x[:-p]
            x_toss = x[-p:]
            x_app = []
            for x in x_toss:
                x_app.append('0')
            x_app = ''.join(x_app)
            x_new = int(x_keep+x_app)
            return x_new
        
        try:
            if self.err < 1:
                p = prec(self)
                self.p = p
                self.pval = round(self.val,p)
                self.perr = round(self.err,p)
                self.pav = round(self.av,p)
                self.pstd = round(self.std,p)
            elif self.err >= 1:
                self.perr = int(self.err)
                p = len(str(self.perr))
                self.p = p
                if p == 1:
                    self.pval = int(self.val)
                    self.pav = int(self.av)
                    self.pstd = int(self.std)
                else:
                    self.pval = int_prec(self.val,p)
                    self.pav = int_prec(self.av,p)
                    self.pstd = int_prec(self.av,p)

        except (TypeError,ValueError):
            try:
                self.pval = np.zeros(len(self.val))
                self.perr = np.zeros(len(self.err))
                self.p = np.zeros(len(self.err))
                self.av_err = np.average(self.err)
                if self.av_err < 1:
                	p_av = int(abs(int(np.log10(abs(self.av_err)))) + 1)
                	self.p_av = p_av
                	self.pav = round(self.av,p_av)
                	self.pav_err = round(self.av_err,p_av)
                	self.pstd = round(self.std,p_av)
                elif self.av_err >= 1:
                	self.pav_err = int(self.av_err)
                	p_av = len(str(self.pav_err))
                	self.p_av = p_av
                	if p_av == 1:
                		self.pav = int(self.av)
                		self.pstd = int(self.std)
                	else:
                		self.pav = int_prec(self.av,p_av)
                		self.pstd = int_prec(self.std,p_av)
                for i in range(len(self.pval)):
                	if self.err[i] < 1:
                		p = int(abs(int(np.log10(abs(self.err[i])))) + 1)
                		self.p[i] = p
                		self.pval[i] = round(self.val[i],p)
                		self.perr[i] = round(self.err[i],p)
                	elif self.err[i] >= 1:
                		self.perr[i] = int(self.err[i])
                		p = len(str(self.perr[i]))
                		self.p[i] = p
                		if p == 1:
                			self.pval[i] = int(self.val[i])
                			self.perr[i] = int(self.err[i])
                		else:
                			self.pval[i] = int_prec(self.val[i],p)
                			self.perr[i] = int_prec(self.err[i],p)
            except TypeError:
            	print("try something else")
        
def printvar(var):
    print(var.name,var.av,var.err,var.units)

def printerr(var):
    print(var.name,var.err,var.err/var.val)

