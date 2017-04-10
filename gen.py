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

def float_prec(x):
    str_x = str(x)
    find_e = str.find(str_x,'e-')
    if find_e != -1:
        p = float(str_x[find_e+2:])
    else:
        str_x = str_x[2:]
        p = 1
        for i in range(len(str_x)):
            if str_x[i] == '0':
                p += 1
            else:
                break
    return p

def int_prec(x,p,err=False):
    if p == 1:
        err=True

    if err == False:
        mult = 10**(p)
        x_new = int(round(x/mult)) * mult
    else:
        mult = 10**(p-1)
        x_new = int(round(x/mult)) * mult

    return x_new

class var:
    def __init__(self,val,err):
        self.val = float(val)
        self.err = float(err)

        try: # try for scalar values in val and err

            if all(( self.err < 1, self.err != 0.0 )):

                p = int(float_prec(self.err))
                self.p = p
                self.pval = round(self.val,p)
                self.perr = round(self.err,p)

            elif self.err == 0.0:
                print("needs work")

                if self.val < 1:
                    p = int(float_prec(self.val))
                    self.p = p
                    self.pval = round(self.val,p)
                    self.perr = .5 * 10**-p
                else:
                    p = len(str(int(self.val)))
                    self.p = p
                    self.pval = int_prec(self.val,p)
                    self.perr = .5 * 10**-p

            elif self.err > 1:
                p = len(str(int(self.err)))
                self.p = p
                self.pval = int_prec(self.val,p)
                self.perr = int_prec(self.err,p,err=True)

        except:
            print("made into an array")

            try: # try for same length array/list in val and err

                self.val = np.array(self.val)
                self.err = np.array(self.err)
                self.p = np.zeros_like(self.err)
                self.pval = np.zeros_like(self.err)
                self.perr = np.zeros_like(self.err)

                for i in range(len(self.err)):

                    if self.err[i] < 1:

                        p = float_prec(self.err[i])
                        self.p[i] = p
                        self.pval[i] = round(self.val[i],p)
                        self.perr[i] = round(self.err[i],p)


                    else:

                        p = len(str(int(self.err[i])))
                        self.p[i] = p
                        self.pval[i] = int_prec(self.val[i],p)
                        self.perr[i] = int_prec(self.err[i],p,err=True)

                self.p = np.array(self.p).astype(int)

            except:
                print("try something else")
