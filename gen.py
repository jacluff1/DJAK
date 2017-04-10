import numpy as np
import pandas as pd
import csv

__all__ = ['CSV','sort','nearest','maxima','minima','extrema',
           'var']

def csv_2_numpy(path_in,title,path_out='',dtype=np.float64,title_row=False):
    """ opens up a .csv file and converts it to nupmpy

    Parameters
    ----------
    path_in:    pathway to csv                          - str
    title:      name of file                            - str
    path_out:   ** if != '', will save to this pathway  - default = ''
    dtype:      ** data type of numpy array             - default = np.float64
    title_row:  ** if == True will skip a row           - default = False

    Returns
    -------
    np.array
    """

    if path_in[-1] != '/':                              path_in  += '/'
    if title[-4] != '.csv':                             title    += '.csv'
    if all(( path_out != '' , path_out[-1] != '/' )):   path_out += '/'

    DATA = []
    with open(path_in+title, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)

    if title_row:
        DATA = np.array(DATA[1:]).astype(dtype)
    else:
        DATA = np.array(DATA).astype(dtype)

    if path_out != '':
        np.save(path_out+title+".npy",DATA)

    return DATA

def csv_2_pandas(path_in,title,path_out='',dtype=np.float64,title_row=False,keys=[]):
    """ opens a .csv file and as a panda data frame

    Parameters
    ----------
    path_in:    pathway to csv                          - str
    title:      name of file                            - str
    path_out:   ** if != '', will save to this pathway  - default = ''
    dtype:      ** data type of numpy array             - default = np.float64
    title_row:  ** if == True will be column names      - default = False
    keys:       ** if len(keys)>0 will be column names  - default = []

    Returns
    -------
    pd.DataFrame
    """

    # check parameters
    assert any(( len(keys) > 0 , title_row )), "either use title row as names or use 'keys'"

    if path_in[-1] != '/':                              path_in     += '/'
    if title[-4:] != '.csv':                            title       += '.csv'
    if all(( path_out != '' , path_out[-1] != '/' )):   path_out    += '/'

    # determine column names
    if title_row == False:
        # assert len(keys) == len(DATA[0,:]), "length of keys must equal number of colums!"
        df = pd.read_csv(path_in+title,names=keys)

    else:
        df = pd.read_csv(path_in+title)

    if path_out != '':
        pd.to_pickle(df,path_out+title)

    return df

# def open_data(path_in,path_out,title,import_type='csv',export_type='pandas',key=[],saveA=True):
#
#     try:
#         if export_type == 'pandas':
#             data = pd.read_pickle(path_out+title+'.pkl')
#         elif export_type == 'numpy':
#             data = np.load(path_out+title+'.npy')
#     except:
#         if all(( import_type == 'csv', export_type == 'pandas')):
#             data = CSV(path_in,path_out,title,saveA=saveA,key=key)
#         elif all(( import_type == 'csv', export_type == 'numpy' )):
#             data = CSV(path_in,path_out,title,saveA=saveA)
#
#     return data



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
    x = str(x)
    x = x[2:]
    p = 1
    for i in range(len(x)):
        if x[i] == '0':
            p += 1
        else:
            break
    return p

def int_prec(x,p,err=False):

    if err == False:
        mult = 10**(p)
        x_new = int(round(x/mult)) * mult
    else:
        mult = 10**(p-1)
        x_new = int(round(x/mult)) * mult

    return x_new

class var:
    def __init__(self,val,err):
        self.val = val
        self.err = err
        # self.units = units

        try: # try for scalar values in val and err

            if self.err < 1:

                p = float_prec(self.err)
                self.p = p
                self.pval = round(self.val,p)
                self.perr = round(self.err,p)

            else:

                p = len(str(int(self.err)))
                self.p = p
                self.pval = int_prec(self.val,p)
                self.perr = int_prec(self.err,p,err=True)

        except:

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

def printvar(var):
    print(var.name,var.av,var.err,var.units)

def printerr(var):
    print(var.name,var.err,var.err/var.val)
