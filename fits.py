import numpy as np
from scipy import optimize
import djak.gen as dg
from importlib import reload
reload(dg)

__all__ = ['SS_xx','SS_xy','m_exp','b_exp','SS_yy','SS_E','SS_R',
           'sig_y','sig_m','sig_b','RMSE','R2','lin_fit','least_square']

def SS_xx(X): #
	x_bar = np.average(X)
	sum = 0
	for x in X:
		sum += (x-x_bar)**2
	return sum

def SS_xy(X,Y):
	if len(X) != len(Y):
		print("X and Y must be same length")
		print(len(X),len(Y))
		raise ValueError

	x_bar = np.average(X)
	y_bar = np.average(Y)

	sum = 0
	N = len(X)
	for n in range(N):
		sum += (X[n]-x_bar)*(Y[n]-y_bar)

	return sum

def m_exp(X,Y):
    """slope of linear fit to data
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        print(len(X),len(Y))
        raise ValueError
    return SS_xy(X,Y)/SS_xx(X)

def b_exp(X,Y):
    """y-intercept of linear fit to data
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        print(len(X),len(Y))
        raise ValueError

    y_bar = np.average(Y)
    x_bar = np.average(X)
    m = m_exp(X,Y)

    return y_bar - m*x_bar

def SS_yy(Y):
    """total sum of squares
    Parameters
    ----------
    Y: range data
    """
    y_bar = np.average(Y)
    sum = 0
    for y in Y:
        sum += (y-y_bar)**2
    return sum

def SS_E(X,Y):
    """ error in sum of squares
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        raise ValueError

    N = len(X)
    m = m_exp(X,Y)
    b = b_exp(X,Y)

    sum = 0
    for n in range(N):
        y_exp = m*X[n] + b
        sum += (Y[n]-y_exp)**2

    return sum

def SS_R(X,Y):
    """ regression sum of squares
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        raise ValueError

    return SS_yy(Y) - SS_E(X,Y)

def sig_y(X,Y):
    """standard deviation of y(x)"""
    if len(X) != len(Y):
        print("X and Y must be same length")
        raise ValueError

    deg_free = len(X)-2 # degrees of freedom, length of array minus the 2 regression params
    assert deg_free != 0, "arrays with length 2 have 0 degrees of freedom. dg.fits.sig_y(X,Y)"
    return np.sqrt(SS_E(X,Y)/deg_free)

def sig_m(X,Y):
    """ standard deviation of fitted slope
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        raise ValueError

    return np.sqrt(sig_y(X,Y)/SS_xx(X))

def sig_b(X,Y):
    """ standard deviation of fitted y-intercept
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        print(len(X),len(Y))
        raise ValueError

    n = len(X)
    x_bar2 = np.average(X)**2

    return np.sqrt(sig_y(X,Y) * ( (1/n) + (x_bar2/SS_xx(X)) ) )

def RMSE(X,Y):
    """ Root Mean Square Error
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    if len(X) != len(Y):
        print("X and Y must be same length")
        raise ValueError

    n = len(X)
    return np.sqrt(SS_E(X,Y)/n)

def R2(X,Y):
    """ Coefficient of Determination - the variability of y_i accounted for by linear model
    Parameters
    ----------
    X: domain data
    Y: range data
    """
    return SS_R(X,Y)/SS_yy(Y)

def lin_fit(X_data,Y_data,X_lin):
    """ complete linear fit of data

    Parameters
    ----------
    X_data: domain of data
    Y_data: range of data
    X_lin:  higher resolution domain to plot smooth linear fit

    Returns
    -------
    Y_lin:  high resolution range of fitted line
    m:      dg.var class object - slope
    b:      dg.var class object - y-intercept
    """

    X_data=np.array(X_data)
    Y_data=np.array(Y_data)
    assert len(X_data) >= 2, "there have to be at least 2 data points to find a fit. df.lin_fit( X_data, Y_data, X_fit )"

    if len(X_data) > 2:
        m = m_exp(X_data,Y_data)
        b = b_exp(X_data,Y_data)
        Y_lin = m*X_lin + b

        m_err = sig_m(X_data,Y_data)
        b_err = sig_b(X_data,Y_data)

        m_var = dg.var(m,m_err)
        b_var = dg.var(b,b_err)

    else:
        m = ( Y_data[1] - Y_data[0] ) / ( X_data[1] - X_data[0] )
        def fit_func(b,x):
            return m*x + b
        def errfunc(b,x,y):
            return y - fit_func(b,x)
        qout,success = optimize.leastsq(errfunc,[1],args=(X_data,Y_data),maxfev=5000)
        b = qout
        Y_lin = m * X_lin + b

        m_var = dg.var(m,0)
        b_var = dg.var(b,0)

    R_2 = R2(X_data,Y_data)

    return {'Y_fit':Y_lin, 'm':m_var, 'b':b_var, 'R2':R_2}

def least_square(fitfunc,par,X_data,Y_data,cycles=5000):
    """ complete least square fit to data

    Parameters
    ----------
    fitfunc: name of fitting function
    par:     list/array of fitting parameters
    X_data:  domain of data
    Y_data:  range of data
    cycles:  ** max iterations of fitting - default=5000

    Returns
    -------
    qout: fitted parameters
    """

    def errfunc(p,x,y):
        return y - fitfunc(p,x)
    xdata = np.array(X_data)
    ydata = np.array(Y_data)

    qout,success = optimize.leastsq(errfunc,par,args=(xdata,ydata),maxfev=cycles)

    if success == False:
        print("OMFG - Failed!!")
    return qout

# http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i

# def fit_leastsq(p0, datax, datay, function):

#     #errfunc = lambda p, x, y: function(x,p) - y
#     def errfunc(p,x,y):
#         return y - function(p,x)

#     pfit, pcov, infodict, errmsg, success = \
#         optimize.leastsq(errfunc, p0, args=(datax, datay), \
#                           full_output=1, epsfcn=0.000001)

#     if (len(datay) > len(p0)) and pcov is not None:
#         s_sq = (errfunc(pfit, datax, datay)**2).sum()/(len(datay)-len(p0))
#         pcov = pcov * s_sq
#     else:
#         pcov = np.inf

#     error = []
#     for i in range(len(pfit)):
#         try:
#           error.append(np.absolute(pcov[i][i])**0.5)
#         except:
#           error.append( 0.00 )
#     pfit_leastsq = pfit
#     perr_leastsq = np.array(error)
#     return pfit_leastsq, perr_leastsq

# def fit_curvefit(p0, datax, datay, function, yerr=err_stdev, **kwargs):

#     pfit, pcov = \
#          optimize.curve_fit(f,datax,datay,p0=p0,\
#                             sigma=yerr, epsfcn=0.0001, **kwargs)
#     error = []
#     for i in range(len(pfit)):
#         try:
#           error.append(np.absolute(pcov[i][i])**0.5)
#         except:
#           error.append( 0.00 )
#     pfit_curvefit = pfit
#     perr_curvefit = np.array(error)
#     return pfit_curvefit, perr_curvefit
