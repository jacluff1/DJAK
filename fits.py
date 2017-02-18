import numpy as np
from scipy import optimize

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

def m_exp(X,Y): # expected slope
	if len(X) != len(Y):
		print("X and Y must be same length")
		print(len(X),len(Y))
		raise ValueError

	return SS_xy(X,Y)/SS_xx(X)

def b_exp(X,Y): # expected y-intercept
	if len(X) != len(Y):
		print("X and Y must be same length")
		print(len(X),len(Y))
		raise ValueError

	y_bar = np.average(Y)
	x_bar = np.average(X)
	m = m_exp(X,Y)

	return y_bar - m*x_bar

def SS_yy(Y): # total sum of squares, Y is np array
	y_bar = np.average(Y)
	sum = 0
	for y in Y:
		sum += (y-y_bar)**2
	return sum

def SS_E(X,Y): # error sum of squares, Y is np array, m and b are linear parameters
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

def SS_R(X,Y): # Regression sum of squares
	if len(X) != len(Y):
		print("X and Y must be same length")
		raise ValueError

	return SS_yy(Y) - SS_E(X,Y)

def sig_y(X,Y): # standard deviation of y(x)
	if len(X) != len(Y):
		print("X and Y must be same length")
		raise ValueError

	deg_free = len(X)-2 # degrees of freedom, length of array minus the 2 regression params
	return np.sqrt(SS_E(X,Y)/deg_free)

def sig_m(X,Y): # standard deviation of m
	if len(X) != len(Y):
		print("X and Y must be same length")
		raise ValueError

	return np.sqrt(sig_y(X,Y)/SS_xx(X))

def sig_b(X,Y): # variance of y-intercept
	if len(X) != len(Y):
		print("X and Y must be same length")
		print(len(X),len(Y))
		raise ValueError

	n = len(X)
	x_bar2 = np.average(X)**2

	return np.sqrt(sig_y(X,Y) * ( (1/n) + (x_bar2/SS_xx(X)) ) )

def RMSE(X,Y): # room mean square error
	if len(X) != len(Y):
		print("X and Y must be same length")
		raise ValueError

	n = len(X)
	return np.sqrt(SS_E(X,Y)/n)

def R2(X,Y): # Coefficient of Determination - the variability of y_i accounted for by linear model
	return SS_R(X,Y)/SS_yy(Y)

def lin_fit(X_data,Y_data,X_lin):
    X_data=np.array(X_data)
    Y_data=np.array(Y_data)
    m = m_exp(X_data,Y_data)
    b = b_exp(X_data,Y_data)
    Y_lin = m*X_lin + b
    
    #m_err = np.sqrt(s_m2(X_data,Y_data))
    #b_err = np.sqrt(s_b2(X_data,Y_data))
    #Yp = (m+m_err)*X_lin + (b+b_err)
    #Ym = (m-m_err)*X_lin + (b-b_err)
    return Y_lin

def least_square(fitfunc,par,X_data,Y_data,cycles=5000):
    
    def errfunc(p,x,y):
        return y - fitfunc(p,x)
    xdata = np.array(X_data)
    ydata = np.array(Y_data)
    
    qout,success = optimize.leastsq(errfunc,par,args=(xdata,ydata),maxfev=cycles)
    
    if success == False:
        print("OMFG - Failed!!")
    return qout



