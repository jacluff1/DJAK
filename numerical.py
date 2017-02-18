import numpy as np 

# DERIVATIVES

def D_fd(y, t, h=1e-9):
    """Forward difference
    calculates numerical derivative

	args
	----
	y: function to take derivative of
	t: argument of y 
	"""
    return (y(t + h) - y(t))/h

def D_cd(y, t, h=1e-8):
    """Central difference
    calculates numerical derivative

	args
	----
	y: function to take derivative of
	t: argument of y 
	"""
    return ( y(t+(h/2)) - y(t-(h/2)) )/h

def D_ep(y, t, h=1e-3):
    """Extended difference
	calculates numerical derivative

	args
	----
	y: function to take derivative of
	t: argument of y 
    """
    return ( 8*(y(t+(h/4)) - y(t-(h/4))) - (y(t+(h/2)) - y(t-(h/2))) ) / (3*h)

 # NUMERICAL INTEGRATORS

 def FOeuler(FOLDE,x0,t0,tf,dt=1e-5):
	"""First Order Euler Engine

	arguments
	---------
	FOLDE: defined function name of first order equation
	x0: initial pos - scalar
	t0: initial time - scalar
	tf: final t - scalar
	**dt: time step value - scalar

	returns
	-------
	T: time values - numpy array
	X: position values - numpy array
	"""

	# Create numpy array of time values
	T = np.arange(t0,tf+dt,dt)

	# Create empty numpy array of postion values
	X = np.zeros(len(T))

	# Add initial position to position array
	X[0] = x0

	# first order engine
	for i in range(len(T[:-1])):
		# filling in Y with euler engine
		X[i+1] = X[i] + FOLDE(T[i])*dt

	return T,X

def SOeuler(SOLDE,x0,v0,t0,tf,dt=1e-5):
	"""Second Order Euler Engine

	arguements
	----------
	SOLDE: defined function for second order equation
	x0: inital position - scalar
	v0: initial velocity - scalar
	t0: initial time - scalar
	tf: final time - scalar
	**dt: time step - scalar

	returns
	-------
	T: time values - numpy array
	X: position values - numpy array
	"""

	# define X, Y, and Y' arrays
	T = np.arange(t0,tf+dt,dt)
	X = np.zeros(len(T))
	V = np.zeros(len(T))

	# initalize Y and Y' arrays
	X[0],V[0] = x0,v0

	#SO engine
	for i in range(len(T[:-1])):
		X[i+1] = X[i] + V[i]*dt
		V[i+1] = V[i] + SOLDE(T[i],X[i],V[i])*dt

	return T,X