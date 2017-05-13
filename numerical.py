import numpy as np

#===============================================================================
""" numerical derivatives
of form f'(f,x) """
#-------------------------------------------------------------------------------

def D_fd(f, x, h=1e-9):
    """Forward difference
    calculates numerical derivative of f(x)

	args
	----
	f: function to take derivative of
	x: argument of f
	"""
    return (f(x + h) - f(x))/h

def D_cd(f, x, h=1e-8):
    """Central difference
    calculates numerical derivative of f(x)

	args
	----
	f: function to take derivative of
	x: argument of f
	"""
    return ( f(x+(h/2)) - f(x-(h/2)) )/h

def D_ep(x, f, h=1e-3):
    """Extended difference
	calculates numerical derivative of f(x)

	args
	----
	f: function to take derivative of
	x: argument of f
    """
    return ( 8*(f(x+(h/4)) - f(x-(h/4))) - (f(x+(h/2)) - f(x-(h/2))) ) / (3*h)

#===============================================================================
""" numerical integrators
1) f - function(t,x) that returns vector (x',x'')
2) t - time t
3) x - vector includes (x,x')
4) h - increment """
#-------------------------------------------------------------------------------

def euler(f,t,x,h):
    """ euler integrator
    3D ok

    arguments
    ---------
    x:  vector includes (x,x')
    f:  function that returns vector (x',x'')
    t:  time t
    h:  increment
    """

    return x + h * f(t, x)

def rk2(f,t,x,h):
    """ Runge-Kutta rk2 integrator
    3D ok

    arguments
    ---------
    x:  vector includes (x,x')
    f:  function that returns vector (x',x'')
    t:  time t
    h:  increment
    """

    k1  =   f(t,            x)
    k2  =   f(t + (1/2)*h,  x + (1/2)*h*k1)
    return x + h*k2

def rk4(f,t,x,h):
    """ Runge-Kutta rk4 integrator
    3D ok

    arguments
    ---------
    x:  vector includes (x,x')
    f:  function that returns vector (x',x'')
    t:  time t
    h:  increment
    """

    k1  =   f(t,            x)
    k2  =   f(t + (1/2)*h,  x + (1/2)*h*k1)
    k3  =   f(t + (1/2)*h,  x + (1/2)*h*k2)
    k4  =   f(t + h,        x + h*k3)
    return x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

def velocity_verlet(f,t,x,h):
    """ velocity verlet integrator
    f0 must be updated each step
    3D ok

    arguments
    ---------
    x:  vector includes (x,x')
    f:  function that returns vector (x',x'')
    t:  time t
    h:  increment
    """

    f1      =   f(t,    x)
    x[1]    +=  (1/2)*h*f1[1]
    x[0]    +=  h*x[1]
    f2      +=  f(t+h,  x)
    x[1]    +=  (1/2)*h*f2[1]
    return y

#===============================================================================
""" integrator engine
an example of how to use the integrators
with a tmax termination condition. """
#-------------------------------------------------------------------------------

def integrator(f,x0=(0,),v0=(0,),h=1e-4,tmax=100,method=rk4):

    Nsteps  =   tmax/h
    dim     =   len(x0)
    assert len(x0) == len(v0), 'x0 and v0 must be same length'
    
    T       =   h * np.arange(Nsteps)
    X       =   np.zeros(( len(T), 2, dim ))
    POS     =   np.zeros(( len(T), dim+1 ))

    # set initial conditions
    X[0,:,:]    =   x0,v0
    POS[0,:]    =   0,*x0

    for i,t in enumerate(T[:-1]):
        X[i+1,:,:]  =   method(f,t,X[i],h)
        POS[i+1,:]  =   t,*X[i+1,0,:]

    return POS
