import sympy as sy
import sympy.physics.vector as spv
import numpy as np

#===============================================================================
""" dictionary results
var:    tuple of symbolic variables
par:    tuple of numerical variable values
y:      symbolic expression
f:      lambda function of symbolic expression
result: f(*par)
"""
#-------------------------------------------------------------------------------

def dic_result(var,par,y,modules=None):
    """ construct dictionary of results

    arguments
    ---------
    var:    tuple of variables
    par:    tuple of numerical 'var' values
    y:      symbolic expression

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    if modules == None:
        f   =   sy.lambdify( [*var] , y)
    else:
        f   =   sy.lambdify( [*var] , y , modules=modules )
    result  =   f(*par)

    return {'var':var, 'par':par, 'y':y, 'f':f, 'result':result}

def test_result(dic):
    """ print a quick test of returned function dictionary

    arguments
    ---------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   dic['var']
    par =   dic['par']
    y   =   dic['y']
    f   =   dic['f']

    print(var)
    print(par)
    print(y)
    print(f(*par))

#===============================================================================
""" symbolic vector functions """
#-------------------------------------------------------------------------------

N           =   spv.ReferenceFrame('N')
u_z         =   np.zeros(3)
u           =   sy.var('u')
q1,q2,q3    =   sy.symbols('q1:4')

def makevector(RF=N):
    """ make a generic vector

    arguments
    ---------
    RF: default = N

    returns
    -------
    sympy.physics.vector.vector.Vector
    """

    return q1*RF.x + q2*RF.y + q3*RF.z

def np2vector(u=u_z,RF=N):
    """ construct a sympy vector from a 1D numpy array ( len = 3 )

    arguments
    ---------
    u:  default = array[ 0 , 0 , 0 ], numpy array to convert
    RF: default = N, standard CC at origin

    returns
    -------
    sympy.physics.vector.vector.Vector
    """

    e1,e2,e3    =   u[0],u[1],u[2]
    vec         =   e1*RF.x + e2*RF.y + e3*RF.z
    return vec

def tup2vector(u1,u2,u3,RF=N):
    """ make a symbolic vector

    arguments
    ---------
    u:  symbol

    returns
    -------
    sympy.physics.vector.vector.Vector
    """

    vec =   u1*RF.x + u2*RF.y + u3*RF.z
    return vec

#===============================================================================
""" Matrix functions """
#-------------------------------------------------------------------------------

A   =   np.ones((3,3))
sym =   'M'
vec =   sy.var('%s_1:4' % sym)
row =   sy.Matrix(1,3,vec)

def np2Matrix(A=A,sym=sym,dim=2):
    """ get a 2D numpy matrix and return both a symbolic and numerical
    sympy Matrix of the same length

    arguments
    ---------
    A:      default = np.ones((3,3))
    sym:    default = 'M'

    returns
    -------
    tuple of 2 sympy.matrices.dense.MutableDenseMatrix
    """

    if dim == 2:
        row =   A.shape[0]
        col =   A.shape[1]

        Mp  =   sy.Matrix(A)

        Mv  =   Mp[:,:]
        for r in range(row):
            for c in range(col):
                Mv[r,c] =   sy.var('%s_%s_%s' % (sym,r,c) )

    elif dim == 1:
        row =   A.shape[0]

        Mp  =   sy.Matrix(A)
        Mv  =   Mp[:,:]
        for r in range(row):
            Mv[r]   =   sy.var('%s_%s' % (sym,r) )

    return Mv,Mp

def symb_norm(Mrow=row):
    tot     =   0
    length  =   len(Mrow)
    for i in range(length):
        tot += Mrow[i]**2
    return sy.sqrt( tot )
