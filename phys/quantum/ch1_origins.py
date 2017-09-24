import scipy.integrate as si
import sympy as sy
import sympy.physics.vector as spv
import numpy as np
# import cmath as cm
# from inspect import isfunction
import djak.constants as dc
import djak.symbolic as ds
from djak.symbolic import dic_result
import pdb

#===============================================================================
""" useful shortcuts """
#-------------------------------------------------------------------------------

SI  =   dc.SI
eV  =   dc.eV

#===============================================================================
""" 1.2 Particle Aspect of Radiation """
#-------------------------------------------------------------------------------

#===============================================================================
""" 1.2.1 Blackbody Radiation """
#-------------------------------------------------------------------------------

def intensity_blackbody(efficiency=1,temperature=1, units=SI):
    """ power/area

    parameters
    ----------
    energy:         default = 1 [ energy ]
    temperature:    default = 1 [ K ]
    units:          default SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('e sigma t')
    par =   efficiency, units['sigma'], temperature

    y   =   e * sigma * t**4

    return dic_result(var,par,y)

def distribution_wein(a0=1,a1=1,frequency=1,temperature=1):
    """ Wien's energy density distribution

    arguments
    ---------
    a1:             default = 1, fitting parameter 1
    a2:             default = 1, fitting parameter 2
    frequency:      default = 1 [ 1 / time ]
    temperature:    default = 1 [K]

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('a:2 nu t')
    par =   a0, a1, frequency, temperature

    y   =   var[0] * nu**3 * sy.exp( - var[1]*nu/t)

    return dic_result(var,par,y)

def distribution_rayleigh_jean(frequency=1,temperature=1, units=SI):
    """ Rayleigh Jean's energy density distribution

    arguments
    ---------
    frequency:      default = 1 [ 1 / time ]
    temperature:    default = 1 [ K ]
    units:          default SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('pi nu c k t')
    par =   np.pi, frequency, units['c'], units['k'], temperature

    y   =   ( 8 * pi * nu**2 * k * t ) / ( c**3 )

    return dic_result(var,par,y)

def distribution_planck_nu(frequency=1,temperature=1, units=SI,printA=False):
    """ Planck energy density distribution WRT nu

    arguments
    ---------
    frequency:      default = 1 [ 1 / time ]
    temperature:    default = 1 [ K ]
    units:          default SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('pi nu h c k t')
    par =   np.pi, frequency, units['h'], units['c'], units['k'], temperature

    y   =   ( 8 * pi * nu**2 * h * nu ) / ( c**3 * (sy.exp(h*nu/k/t) - 1) )

    return dic_result(var,par,y)

def distribution_planck_lambda(wavelength=1,temperature=1, units=SI,printA=False):
    """ Planck energy density distribution WRT lambda

    arguments
    ---------
    wavelength:     default = 1 [ length ]
    temperature:    default = 1 [ K ]
    units:          default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('pi h c l k t')
    par =   np.pi, units['h'], units['c'], wavelength, units['k'], temperature

    y   =   ( 8 * pi * h * c ) / l**5 / ( sy.exp(h*c/l/k/t) - 1 )

    return dic_result(var,par,y)

def displacement_wein(temperature=1, units=SI):
    """ wavelength of blackbody peak energy density

    arguments
    ---------
    temperature:    default = 1 [K]
    units:          default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('h c k t')
    par = units['h'], units['c'], units['k'], temperature

    y   =   h * c / 4.9663 / k / t

    return dic_result(var,par,y)

#===============================================================================
""" 1.2.2 Photoelectric Effect """
#-------------------------------------------------------------------------------

def frequency_threshold(workFunction=1, units=eV):
    """ threshold frequency of material

    arguments
    ---------
    workFunction:   default = 1 [ energy ]
    units:          default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('W h')
    par =   workFunction, units['h']

    y   =   W / h

    return dic_result(var,par,y)

def KE_ejected_electron(frequency=1,workFunction=1, units=eV):
    """ kinetic energy of ejected photon

    arguments
    ---------
    frequency:      default = 1 [ 1 / time ]
    workFunction:   default = 1 [ energy ]
    units:          default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('h nu W')
    par =   units['h'], frequency, workFunction

    y   =   h*nu - W

    return dic_result(var,par,y)

def potential_stopping(frequency=1,workFunction=1, units=eV):
    """ stopping potential WRT frequency

    arguments
    ---------
    frequency:      default = 1 [ 1 / time ]
    workFunction:   default = 1 [ energy ]
    units:          default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('h e nu W')
    par =   units['h'], units['e'], frequency, workFunction

    y   =   (h/e) * nu - (W/e)

    return dic_result(var,par,y)

#===============================================================================
""" 1.2.3 Compton Effect """
#-------------------------------------------------------------------------------

def wavelength_compton_shift(angle=0, units=SI):
    """ shift in wavelength from compton scatter

    arguments
    ---------
    angle:  default = 0 [degree], scattering angle
    units:  default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('h m_e c theta')
    par =   units['h'], units['m_e'], units['c'], np.deg2rad(angle)

    y   =   h / (m_e * c) * ( 1 - sy.cos(theta) )

    return dic_result(var,par,y)

def energy_scattered_electron(photonEnergy=1, units=eV):
    """ energy of scattered electron from compton scatter

    arguments
    ---------
    photonEnergy:   default = 1 [ energy ], energy of photon before interaction
    units:          default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('m_e c E')
    par =   units['m_e'], units['c'], photonEnergy

    m0  =   m_e * c**2
    y   =   m0 / ( (m0/E) + 2 )

    return dic_result(var,par,y)

#===============================================================================
""" 1.2.4 Pair Production """
#-------------------------------------------------------------------------------

def energy_photon_pair_production(k_electron=1,k_proton=1,k_nucleus=1,m_nucleus=1, units=eV):
    """

    arguments
    ---------
    k_electron: default = 1 [ energy ], kinetic energy of electron
    k_proton:   default = 1 [ energy ], kinetic energy of positron
    k_nucleus:  default = 1 [ energy ], kinetic energy of nucleus
    m_nucleus:  default = 1 [ mass ], mass of nucleus
    units:      default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('m_e c k_e k_p k_n m_n ')
    par =   units['m_e'], units['c'], k_electron, k_proton, k_nucleus, m_nucleus

    y   =   (2 * m_e + m_n) * c**2 + k_e + k_p + k_n

    return dic_result(var,par,y)

#===============================================================================
""" 1.3 Wave Aspect of Particles """
#-------------------------------------------------------------------------------

#===============================================================================
""" 1.3.1 de Broglie's Hypothesis: Matter Waves """
#-------------------------------------------------------------------------------

def wavelength_deBroglie(momentum=1, units=SI):
    """ deBroglie wavelength

    arguments
    ---------
    momentum:   default = 1 [ mass velocity ], magnitude of momentum
    units:      default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('h p')
    par =   units['h'], momentum

    y   =   h / p
    return dic_result(var,par,y)

def waveVector_deBroglie(momentum=1, units=SI):
    """ deBroglie wave vector

    arguments
    ---------
    momentum:   default = 1 [ mass velocity ], magnitude of momentum
    units:      default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('p hbar')
    par =   momentum, units['hbar']

    y   =   p / hbar
    return dic_result(var,par,y)

#===============================================================================
""" 1.4 Particles verses Waves """
#-------------------------------------------------------------------------------

# #===============================================================================
# """ 1.4.1 Classical View of Partiles and Waves """
# #-------------------------------------------------------------------------------
#
# k_z     =   np.zeros(3)
# r_z     =   np.array(k_z, copy=True)
# psi_z   =   np.zeros((2,3))
#
# def psi_3Dplane(amplitude=1,waveVector=k_z,position=r_z,afrequency=1,time=0,conj=False):
#     """ 3D plane wave
#
#     arguments
#     ---------
#     amplitude:  default = 1 [ length ]
#     waveVector: default = array[0,0,0] [ 1 / length ]
#     position:   default = array[0,0,0] [ length ]
#     afrequency: default = 0 [ radians / time ]
#     time:       default = 0 [ time ]
#     conj:       default = False
#
#     returns
#     -------
#     dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
#     """
#
#     var1    =   sy.var('A omega t',real=True)
#     kvar    =   sy.var('k1:4', real=True)
#     rvar    =   sy.var('r1:4', real=True)
#     var     =   var1 + kvar + rvar
#
#     par1    =   amplitude, afrequency, time
#     kpar    =   tuple( [waveVector[n] for n in range(3)] )
#     rpar    =   tuple( [position[n] for n in range(3)] )
#     par     =   par1 + kpar + rpar
#
#     k       =   ds.tup2vector(*kvar)
#     r       =   ds.tup2vector(*rvar)
#     y       =   A * sy.exp( sy.I * (spv.dot(k,r) - omega*t) )
#
#     if conj == True:
#         y = sy.conjugate(y)
#
#     return dic_result(var,par,y, modules='numpy')
#
# def intensity_particle_source(psi=psi_z):
#     """ find the total intensity from particle sources
#
#     arguments
#     ---------
#     psi:    2D array of psi vectors ( wave , e )
#
#     returns
#     -------
#     dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
#     """
#
#     N           =   len(psi[:,0])
#
#     Mv,Mp       =   ds.np2Matrix(psi,sym='psi')
#     var         =   tuple(Mv)
#     par         =   tuple(Mp)
#
#     y   =   0
#     for n in range(N):
#         psi_n   =   Mv.row(n)
#         y       +=  ds.symb_norm(Mrow=psi_n)
#
#     return dic_result(var,par,y)
#
# def intensity_wave_source(psi=psi_z):
#     """ find the total intensity from wave sources
#
#     arguments
#     ---------
#     psi:    2D array of psi vectors ( wave , e )
#
#     returns
#     -------
#     dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
#     """
#     N       =   psi.shape[0]
#     L       =   psi.shape[1]
#
#     Mv,Mp   =   ds.np2Matrix(A=psi,sym='psi')
#     var     =   tuple(Mv)
#     par     =   tuple(Mp)
#
#     y1  = 0
#     for n in range(N):
#         y   +=  Mv.row(n)
#     y   =   sy.sqrt( tot )
#
#     return dic_result(var,par,y)
#
# #===============================================================================
# """ 1.4.4 Principle of Linear Superposition """
# #-------------------------------------------------------------------------------
#
# alpha   =   np.ones(5)
# psi     =   np.ones((5,3))
#
# def psi_linear(alpha=alpha,psi=psi,conj=False):
#     """ find wave function from linear combination of wave functions
#
#     arguments
#     ---------
#     alpha:  1D array of coefficients ( wave )
#     psi:    2D array of psi vectors ( wave , e )
#     conj:   ** default = False
#
#     returns
#     -------
#     vector
#     """
#
#     psiSum  =   0
#     R       =   psi.shape[0]
#     C       =   psi.shape[1]
#     assert alpha.shape[0] == alpha.shape[0], "check alpha and psi dimentions"
#
#     Mv,Mp   =   ds.np2Matrix(A=psi,sym='psi')
#     Av,Ap   =   ds.np2Matrix(A=alpha,sym='alpha',dim=1)
#
#     var =   tuple(Av) + tuple(Mv)
#     par =   tuple(Ap) + tuple(Mp)
#
#     y   =   sum( Av.T * Mv )
#
#     if conj == True:
#         y   =   y.conjugate()
#
#     return dic_result(var,par,y)

#===============================================================================
""" 1.5 Indeterministic Nature of the Microphysical World """
#-------------------------------------------------------------------------------

#===============================================================================
""" 1.5.1 Heisenberg's Uncertainty Principle """
#-------------------------------------------------------------------------------

deltap  =   'Delta_p'

def uncertainty_heisenberg(uncertainty=1,symb=deltap, units=SI):
    """ finds the uncertainty between two conjugate paired variables

    arguments
    ---------
    var:    default = dp, uncertainty of known conjugate variable
    units:  default = SI

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var1    =   sy.var(symb)
    var2    =   sy.var('hbar')
    var     =   var1,var2
    par     =   uncertainty,units['hbar']

    y   =   var2 / var1

    return dic_result(var,par,y)

#===============================================================================
""" 1.5.2 Probabilistic Interpretation """
#-------------------------------------------------------------------------------

def probability_density(dic):
    """ probability density psi* psi
    requires psi function to have conjugate keyword argument

    arguments
    ---------
    dic:    returned dictionary from desired psi function

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   dic['var']
    par =   dic['par']
    y1  =   dic['y']
    y   =   y1.conjugate() * y
    return dic_result(var,par,y)

#===============================================================================
""" 1.6 Atomic Transitions and Spectroscopy """
#-------------------------------------------------------------------------------

#===============================================================================
""" 1.6.2 Bohr Model of the Hydrogen Atom """
#-------------------------------------------------------------------------------

Z1  =   1
mn1 =   eV['m_n']

def radius_bohrlike(atomicZ=Z1,mass_n=mn1,num=1,units=eV):
    """ find radius of bohr-like atom

    arguments
    ---------
    atomicZ:    default = 1, atomic number
    mass_N:     default = eV['m_n'], mass of nucleus
    num:        default = 1, orbital number
    units:      default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('me mn a0 Z n')
    par =   units['m_e'], mass_n, units['a_0'], atomicZ, num

    y   =   (1 + me/mn) * (a0/Z) * n**2

    return dic_result(var,par,y)

def energy_bohr_orbital(atomicZ=Z1,mass_n=mn1,num=1,units=eV):
    """ find orbital energy of bohr-like atom

    arguments
    ---------
    atomicZ:    default = 1, atomic number
    mass_N:     default = eV['m_n'], mass of nucleus
    num:        default = 1, orbital number
    units:      default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('Z me mn R n')
    par =   atomicZ, units['m_e'], mass_n, units['Rydberg'], num

    y   =   - Z**2 / (1 + (me/mn)) * (R/n**2)
    return dic_result(var,par,y)

def energy_bohr_photon(atomicZ=Z1,mass_n=mn1,num1=1,num2=2,units=eV):
    """ energy of absorbed/emitted photon from bohr-like atom

    arguments
    ---------
    atomicZ:    default = 1, atomic number
    mass_N:     default = eV['m_n'], mass of nucleus
    num1:       default = 1, initial orbital number
    num2:       default = 2, final orbital number
    units:      default = eV

    returns
    -------
    dictionary: keys('var', 'par', 'y' , 'f' , 'result' )
    """

    var =   sy.var('R Z me mn m n')
    par =   units['Rydberg'], atomicZ, units['m_e'], mass_n, num2, num1

    y   =   R * Z**2 / ( (me/mn) + 1) * ( 1/n**2 - 1/m**2 )
    return dic_result(var,par,y)

#===============================================================================
""" 1.8 Wave Packets """
#-------------------------------------------------------------------------------

#===============================================================================
""" 1.8.1 Localized Wave Packets """

x1          =   0
omega1      =   0
t1          =   0

#-------------------------------------------------------------------------------

def example_phik():

    var =   sy.var('a pi k k0')
    y   =   (a**2 / (2*pi))**sy.Rational(1,4) * sy.exp( -a**2 * (k-k0)/4 )
    return {'var':var, 'y':y}

def psi_packet_k(pos=x1,angular_f=omega1,time=t1,phi_dic=example_phik()):
    """ find psi(x,t) from phi(k)

    arguments
    ---------
    pos:        default = 0, 1D position
    angular_f:  default = 0, angular frequency
    time:       default = 0, time
    phi_dic:    dictionary result of phi function(k)

    returns
    -------
    scalar
    """

    phi_var =   phi_dic['var']
    phi_par =   1, np.pi, 
    phi_y   =   phi_dic['y']

    psi_var =   sy.var('x omega t i')
    psi_par =   pos, angular_f, time, sy.I

    # var     =   phi_var + psi_var
    # par     =   phi_par + psi_par
    var     =   psi_var
    par     =   psi_var

    k       =   sy.symbols('k')
    y1      =   phi_y * sy.exp( i * (k*x - omega*t) )

    return y1

    # def integrand(k,x,omega,t):
    #     return phi(k) * np.exp( 1j * (k*x - omega*t) )
    #
    # I   =   si.quad(integrand, -np.inf, np.inf, args=(x,omega,t) )
    #
    # return 1/np.sqrt( 2 * np.pi ) * ( I[0] - I[1] )
#
# def phi_packet_k(k,psi):
#     """ find phi(k) from psi(x,0)
#
#     arguments
#     ---------
#     k:      scalar  -   1D wave number
#     psi:    function(x,0)
#
#     returns
#     -------
#     scalar
#     """
#
#     def integrand(x,k):
#         return psi(x,0) * np.exp( -1j * k * x )
#
#     I   =   si.quad(integrand, -np.inf, np.inf, args=k )
#
#     return 1/np.sqrt( 2 * np.pi ) * ( I[0] - I[1] )
#
# def psi_packet_p(x,E,t,phi, units=SI):
#     """ find psi(x,t) from phi(p)
#
#     arguments
#     ---------
#     x:      scalar  -   1D position
#     E:      scalar  -   energy
#     t:      scalar  -   time
#     phi:    function(p)
#     units:  ** default = SI
#
#     returns
#     -------
#     scalar
#     """
#
#     def integrand(p,x,E,t):
#         return phi(p) * np.exp( 1j * (p*x - E*t)/units['hbar'] )
#
#     I   =   si.quad(integrand, -np.inf, np.inf, args=(x,E,t) )
#
#     return 1/np.sqrt( 2 * np.pi * units['hbar'] ) * ( I[0] - I[1] )
#
# def phi_packet_p(p,psi, units=SI):
#     """ find phi(p) from psi(x,0)
#
#     arguments
#     ---------
#     p:      scalar  -   1D momentum
#     psi:    function(x,0)
#     units:  ** default = SI
#
#     returns
#     -------
#     scalar
#     """
#
#     def integrand(x,p):
#         return psi(x,0) * np.exp( -1j * p * x / units['hbar'] )
#
#     I   =   si.quad(integrand, -np.inf, np.inf, args=p)
#
#     return 1/np.sqrt( 2 * np.pi * units['hbar'] ) * ( I[0] - I[1] )
#
# #===============================================================================
# """ 1.8.3 Motion of Wave Packets """
# #-------------------------------------------------------------------------------
#
# def find_phase_velocity(omega,k):
#     """ find the phase speed from omega(k)
#
#     arguments
#     ---------
#     omega:  sympy symbolic function
#     k:      scalar  -   wavenumber
#
#     returns
#     -------
#     scalar
#     """
#
#     try:
#         sy.diff(omega,k)

# def test_functions():
