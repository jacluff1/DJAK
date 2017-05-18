import numpy as np
from djak.constants import SI

#===============================================================================
""" Jean Cloud """
#-------------------------------------------------------------------------------

def Jean_mass_density(rho0,mu,T,units=SI):
    """ Jeans mass in terms of initial density

    parameters
    ----------
    rho0:   initial density
    mu:     average particle mass
    T:      temperature of cloud (K)
    units:  ** dictionary of units, default = SI

    returns
    -------
    scalar
    """

    return ( (5*units['k'].value*T) / (units['G'].value*mu) )**(3/2) * ( 3 / (4*np.pi*rho0) )**(1/2)

def Jean_radius_density(rho0,mu,T,units=SI):
    """ Jeans radius in terms of density

    parameters
    ----------
    rho0:   initial density
    mu:     average particle mass
    T:      temperature of cloud (K)
    units:  ** dictionary of units, default = SI

    returns
    -------
    scalar
    """

    return ( (15*units['k'].value*T) / (4*np.pi*units['G'].value*mu*rho0) )**(1/2)

def Jean_mass_R(R,mu,T,units=SI):
    """ Jeans mass in terms of initial cloud radius

    parameters
    ----------
    R:      initial cloud radius
    mu:     average particle mass
    T:      temperature of cloud (K)
    units:  ** dictionary of units, default = SI

    returns
    -------
    scalar
    """

    return (5*units['k'].value*T*R) / (mu*units['G'].value)

def Jean_radius_M(M,mu,T,units=SI):
    """ Jeans radius in terms of cloud mass

    parameters
    ----------
    M:      cloud mass
    mu:     average particle mass
    T:      temperature of cloud (K)
    units:  ** dictionary of units, default = SI

    returns
    -------
    scalar
    """

    return (M*mu*units['G'].value) / (5*units['k'].value*T)
