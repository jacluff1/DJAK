import numpy as np

#===============================================================================
""" Sherical Polar Coordinates (SPC) """
#-------------------------------------------------------------------------------

def SPC2CC(radius,theta,phi):
    """ spherical polar coordinates to cartesian coordinates

    args
    ----
    radius: obvious
    theta:  colatidudnal angle
    phi:    azimuthal angle
    """
    x   =   radius * np.sin(theta) * np.cos(phi)
    y   =   radius * np.sin(theta) * np.sin(phi)
    z   =   radius * np.cos(theta)

    return x,y,z
