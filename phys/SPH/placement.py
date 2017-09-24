import numpy as np
import djak.coordinates as dc

O   =   np.zeros(3)

def sphere_random(N,radius,R=O):
    """ constructs particles into a randimized sphere

    args
    ----
    N:      number of particles
    radius: radius of sphere
    R:      ** COM position vector - default = Origin
    """
    r   =   radius * np.random.rand(N)
    t   =   np.pi * np.random.rand(N)
    p   =   2 * np.pi * np.random.rand(N)

    x,y,z   =   dc.SPC2CC(r,t,p)

    return x + R[0], y + R[1], z + R[2]

def box_random(N,length,R=O):
    """ constructs particles into randomized box

    args
    ----
    N:      number of particles
    length: length of box
    R:      ** COM position vector - default = Origin
    """

    # make sure N is divisible by 8
    N   +=  divmod(N,8)[1]

    def octant():
        return (length/2) * np.random.rand(int(N/8),3)

    # form octants
    O1  =   octant()
    O2  =   octant()
    O3  =   octant()
    O4  =   octant()
    O5  =   octant()
    O6  =   octant()
    O7  =   octant()
    O8  =   octant()

    # octants 2,3,6, and 8 have -x
    O2[:,0] *= -1
    O3[:,0] *= -1
    O6[:,0] *= -1
    O8[:,0] *= -1

    # octants 3,4,7, and 8 have -y
    O3[:,1] *= -1
    O4[:,1] *= -1
    O7[:,1] *= -1
    O8[:,1] *= -1

    # octants 5,6,7, and 8 have -z
    O5[:,2] *= -1
    O6[:,2] *= -1
    O7[:,2] *= -1
    O8[:,2] *= -1

    M   =   np.vstack((O1,O2,O3,O4,O5,O6,O6,O8))
    x   =   M[:,0]
    y   =   M[:,1]
    z   =   M[:,2]

    return x + R[0], y + R[1], z + R[2]

def sphere_uniformish(N,radius,R=O):
    """ constructs random box, then reallocates particles outside sphere

    args
    ----
    N:      number of particles
    radius: cloud radius
    R:      ** COM position vector
    """

    # make sure N is divisible by 8
    N       +=  divmod(N,8)[1]

    # make box of particles
    x,y,z   =   box_random(N,radius*2)

    # if particles are outside sphere, relocate them in random positions within sphere
    for i in range(N):
        if np.sqrt( x[i]**2 + y[i]**2 + z[i]**2 ) > radius:
            r       =   radius * np.random.rand(1)[0]
            t       =   np.pi * np.random.rand(1)[0]
            p       =   2 * np.pi * np.random.rand(1)[0]

            x[i]    =   r * np.sin(t) * np.cos(p)
            y[i]    =   r * np.sin(t) * np.sin(p)
            z[i]    =   r * np.cos(t)

    return x + R[0], y + R[1], z + R[2]


    return
