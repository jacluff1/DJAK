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