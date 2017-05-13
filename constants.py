from numpy import pi

#===============================================================================
""" unit class """
#-------------------------------------------------------------------------------

class conv:
	def __init__(self,value,units):
		self.value	=	value
		self.units 	= 	units

#===============================================================================
""" conversion factors
- multiply by conversion factor to convert from SI
- divide by conversion factor to convert to SI
- time and astronomical distances based on Julian year """
#-------------------------------------------------------------------------------

# length	1 / m
m		=	conv(1, 'm')
cm 		=	conv(100, 'cm')
ly 		=	conv(1 / 9.460730472e15, 'ly')
pc		=	conv(1 / 3.0856776e16, 'pc')
kpc		= 	conv(1 / (pc.value*1e3), 'kpc')
Mpc 	=	conv(1 / (pc.value*1e6), 'Mpc')
AU		=	conv(1 / 1.4959787066e11, 'AU')
solarR	=	conv(1 / 6.95508e8, 'R$_\odot$')
earthR 	=	conv(1 / 6.378136e6, 'R$_\oplus$')

# time		1 / s
s 		=	conv(1 , 's')
yr 		=	conv(1 / (60*60*24*365.25), 'yr')
Myr		=	conv(1 / (yr.value*1e6), 'Myr')
Gyr 	=	conv(1 / (yr.value*1e9), 'Gyr')

# mass 		1 / kg
kg 		=	conv(1, 'kg')
gr		=	conv(1000, 'g')
solarM 	=	conv(1 / 1.9891e30, 'M$_\odot$')
earthM 	= 	conv(1 / 5.9736e24, 'M$_\oplus')

# misc
solarI 	=	conv(1 / 1.365e3, 'S')
solarL 	= 	conv(1 / 3.839e26, 'L$_\odot$')


#===============================================================================
""" unit dictionary
- a function that returns a dictionary of constants with desired units
- cgs units are currently 'practical' units for electromagnetism
- plan to include gaussian functionality for cgs in the future"""
#-------------------------------------------------------------------------------

def unit_conversion(length=m,mass=kg,time=s):
	""" returns a dictionary of constants based off
	desired units of length, mass, and time.

	keyword parameters
	------------------
	length_object:	'conv' class object (value,units), default = m
	mass_object:	'conv' class object (value,units), default = kg
	time_object:	'conv' class object (value,units), default = s
	"""

	if all(( length == cm , mass == gr , time == s )):
		vel 	=	conv( length.value / time.value, '(cm/s)')
		acc 	=	conv( vel.value / time.value, 'Gal')
		force 	=	conv( mass.value * acc.value, 'dyn')
		energy 	=	conv( force.value * length.value, 'erg')
		power 	=	conv( energy.value / time.value, '(erg/s)')
		press 	=	conv( mass.value / length.value / time.value**2, 'Ba')
	elif all(( length == m , mass == kg , time == s )):
		vel 	=	conv( length.value / time.value, '(m/s)')
		acc 	=	conv( vel.value / time.value, '(m/s$^2$)')
		force 	=	conv( mass.value * acc.value, 'N')
		energy 	=	conv( force.value * length.value, 'J')
		power 	=	conv( energy.value / time.value, 'W')
		press 	=	conv( mass.value / length.value / time.value**2, 'Pa')
	else:
		vel 	=	conv( length.value / time.value, '(%s/%s)' % (length.units, time.units) )
		acc 	=	conv( vel.value / time.value, '(%s/%s$^2$)' % (length.units, time.units) )
		force 	=	conv( mass.value * acc.value, '(%s %s / %s$^2$)' % (mass.units, length.units, time.units) )
		energy 	=	conv( force.value * length.value, '(%s %s$^2$ / %s$^2$)' % (mass.units, length.units, time.units) )
		power 	=	conv( energy.value / time.value, '(%s %s$^2$ / %s$^3$)' % (mass.units, length.units, time.units) )
		press 	=	conv( mass.value / length.value / time.value**2, '%s / %s / %s$^2$' % (mass.units, length.units, time.units) )


	units 	=	{}

	# conversion factors
	units['velocity']		=	vel
	units['acceleration']	=	acc
	units['force']			=	force
	units['energy']			=	energy
	units['power']			=	power
	units['pressure']		=	press

	# conversions
	units['G']			=	conv(6.673e-11 * force.value * length.value**2 / mass.value**2, '%s %s$^2$ / %s$^2$' % (force.units, length.units, mass.units) )
	units['c']			=	conv(2.99792458e8 * vel.value, vel.units)
	units['mu_0']		=	conv(4*pi*1e-7 * force.value, '%s / A$^2$' % force.units)
	units['ep_0']		=	conv(1/units['mu_0'].value/units['c'].value**2, 'A$^2$ / %s / %s$^2$' % (force.units, vel.units) )
	units['e']			=	conv(1.602176462e-19, 'C')
	units['eV']			=	conv(1.602176462e-19 * energy.value, energy.units )
	units['h']			=	conv(6.62606876e-34 * energy.value * time.value, '%s %s' % (energy.units, time.units) )
	units['h_eV']		=	conv(4.13566727e-15 * time.value, 'eV %s' % time.units )
	units['hbar']		=	conv(units['h'].value/(2*pi), units['h'].units)
	units['hbar_eV']	=	conv(units['h_eV'].value/(2*pi), units['h_eV'].units)
	units['hc']			=	conv(units['h'].value * units['c'].value, '%s %s' % (units['h'].units, units['c'].units) )
	units['hc_eV']	 	=	conv(units['h_eV'].value * units['c'].value, '%s %s' % (units['h_eV'].units, units['c'].units) )
	units['k']			=	conv(1.3806503e-23 * energy.value, '%s / K' % energy.units)
	units['sigma']		=	conv(5.670400e-8 * power.value / length.value**2, '%s / %s$^2$ / K$^4$' % (power.units, length.units) )
	units['a']			=	conv(7.565767e-16 * energy.value / length.value**3, '%s / %s$^3$ / K$^4$' % (energy.units, length.units) )
	units['u']			=	conv(1.66053873e-27 * mass.value, mass.units)
	units['m_e']	 	=	conv(9.10938188e-31 * mass.value, mass.units)
	units['m_p']		=	conv(1.67262158e-27 * mass.value, mass.units)
	units['m_n']		=	conv(1.67492716e-27 * mass.value, mass.units)
	units['m_H']		=	conv(1.673532499e-27 * mass.value, mass.units)
	units['N_A']		=	conv(1000*6.02214199e23 / mass.value, '# / %s' % mass.units)
	units['R']			=	conv(1000*8.314472 * energy.value / mass.value, '%s / %s / K' % (energy.units, mass.units) )
	units['a_0']		=	conv(5.291772083e-11 * length.value, length.units)
	units['Ry_cons']	=	conv(1.0973731568549e7 / length.value, '1/%s' % length.units )
	units['Ry_energy']	=	conv(units['Ry_cons'].value * units['hc'].value, energy.units)
	units['Ry_eV']		=	conv(units['Ry_energy'].value/units['eV'].value, 'eV')

	# if all(( length == cm , mass == gr , time == s )):
	# 	units['e'] 		=

	return units
