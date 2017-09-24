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

	const 	=	{}
	units 	=	{}

	# conversion factors
	const['velocity']		=	vel.value
	const['acceleration']	=	acc.value
	const['force']			=	force.value
	const['energy']			=	energy.value
	const['power']			=	power.value
	const['pressure']		=	press.value

	units['velocity']		=	vel.units
	units['acceleration']	=	acc.units
	units['force']			=	force.units
	units['energy']			=	energy.units
	units['power']			=	power.units
	units['pressure']		=	press.units

	# conversions
	const['G']			=	6.673e-11 * force.value * length.value**2 / mass.value**2
	const['c']			=	2.99792458e8 * vel.value
	const['mu_0']		=	4*pi*1e-7 * force.value
	const['ep_0']		=	1 / const['mu_0']/const['c']**2
	const['e']			=	1.602176462e-19
	const['eV']			=	1.602176462e-19 * energy.value
	const['h']			=	6.62606876e-34 * energy.value * time.value
	const['h_eV']		=	4.13566727e-15 * time.value
	const['hbar']		=	const['h'] / (2*pi)
	const['hbar_eV']	=	const['h_eV'] / (2*pi)
	const['hc']			=	const['h'] * const['c']
	const['hc_eV']	 	=	const['h_eV'] * const['c']
	const['k']			=	1.3806503e-23 * energy.value
	const['sigma']		=	5.670400e-8 * power.value / length.value**2
	const['a']			=	7.565767e-16 * energy.value / length.value**3
	const['u']			=	1.66053873e-27 * mass.value
	const['m_e']	 	=	9.10938188e-31 * mass.value
	const['m_p']		=	1.67262158e-27 * mass.value
	const['m_n']		=	1.67492716e-27 * mass.value
	const['m_H']		=	1.673532499e-27 * mass.value
	const['m_He']		=	6.6464764e-27 * mass.value
	const['N_A']		=	1000*6.02214199e23 / mass.value
	const['R']			=	1000*8.314472 * energy.value / mass.value
	const['a_0']		=	5.291772083e-11 * length.value
	const['Ry_cons']	=	1.0973731568549e7 / length.value
	const['Rydberg']	=	const['Ry_cons'] * const['hc']
	const['Ry_eV']		=	const['Rydberg']/const['eV']

	units['G']			=	'%s %s$^2$ / %s$^2$' % (force.units, length.units, mass.units)
	units['c']			=	vel.units
	units['mu_0']		=	'%s / A$^2$' % force.units
	units['ep_0']		=	'A$^2$ / %s / %s$^2$' % (force.units, vel.units)
	units['e']			=	'C'
	units['eV']			=	energy.units
	units['h']			=	'%s %s' % (energy.units, time.units)
	units['h_eV']		=	'eV %s' % time.units
	units['hbar']		=	units['h']
	units['hbar_eV']	=	units['h_eV']
	units['hc']			=	'%s %s' % (units['h'], units['c'])
	units['hc_eV']	 	=	'%s %s' % (units['h_eV'], units['c'])
	units['k']			=	'%s / K' % energy.units
	units['sigma']		=	'%s / %s$^2$ / K$^4$' % (power.units, length.units)
	units['a']			=	'%s / %s$^3$ / K$^4$' % (energy.units, length.units)
	units['u']			=	mass.units
	units['m_e']	 	=	mass.units
	units['m_p']		=	mass.units
	units['m_n']		=	mass.units
	units['m_H']		=	mass.units
	units['m_He']		=	mass.units
	units['N_A']		=	'# / %s' % mass.units
	units['R']			=	'%s / %s / K' % (energy.units, mass.units)
	units['a_0']		=	length.units
	units['Ry_cons']	=	'1/%s' % length.units
	units['Rydberg']	=	energy.units
	units['Ry_eV']		=	'eV'

	return const, units

#===============================================================================
""" eV dictionary """
#-------------------------------------------------------------------------------

def eV_dictionary(Escale='eV', Tscale='s' , Lscale='m'):
	""" write dictionary for energy units

	keyword parameters
	------------------
	Escale:	decide what scaling for eV: 	( eV , keV , MeV , GeV  )
	Tscale: decide what scaling for time: 	( s , ms , micros , ns )
	Mscale: decide what scaling for length: ( m , mm , microm , nm )
	"""

	if Escale == 'eV':
		Efactor = 1
	elif Escale == 'keV':
		Efactor = 1 / 1e3
	elif Escale == 'MeV':
		Efactor = 1 / 1e6
	elif Escale == 'GeV':
		Efactor = 1 / 1e9

	if Tscale == 's':
		Tfactor = 1
	elif Tscale == 'ms':
		Tfactor = 1 / 1e-3
	elif Tscale == 'micros':
		Tfactor = 1 / 1e-6
	elif Tscale == 'ns':
		Tfactor = 1 / 1e-9

	if Lscale == 'm':
		Lfactor = 1
	elif Lscale == 'mm':
		Lfactor = 1 / 1e-3
	elif Lscale == 'microm':
		Lfactor = 1 / 1e-6
	elif Lscale == 'nm':
		Lfactor = 1 / 1e-9

	const 	=	{}
	units 	=	{}

	const['c'] 		=	299792458 * Lfactor / Tfactor**2
	const['h'] 		=	4.135667662e-15 * Efactor * Tfactor
	const['hbar']	=	6.582119514e10-16 * Efactor * Tfactor
	const['hc']		=	const['h'] * const['c']
	const['u']		=	931.494061e6 * Efactor / (Lfactor/Tfactor**2)
	const['m_e']	=	0.51100e6 * Efactor / (Lfactor/Tfactor**2)
	const['m_p']	=	938.27 * Efactor / (Lfactor/Tfactor**2)
	const['m_n']	=	939.57 * Efactor / (Lfactor/Tfactor**2)
	const['m_H']	=	938.790 * Efactor / (Lfactor/Tfactor**2)
	const['m_He']	=	4.00260325 * const['u']

	units['c']		=	'%s / %s$^2$' % ( Lscale , Tscale )
	units['h'] 		=	'%s %s' % ( Escale , Tscale )
	units['hbar']	=	'%s %s' % ( Escale , Tscale )
	units['hc']		=	'( %s ) ( %s )' % ( units['h'] , units['c'] )
	units['u']		=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )
	units['m_e']	=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )
	units['m_p']	=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )
	units['m_n']	=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )
	units['m_H']	=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )
	units['m_He']	=	'( %s ) / ( %s )$^2$' % ( Escale , units['c'] )

	return const, units

#===============================================================================
""" completed dictionaries """
#-------------------------------------------------------------------------------

SI,SI_units 	=	unit_conversion()
CGS,CGS_units	=	unit_conversion(length=cm,mass=gr,time=s)
eV,eV_units 	=	unit_conversion()
