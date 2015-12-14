	#
	#	Script for calculation of far-field pattern in the liquid case (prolate)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	# Physical parameters
	rho10 = 5.00 ; # Density ratio
	k_0 = 4 ; # Wave number in media 0
	k_1 = 6 ; # Wave number in media 1
	a = 1 ; # Major semiaxis spheroid
	b = 0.99 ; # Minor semiaxis spheroid
	theta_inc = pi/4 ; # Incidence angle in radians
	
	# Software parameters
	M = 4 ;
	method = 2 ;
	delta_eta = 0.005 ;
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ; 
	
	# Coefficient evaluation
	(Amn, Bmn) = Coeff_Pro( rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
 	
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidPro( k_0, a, b, M, method, eta_inc, Amn, delta_eta ) ;
	
	# Saving to disk
	writedlm( "Out.Pattern.dat", Pattern , '\t' ) ;			