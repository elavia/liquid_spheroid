	#
	#	Script for calculation of far-field pattern in the Dirichlet case (oblate)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	# Physical parameters
	k_0 = 4 ; # Wave number in media 0
	a = 1 ; # Major semiaxis spheroid
	b = 0.8 ; # Minor semiaxis spheroid
	theta_inc = pi/4 ; # Incidence angle in radians
	
	# Software parameters
	M = 8 ;
	method = 2 ;
	delta_eta = 0.005 ;
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ; 
	# Coefficient evaluation
	Amn = Coeff_OblDN( k_0, a, b, M, method, eta_inc, "D" ) ;
 		
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidObl( k_0, a, b, M, method, eta_inc, Amn, delta_eta ) ;
	
	# Saving to disk
	writedlm( "Out.Pattern.dat", Pattern , '\t' ) ;	
	