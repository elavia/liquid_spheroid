	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Physical parameters
	rho10 =1000 ; # Density ratio
	k_0 = 5 ; # Wave number in media 0
	k_1 = 3 ; # Wave number in media 1
	a = 1 ; # Major semiaxis spheroid
	b = 0.25 ; # Minor semiaxis spheroid
	theta_inc = 2*pi/3 ; # Incidence angle in radians
	theta_incdeg= round(Int,180*theta_inc/pi); 
	# Software parameters
	M = 8 ;
	method = 2 ;
	delta_eta = 0.00390625 ;  # safe value: 1/2^8 
    
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conversion theta to eta
	eta_inc = cos( theta_inc ) ; 
	
	# Coefficient evaluation
	(Amn, Bmn) = Coeff_Obl( rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
	
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidObl(k_0, a, b, M, method, eta_inc, Amn, delta_eta ) ;
	
	# Saving to disk
	fileName = string("Obl_a_",a,"_b_",b,"_k0_",k_0,"_k1_",k_1,"_rho1rho0_",rho10,"_inc_",theta_incdeg,".dat");
	writedlm(fileName, Pattern , '\t' ) ;

	# Coefficient evaluation
 	
	 AmnD=Coeff_OblDN( k_0, a, b, M, method, eta_inc, "N")

	# Far-Field pattern calculation
	Pattern = Pattern_LiquidObl( k_0, a, b, M, method, eta_inc, AmnD, delta_eta ) ;

	# Saving to disk
	fileName = string("Obl_a_",a,"_b_",b,"_k0_",k_0,"_k1_",k_1,"_rho1rho0_",rho10,"_inc_",theta_incdeg,"Neumann.dat");
	writedlm(fileName, Pattern , '\t' ) ;
	
	# caso comun 
	rho10=3;
	
	# Coefficient evaluation
	(Amn, Bmn) = Coeff_Obl(rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
 	
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidObl(k_0, a, b, M, method, eta_inc, Amn, delta_eta ) ;
	
	# Saving to disk
	fileName = string("Obl_a_",a,"_b_",b,"_k0_",k_0,"_k1_",k_1,"_rho1rho0_",rho10,"_inc_",theta_incdeg,".dat");
	writedlm(fileName, Pattern , '\t' ) ;

	
	
	
	
	
	
	
	
	
	
	