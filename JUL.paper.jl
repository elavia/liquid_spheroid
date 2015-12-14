

# 	using Base.LinAlg.BLAS ;
# 	using Base.LinAlg.LAPACK ;
# 	include("JUL.obladas.LIQ.paper.jl") ;
# 
# 	function Paper( cuerpo::AbstractString )
# 		
# 		const rho10 = 1000 ; # 1/1000 ;
# 		const k_0 = 5 ;
# 		const k_1 = 3  ;
# 		const a = 1 ;
# 		const b = 0.25 ; # 0.99 ; 
# 		const M = 12 ;
# 		const method = 2 ;
# 		const delta_eta = 0.005 ;
# 		
# 		# Las incidencias son: (1 , pi) y (0 , pi*3/2) o bien ( -1 , 0)
# 		
# 		const eta_0 = -0.5 ; # 0 ; # 1 ;
# 		const theta_inc = pi*2/3 ; # pi*3/2 ; #  pi ; # pi*3/2 ;
# 		
# 		if cuerpo == "prolate"
# 			# Coefficient calculation
# #  			(Amn, Bmn) = Coeff_pro( rho10, k_0, k_1, a, b, M, method, eta_0 ) ;
# # 			Amn = Coeff_pro_DN( k_0, a, b, M, method, eta_0, "D" ) ;
#  			Amn = Coeff_pro_DN( k_0, a, b, M, method, eta_0, "N" ) ;
#  			Bmn = Amn ;
# 			# Far-Field pattern calculation
# 			Pattern = PatternProlatePaper( rho10, k_0, k_1, a, b, M, method, eta_0, Amn, Bmn, delta_eta ) ;
# 			# Saving matrix to disk
# 			writedlm( "Out.Pattern.dat", Pattern , '\t' );
# 		elseif cuerpo == "sphere"
# 			Pattern = PatternSpherePaper( rho10, k_0, k_1, a, 200, 150, theta_inc ) ;
# 			# Saving matrix to disk
# 			writedlm( "Out.Pattern.dat", Pattern , '\t' );				
# 		elseif cuerpo == "oblate"
# 			# Coefficient calculation
# #  			(Amn, Bmn) = Coefficientes( rho10, k_0, k_1, a, b, M, method, eta_0 ) ;
# #  			Amn = Coeff_obl_DN( k_0, a, b, M, method, eta_0, "D" ) ;
#   			Amn = Coeff_obl_DN( k_0, a, b, M, method, eta_0, "N" ) ;
#   			Bmn = Amn ;			
# 			# Far-Field pattern calculation
# 			Pattern = PatternOblatePaper( rho10, k_0, k_1, a, b, M, method, eta_0, Amn, Bmn, delta_eta ) ;
# 			# Saving matrix to disk
# 			writedlm( "Out.Pattern.dat", Pattern , '\t' );
# 		else 
# 			println("Error. Ninguno lo satisfizo") ;
# 		end 
# 		
# 	end


# 	using Base.LinAlg.BLAS ;
# 	using Base.LinAlg.LAPACK ;
# 	include("JUL.obladas.LIQ.paper.jl") ;

	function Paper( cuerpo::AbstractString )
		
		const rho10 = 1000 ; # 1/1000 ;
		const k_0 = 5 ;
		const k_1 = 3  ;
		const a = 1 ;
		const b = 0.75 ; # 0.99 ; 
		const M = 8 ;
		const method = 2 ;
		const delta_eta = 0.005 ;
		
		# Las incidencias son: (1 , pi) y (0 , pi*3/2) o bien ( -1 , 0)
		
# 		const eta_0 = -0.5 ; # 0 ; # 1 ;
		const theta_inc = 4.5 - pi ; # pi*2/3 ; # pi*3/2 ; #  pi ; # pi*3/2 ;
		const eta_inc = cos( theta_inc ) ;
		
		if cuerpo == "prolate"
			# Coefficient calculation
  			(Amn, Bmn) = Coeff_Pro( rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
			# Far-Field pattern calculation
			Pattern = Pattern_LiquidPro( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );
		elseif cuerpo == "prolateD"
			# Coefficient calculation
			Amn = Coeff_ProDN( k_0, a, b, M, method, eta_inc, "D" ) ;
			Pattern = Pattern_LiquidPro( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );
		elseif cuerpo == "prolateN"
			# Coefficient calculation
			Amn = Coeff_ProDN( k_0, a, b, M, method, eta_inc, "N" ) ;
			Pattern = Pattern_LiquidPro( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );
		elseif cuerpo == "sphere"
			Pattern = Pattern_LiquidSphere( rho10, k_0, k_1, a, 200, 150, theta_inc ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );				
		elseif cuerpo == "oblate"
			# Coefficient calculation
  			(Amn, Bmn) = Coeff_Obl( rho10, k_0, k_1, a, b, M, method, eta_inc ) ;
			# Far-Field pattern calculation
			Pattern = Pattern_LiquidObl( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );
		elseif cuerpo == "oblateD"
			# Coefficient calculation
			Amn = Coeff_OblDN( k_0, a, b, M, method, eta_inc, "D" ) ;
			Pattern = Pattern_LiquidObl( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );
		elseif cuerpo == "oblateN"
			# Coefficient calculation
			Amn = Coeff_OblDN( k_0, a, b, M, method, eta_inc, "N" ) ;
			Pattern = Pattern_LiquidObl( k_0, k_1, a, b, M, method, eta_inc, Amn, delta_eta ) ;
			# Saving matrix to disk
			writedlm( "Out.Pattern.dat", Pattern , '\t' );			
		else 
			println("Error. Ninguno lo satisfizo") ;
		end 
		
	end