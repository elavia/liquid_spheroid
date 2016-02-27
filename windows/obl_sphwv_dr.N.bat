	:: Bash script for calculation of 'dr' coefficients and normalization constant 'N'. The inputs are 'c','m','n' 'name' 
	:: Sintaxis:
	:: 	obl_sphwv_dr.N.sh	c	m	n name 
		
	
	:: Spheroidal oblate function parameters
	 call obl.parameters.bat
	
	:: Calculation
	obl_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w lambda
	obl_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w dr -n_dr %n_dr% -dr_min %dr_min%
	copy data\obl_%4_dr.txt Out_dr.dat	
	obl_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w dr_neg -n_dr_neg %n_dr_neg% -dr_neg_min %dr_neg_min%
	obl_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w N
	copy data\obl_%4_N.txt Out_S_N.dat	
		