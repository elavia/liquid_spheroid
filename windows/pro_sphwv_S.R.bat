	:: BashShell script for calculation of S and R in one point. The inputs are 'c', 'm', 'n', 'eta', 'xi' 
	:: Sintaxis:
	:: 	pro_sphwv_S.R.sh	c	m	n	eta	xi
	::				%1	%2	%3	%4	%5
	
	:: Spheroidal prolate function parameters
	
	call pro.parameters.bat

	:: Preliminaries
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w lambda
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w dr -n_dr %n_dr% -dr_min %dr_min%
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w dr_neg -n_dr_neg %n_dr_neg% -dr_neg_min %dr_neg_min%
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w N
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w F
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w k1
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w k2
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w c2k -n_c2k %n_c2k% -c2k_min %c2k_min%
	:: Spheroidal wave function S
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w S1 -a %4 -b %4 -d 0.01 -arg_type eta -p 34>Out_S.dat
	:: Spheroidal wave function R
	pro_sphwv -max_memory %max_mem% -precision %prec% -verbose n -c %1 -m %2 -n %3 -w R -a %5 -b %5 -d 0.1 	-arg_type xi -which R1_1,R1_2,R2_1,R2_2 -p 34>Out_R.dat
	
	