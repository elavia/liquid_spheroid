	# !/bin/bash

	# Shell script for calculation of S (grid) and S incident.
	# The inputs are 'c','m','n','eta_inc','grid_step'
	# Sintaxis:
	#	pro_sphwv_Sgrid.Sinc.sh	c	m	n	eta_inc		grid_step
	#				$1	$2	$3	$4		$5		

	# Spheroidal prolate function parameters
	source ./pro.parameters	
		
	# Preliminaries
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w lambda
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr -n_dr $n_dr -dr_min $dr_min
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr_neg -n_dr_neg $n_dr_neg -dr_neg_min $dr_neg_min
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w N
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w F
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w k1
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w k2
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w c2k -n_c2k $n_c2k -c2k_min $c2k_min	
 	
 	# Spheroidal wave function S incident
	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w S1 -a $4 -b $4 -d 0.01 \
 	-arg_type eta -p 34
	cp Salida_S.dat Out_Sinc.dat
	
 	# Spheroidal wave function S
 	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w S1 -a -1.0 -b 1.0000001 -d $5 \
 	-arg_type eta -p 34 
 	cp Salida_S.dat Out_S.dat
