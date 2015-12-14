	# !/bin/bash

	# Shell script for calculation of 'dr' coefficients and normalization constant 'N'.
	# The inputs are 'c','m','n','name' (variable part of the name)
	# Sintaxis:
	# 	pro_sphwv_dr.N.sh	c	m	n	name
	#				$1	$2	$3	$4
	
	# Spheroidal prolate function parameters
	source ./pro.parameters	

	# Calculation 
 	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w lambda
 	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr -n_dr $n_dr -dr_min $dr_min
   	cp data/pro_"$4"_dr.txt ./Out_dr.dat
 	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr_neg -n_dr_neg $n_dr_neg -dr_neg_min $dr_neg_min
  	./pro_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w N
 	cp data/pro_"$4"_N.txt ./Out_S_N.dat