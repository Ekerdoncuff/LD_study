#!/bin/sh

# Made to run D_stats_study.c with different output name and parameters
# Requires D_stats_study.c D_etude_lib.c and D_etude_lib.h

# Create executable
gcc D_stats_study_double_treeweight.c D_etude_lib.c -o D_etude_study_treeweight.exe

# Define output file
outfile="d8_stats_outfile_n10_ite10000_treeweight.txt"

# Create header for outfile
echo "stats_D\tN\tn\tTCHANG\tCONTRAC\tmean_D\tvar_D\tmean_Dp\tvar_Dp\tmean_r\tvar_r\tmean_r2\tvar_r2_global\tmean_sd\tmean_sd_global_newweigth\tmean_CP0_global\tvar_CP0_global\tmean_CP1_global\tvar_CP1_global\tmean_CP2_global\tvar_CP2_global\tmean_CP3_global\tvar_CP3_global" > ${outfile}

# Add one line in the outfile with the chosen parameters
#./D_etude_study.exe -N 10 -T 0 -C 1 -A 0 -D 0 -n 100 -R 0 ${outfile}

# Example of loop on different strength of contraction
for log_contract in $(seq -2 0.01 2)
do
    contract=$(echo "e(${log_contract}*l(10))" | bc -l)
    echo $contract
    ./D_etude_study_treeweight.exe -N 10 -T 0.5 -C ${contract} -A 0 -D 8 -n 10000 -R 0 ${outfile}
done

###
# For the different type of recombination study:
###

#recombi_outfile="rec_d_stats_outfile.txt"

# Create headers for recombi outputs
#echo "stats_D\tN\tn\tTCHANG\tCONTRAC\tmean_D\tvar_D\tmean_Dp\tvar_Dp\tmean_r\tvar_r\tmean_r2\tvar_r2_global" > ${recombi_outfile}
#echo "N\tTCHANG\tCONTRAC\tpercentage_type0\tpercentage_type1\tpercentage_type2\tpercentage_type3" > type_recombi_${recombi_outfile}

# Headers for different recombi types
#echo "stats_D\t0\tN\tnb_recombi_type0\tTCHANG\tCONTRAC\tmean_D_tr0\tvar_D_tr0\tmean_D_p_tr0\tvar_D_p_tr0\tmean_r_tr0\tvar_r_tr0\tmean_r2_tr0\tvar_r2_tr0" > recombi_type0_${recombi_outfile}
#echo "stats_D\t1\tN\tnb_recombi_type1\tTCHANG\tCONTRAC\tmean_D_tr1\tvar_D_tr1\tmean_D_p_tr1\tvar_D_p_tr1\tmean_r_tr1\tvar_r_tr1\tmean_r2_tr1\tvar_r2_tr1" > recombi_type1_${recombi_outfile}
#echo "stats_D\t2\tN\tnb_recombi_type2\tTCHANG\tCONTRAC\tmean_D_tr2\tvar_D_tr2\tmean_D_p_tr2\tvar_D_p_tr2\tmean_r_tr2\tvar_r_tr2\tmean_r2_tr2\tvar_r2_tr2" > recombi_type2_${recombi_outfile}
#echo "stats_D\t3\tN\tnb_recombi_type3\tTCHANG\tCONTRAC\tmean_D_tr3\tvar_D_tr3\tmean_D_p_tr3\tvar_D_p_tr3\tmean_r_tr3\tvar_r_tr3\tmean_r2_tr3\tvar_r2_tr3" > recombi_type3_${recombi_outfile}

# Run with recombination type parameters
#./D_etude_study.exe -N 10 -T 0 -C 1 -A 0 -D 1 -n 100 -R 1 ${recombi_outfile}

# Loop on different strength of contraction
#for log_contract in $(seq -2 1 2)
#do
#    contract=$(echo "10^${log_contract}" | bc -l)
#    ./D_etude_study.exe -N 10 -T 0.5 -C ${contract} -A 0 -D 0 -n 100 -R 0 ${outfile}
#done
