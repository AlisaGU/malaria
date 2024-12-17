#!/bin/bash
#SBATCH -n 10
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/paper_ZbynekBozdech_Tstage"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/gene_6mA_density_signal/drug_resistance"

prefixs="s16_S1A_6AR_vs_6AS.up s16_S1A_6AR_vs_6AS.down s16_S1B_11CR_vs_11CS.up s16_S1B_11CR_vs_11CS.down s17_S2A_6AR_vs_6AS_ART.up s17_S2A_6AR_vs_6AS_ART.down s17_S2B_11CR_vs_11CS_ART.up s17_S2B_11CR_vs_11CS_ART.down"
WT_chip="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/3D7-T3_ChIP.bam"
WT_input="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/3D7-T3_Input.bam"
KD_chip="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/6mAKD-T3_ChIP.bam"
KD_input="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD/6mAKD-T3_Input.bam"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:

for prefix in $(echo $prefixs | tr " " "\n"); do
    # sed -e "s/prefix_template/$prefix/g" -e "s@working_dir_template@$working_dir@g" -e "s@WT_chip_template@$WT_chip@g" -e "s@WT_input_template@$WT_input@g" -e "s@KD_chip_template@$KD_chip@g" -e "s@KD_input_template@$KD_input@g" $code_dir/s1.computeDensity.template.sh >$code_dir/s1.computeDensity.$prefix.sh
    # sbatch $code_dir/s1.computeDensity.$prefix.sh

    $code_dir/s2_get_line_Tstage.R $working_dir $prefix
done
