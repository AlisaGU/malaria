#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/gene_6mA_density_signal/KD"
genome_size="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length"

macs2="/home/gushanshan/anaconda3/envs/vscode_r/bin/macs2"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
# $macs2 callpeak -t 6mAKD-T3_ChIP.bam -c 6mAKD-T3_Input.bam -n 6mAKD -B --SPMR -f "BAM" -g 2.3e7 --extsize 147 --nomodel --call-summits --outdir $working_dir -p 0.05 --bw 50
$bedtools unionbedg -i 6mAKD_treat_pileup.bdg 6mAKD_control_lambda.bdg -empty -g $genome_size -header -names treat control | awk 'NR>1{fe=($4*8.323738+1)/($5*8.323738+1);print $1"\t"$2"\t"$3"\t"fe}' >merged.bed
