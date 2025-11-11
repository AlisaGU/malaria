#!/bin/bash
#SBATCH -n 30
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/conserved_nonconserved_6mA"
bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd"
variableGenomicRegion="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7.noncore.bed"
conservedGenomicRegion="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7.core.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd $working_dir
$bedtools intersect -a $variableGenomicRegion -b $bam_dir/3D7-T3_ChIP.bam -c >variableGenomicRegion.rawReadCount.chip &

$bedtools intersect -a $variableGenomicRegion -b $bam_dir/3D7-T3_Input.bam -c >variableGenomicRegion.rawReadCount.input &

$bedtools intersect -a $conservedGenomicRegion -b $bam_dir/3D7-T3_ChIP.bam -c >conservedGenomicRegion.rawReadCount.chip &

$bedtools intersect -a $conservedGenomicRegion -b $bam_dir/3D7-T3_Input.bam -c >conservedGenomicRegion.rawReadCount.input &
