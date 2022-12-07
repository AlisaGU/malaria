#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
histone_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/histone_methy"
histone_type="3D7-T3-H3K9Ac"
macs2="/home/gushanshan/anaconda3/envs/vscode_r/bin/macs2"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
cd $histone_global_dir/$histone_type
$macs2 callpeak -t ${histone_type}_ChIP.bam -c ${histone_type}_Input.bam -n ${histone_type} --nomodel -f "BAM" -B --SPMR --extsize 147 --broad -g 2.3e7 --broad-cutoff 0.1 --outdir ./BAM_broad_nomodel_dup1/

# $macs2 callpeak -t ${histone_type}_ChIP.bam -c ${histone_type}_Input.bam -n ${histone_type} --nomodel -f "BAM" -B --SPMR --extsize 147 -g 2.3e7 --call-summits --outdir ./BAM_nomodel_dup1/

$macs2 callpeak -t ${histone_type}_ChIP.bam -c ${histone_type}_Input.bam -n ${histone_type} --keep-dup all --nomodel -f "BAM" -B --SPMR --extsize 147 --broad -g 2.3e7 --broad-cutoff 0.1 --outdir ./BAM_broad_nomodel_dupall/

# $macs2 callpeak -t ${histone_type}_ChIP.bam -c ${histone_type}_Input.bam -n ${histone_type} --keep-dup all --nomodel -f "BAM" -B --SPMR --extsize 147 -g 2.3e7 --call-summits --outdir ./BAM_nomodel_dupall/

for i in $(find ./ -name *xls); do
    grep "^Pf3D7" $i | awk '{print $1"\t"$2"\t"$3}' >$(echo $i | sed "s/xls/bed/g")
done
