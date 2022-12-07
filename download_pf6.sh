#!/bin/bash
#SBATCH -n 10
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
axel="/home/gushanshan/anaconda3/envs/vscode_r/bin/axel"
out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
ftp_prefix="ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf/"
for chro in `echo "02 03 04 05 06 07 08 09 10 11 12 13 14"|tr " " "\n"`
do
    $axel -n 10 -o $out_dir $ftp_prefix/"Pf_60_public_Pf3D7_"$chro"_v3.final.vcf.gz"
    $axel -n 10 -o $out_dir $ftp_prefix/"Pf_60_public_Pf3D7_"$chro"_v3.final.vcf.gz.md5"
    $axel -n 10 -o $out_dir $ftp_prefix/"Pf_60_public_Pf3D7_"$chro"_v3.final.vcf.gz.tbi"
done