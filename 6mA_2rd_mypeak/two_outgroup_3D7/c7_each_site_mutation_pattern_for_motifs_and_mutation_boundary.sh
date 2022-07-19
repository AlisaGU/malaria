#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
# motifs="GAACAA GAAGAA GAAGAT GACCAA GACGAA GACGAT GATGAA GATGAT"
# motifs="GAMSAA GAHGAW GAWGAW GAatGAat_1 GAatGAat_3 GAatGAat_4"
# motifs="GAatGAat_1 GAatGAat_2 GAatGAat_3 GAatGAat_4"
# motifs="GAASAA GAAGAA GAACAA random_seq"
# motifs="GACGAA GACCAA"
# motifs="motifs_summary"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# 计算motif各位点在基因组的位置
# for motif_seq in $(echo $motifs | tr " " "\n"); do
#     sed "s/motif_seq_template/$motif_seq/g" $code_dir/s7_each_site_genome_index_for_motifs.sh >$code_dir/s7_each_site_genome_index_for_motifs.sh_$motif_seq
#     sbatch $code_dir/s7_each_site_genome_index_for_motifs.sh_$motif_seq
# done

# 计算motif位点的突变速率
# for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
#     for motif_seq in $(echo $motifs | tr " " "\n"); do
#         sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_mutation_load_for_each_site_of_motif.sh | sed "s/motif_seq_template/${motif_seq}/g" >$code_dir/s7_compute_mutation_load_for_each_site_of_motif.sh_${popu_symbol}_${motif_seq}
#         sbatch $code_dir/s7_compute_mutation_load_for_each_site_of_motif.sh_${popu_symbol}_${motif_seq}
#     done
# done

# 计算flank site在基因组的位置，以及对应的突变速率
# for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
#     for motif_seq in $(echo $motifs | tr " " "\n"); do
#         sed "s/motif_seq_template/$motif_seq/g" $code_dir/s7_determine_mutation_boundary.sh | sed "s/popu_symbol_template/${popu_symbol}/g" >$code_dir/s7_determine_mutation_boundary.sh_$motif_seq
#         sbatch $code_dir/s7_determine_mutation_boundary.sh_$motif_seq
#     done
# done

# 计算A、T、C、G四种碱基在信号区的基底速率
# for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
#     sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh >$code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh_$popu_symbol
#     sbatch $code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh_$popu_symbol
# done

# 计算各位点碱基比例
# for motif_seq in $(echo $motifs | tr " " "\n"); do
#     sed "s/motif_template/$motif_seq/g" $code_dir/s7_compute_each_site_base_propor.sh >$code_dir/s7_compute_each_site_base_propor.sh_$motif_seq
#     sbatch $code_dir/s7_compute_each_site_base_propor.sh_$motif_seq
# done

#
#
#
#
#
#
#
#
# cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern
# cat $(find */motif_bed/ -name Pf*.peak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >all_motifs_position.peak.bed
# cat $(find */motif_bed/ -name Pf*.nopeak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >all_motifs_position.nopeak.bed

# cat $(find GAMSAA/motif_bed/ -name Pf*.peak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >GAMSAA.peak.bed
# cat $(find GAMSAA/motif_bed/ -name Pf*.nopeak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >GAMSAA.nopeak.bed

# cat $(find GAHGAW/motif_bed/ -name Pf*.peak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >GAHGAW.peak.bed
# cat $(find GAHGAW/motif_bed/ -name Pf*.nopeak.*) | awk '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t"$4}' | sort -k1,1 -k2,2n | uniq >GAHGAW.nopeak.bed
