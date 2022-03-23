#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
extract_vcf() {
    for chrom in $(echo "01 02 03 04 05 06 07 08 09 10 11 12 13 14" | tr " " "\n"); do
        $code_dir/s4_extract_vcf.sh $chrom
    done
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_narrow_peak_out=$macs2_out_dir/"3D7_2rd_peaks.narrowPeak"
ref_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7"
chrom_len=$ref_dir/"genome/chrom.length"
ref_anno=$ref_dir/"anno/PlasmoDB-36_Pfalciparum3D7.gff"

peak_bed=$macs2_out_dir/"peak.bed"
nopeak_bed=$macs2_out_dir/"nopeak.bed"
genes_bed=$ref_dir/"anno"/genes.bed
promoters_bed=$ref_dir/"anno"/promoters.bed
nonintergenic_bed=$ref_dir/"anno"/nonintergenic.bed
intergenic_bed=$ref_dir/"anno"/intergenic.bed
intergenic_peak_bed=$macs2_out_dir/"intergenic_peak.bed"
intergenic_nopeak_bed=$macs2_out_dir/"intergenic_nopeak.bed"
intergenic_nopeak_closest_to_peak_bed=$macs2_out_dir/"intergenic_nopeak_closest_to_peak.bed"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
bedmap="/picb/evolgen/users/gushanshan/software/bedops/bin/bedmap"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# $bedtools complement  -i $peak_bed -g $chrom_len -L > $nopeak_bed # 非peak区位置
# $code_dir/s4_create_bed_for_genomic_features.R $ref_anno $ref_dir/"anno"
# $bedtools flank -i $genes_bed -g $chrom_len -l 200 -r 0 -s > $promoters_bed
# cat $promoters_bed $genes_bed | $bedtools sort -i > $nonintergenic_bed
# $bedtools complement -i $nonintergenic_bed -g $chrom_len -L > $intergenic_bed

# $bedtools intersect -a $peak_bed -b $intergenic_bed | awk '{name="peak"NR;print $1"\t"$2"\t"$3"\t",name,"\t",".","\t","."}' | sort -k1,1 -k2,2n >$intergenic_peak_bed
# $bedmap --echo --delim '\t' peak_bed $intergenic_bed >$macs2_out_dir/a
# $bedtools intersect -a $nopeak_bed -b $intergenic_bed | awk '{name="nopeak"NR;print $1"\t"$2"\t"$3"\t",name,"\t",".","\t","."}' | sort -k1,1 -k2,2n >$intergenic_nopeak_bed

# $bedtools closest -a $intergenic_peak_bed -b $intergenic_nopeak_bed -io -k 10 -d | awk '$13<=10000{print $0"\t"$9-$8}' >$macs2_out_dir/a
# $code_dir/s4_choose_best_nopeak.R $macs2_out_dir/a $macs2_out_dir/b
# sort -k1,1 -k2,2n $macs2_out_dir/b >$intergenic_nopeak_closest_to_peak_bed
# rm -rf $macs2_out_dir/a $macs2_out_dir/b

extract_vcf
