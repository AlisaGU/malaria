#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err

# FUNCTIONS TODO:
get_two_outgroup_not_consistent_region_in_core() {
    # $code_dir/s2_1_get_two_outgroup_strain_inconsistent_region_in_core.R $inconsistent_region_in_core # 放弃

    # for chrom in $(echo $chroms | tr " " "\n"); do
    #     sed "s/chrom_template/$chrom/g" $code_dir/s2_get_two_outgroup_not_consistent_region_in_core.sh >$code_dir/s2_get_two_outgroup_not_consistent_region_in_core.sh_${chrom}
    #     sbatch $code_dir/s2_get_two_outgroup_not_consistent_region_in_core.sh_${chrom}
    # done

    cat $two_outgroup_working_dir/*.inconsistent_region_in_core.bed | sort -k1,1 -k2,2n >$inconsistent_region_in_core
    rm -rf *.inconsistent_region_in_core.bed
}

get_two_outgroup_consistent_region_in_core() {
    $bedtools complement -i $inconsistent_region_in_core -g $genome_size >$complementary_to_inconsistent_region_in_core
    $bedtools intersect -a $complementary_to_inconsistent_region_in_core -b $core_bed >$consistent_region_in_core
}

get_peak_nopeak_in_core_consistent_region() {
    $bedtools intersect -a $consistent_region_in_core -b $peak_bed >$core_consistent_peak
    $bedtools intersect -a $consistent_region_in_core -b $nopeak_bed >$core_consistent_nopeak

    cd $macs2_two_outgroup_consistent_dir
    awk -F "\t" '{print>$1}' $core_consistent_peak
    for file in Pf3D7_*_v3; do
        mv $file $file".peak.bed"
    done

    awk -F "\t" '{print>$1}' $core_consistent_nopeak
    for file in Pf3D7_*_v3; do
        mv $file $file".nopeak.bed"
    done
}

extract_peak_seq_in_core_consistent_region() {
    mkdir -p $core_peak_seq_dir
    cd $core_peak_seq_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        $seqkit subseq --bed $macs2_two_outgroup_consistent_dir/${autosome}.peak.bed $genome -o ${autosome}.peak.core.consistent.fa
    done

    cd -
}

locate_motif() {
    mkdir -p $motif_bed_dir
    cd $motif_bed_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_peak_seq_dir/${autosome}.peak.core.consistent.fa | $seqkit locate -i -d -f $motif_pattern --bed >${autosome}.motif.bed
    done

    cd -
}

get_block_bed_for_motif() {
    mkdir -p $region_bed_dir
    for bed in $(/usr/bin/ls $motif_bed_dir); do
        # $code_dir/s2_get_block_bed_for_motif.R ${motif_bed_dir}/$bed $region_bed_dir
        sed "s/bed_template/$bed/g" $code_dir/s2_get_block_bed_for_motif_one_chrom.sh >$code_dir/s2_get_block_bed_for_motif_one_chrom.sh_${bed}
        sbatch $code_dir/s2_get_block_bed_for_motif_one_chrom.sh_${bed}
    done
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
reference_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7"
anno_dir=${reference_dir}/anno
two_outgroup_working_dir=$anno_dir/3D7_distant_africa_strain_consistent_region_in_core
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
core_peak_seq_dir=$macs2_two_outgroup_consistent_dir/peak_seq
motif_bed_dir=$macs2_two_outgroup_consistent_dir/motif_loc
region_bed_dir=$macs2_two_outgroup_consistent_dir/region_loc

genome=${reference_dir}/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta
genome_size=${reference_dir}/genome/chrom.length
peak_bed=$macs2_out_dir/"peak.bed"
nopeak_bed=$macs2_out_dir/"nopeak.bed"
core_bed=$reference_dir/anno/3D7.core.bed
inconsistent_region_in_core=$two_outgroup_working_dir/inconsistent_region_in_core.bed
complementary_to_inconsistent_region_in_core=$two_outgroup_working_dir/complementary_to_inconsistent_region_in_core.bed
consistent_region_in_core=$two_outgroup_working_dir/consistent_region_in_core.bed
core_consistent_peak=$macs2_two_outgroup_consistent_dir/core_consistent_peak.bed
core_consistent_nopeak=$macs2_two_outgroup_consistent_dir/core_consistent_nopeak.bed
motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition_jiang/intergenic_motif/motif_pattern"

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
chroms="01 02 03 04 05 06 07 08 09 10 11 12 13 14"
bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# mkdir -p $two_outgroup_working_dir $macs2_two_outgroup_consistent_dir
# get_two_outgroup_not_consistent_region_in_core
# get_two_outgroup_consistent_region_in_core
get_peak_nopeak_in_core_consistent_region
extract_peak_seq_in_core_consistent_region
locate_motif
get_block_bed_for_motif
