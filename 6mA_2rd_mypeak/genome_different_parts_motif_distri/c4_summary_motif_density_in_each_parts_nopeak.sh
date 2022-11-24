#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_average_motif_count_and_base_count() {
    density_file=$motif_density_dir/nopeak_motif_density
    rm -rf $density_file

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for region in $(echo "intergenic TSS_2KB TTS_2KB genes exons introns" | tr " " "\n"); do
            region_bed=$genome_parts_dir/${region}.bed
            region_bed_chrom=$(awk -v chrom=$chrom '$1==chrom' $region_bed | $bedtools intersect -a stdin -b $chrom_core_genome_two_outgroup_nopeak_dir/${chrom}.nopeak.bed)

            motif_base_count=$(echo "$region_bed_chrom" | $bedtools intersect -a $motif_core_genome_two_outgroup_nopeak_dir/${chrom}.nopeak.motif.right.bed -b stdin | sort -k1,1 -k2,2n | $bedtools merge -i stdin | awk '{print $3-$2}' | awk '{s+=$1}END{print s}')

            if [ $motif_base_count -gt 0 ] 2>/dev/null; then
                echo ""
            else
                motif_base_count=0
            fi

            region_len=$(echo "$region_bed_chrom" | awk '{s+=$3-$2}END{print s}')
            if [ $region_len -gt 0 ] 2>/dev/null; then
                echo ""
            else
                region_len=0
            fi

            value+=" $(echo $motif_base_count" "$region_len)"
        done
        echo $value >>$density_file
    done
}
# VARIABLE NAMING TODO:
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri"
motif_density_dir=$global_output_dir/motif_density_diff_denominator_core_genome_two_outgroup
genome_parts_dir=$global_output_dir/genome_parts
motif_core_genome_two_outgroup_nopeak_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_GAWGAW_nopeak/motif_loc"
chrom_core_genome_two_outgroup_nopeak_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
mkdir -p $motif_density_dir
cd $motif_density_dir
compute_average_motif_count_and_base_count
