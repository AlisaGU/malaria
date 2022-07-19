#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
validate_6mA_sliding_window_distribution_in_3th_data() {
    mkdir -p sliding_window
    cd sliding_window
    m6mA_distri="sliding_window_6mA_distri${prefix}"
    # m6mA_distri="sliding_window_6mA_distri_ATcountScaled${prefix}"

    rm -rf $m6mA_distri
    for chrom in $(echo $autosomes | tr " " "\n"); do
        regions=$m6mA_2rd_whole_genome_nopeak_dir/${chrom}.nopeak.two_outgroup_consistent_region.bed
        regions+=" $m6mA_2rd_whole_genome_nopeak_dir/${chrom}.nopeak500.two_outgroup_consistent_region.bed"
        for w in $(seq 1 10); do
            regions+=" $m6mA_2rd_whole_genome_peak_sliding_window_dir/$chrom/w${w}.two_outgroup_consistent_region.bed"
        done

        value=""
        for region in $(echo $regions | tr " " "\n"); do
            value+=" $(get_info $working_dir/$m6mA_3th_bed $region)"
        done
        echo $value >>$m6mA_distri
    done
    cd -
}

get_info() {
    a=$1
    b=$2
    m6mA_3th_count=$($bedtools intersect -a $a -b $b | wc -l)
    if [ $m6mA_3th_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        m6mA_3th_count=0
    fi
    region_len=$(awk '{print $3-$2}' $b | awk '{s+=$1}END{print s}')
    # region_len=$($seqkit subseq --quiet --bed $b $genome | $seqkit fx2tab -i -H -n -C "A" -C "T" | grep -v "^#" | awk '{a+=$2;t+=$3}END{print a+t}')
    echo $m6mA_3th_count $region_len
}

validate_6mA_motif_flank_pattern_in_3th_data() {
    mkdir -p motif_flank_pattern
    cd motif_flank_pattern

    m6mA_distri="motif_flank_pattern_6mA_distri${prefix}"
    # m6mA_distri="motif_flank_pattern_6mA_distri_ATcountScaled${prefix}"

    rm -rf $m6mA_distri
    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for site_number in $(seq 1 6); do
            value+=" $(get_motif_3th_6mA_count $working_dir/$m6mA_3th_bed site$site_number $chrom)"
        done

        for flank_number in $(seq 1 10); do
            value+=" $(get_info $working_dir/$m6mA_3th_bed $motif_flank_dir/flank_nomatter_strand.remove_close_overlap/flank_${flank_number}/${chrom}.peak.bed)"
        done
        echo $value >>$m6mA_distri
    done
    cd -
}

get_motif_3th_6mA_count() {
    m6mA_3th=$1
    site=$2
    chrom=$3

    m6mA_3th_count=$(cat $motif_flank_dir/$site/${chrom}.peak.* | sort -k1,1 -k2,2n | $bedtools intersect -a $m6mA_3th -b stdin | wc -l)
    if [ $m6mA_3th_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        m6mA_3th_count=0
    fi
    region_len=$(cat $motif_flank_dir/$site/${chrom}.peak.* | wc -l)

    # cat $motif_flank_dir/$site/${chrom}.peak.* >intermediate
    # region_len=$($seqkit subseq --quiet --bed intermediate $genome | $seqkit fx2tab -i -H -n -C "A" -C "T" | grep -v "^#" | awk '{a+=$2;t+=$3}END{print a+t}')
    # rm -rf intermediate

    echo $m6mA_3th_count $region_len
}
# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public"
# m6mA_3th_bed="3D7-6mA-9Cell.bed"

m6mA_2rd_whole_genome_nopeak_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/histone_methy/6mA/BAM_broad_nomodel_dup1/chrom"
m6mA_2rd_whole_genome_peak_sliding_window_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/histone_methy/6mA/BAM_broad_nomodel_dup1/chrom/sliding_window/relative.0.05"
motif_flank_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern/GAWGAW"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
cd $working_dir
prefix=""
m6mA_3th_bed="3D7-6mA-9Cell${prefix}.bed"
validate_6mA_sliding_window_distribution_in_3th_data
validate_6mA_motif_flank_pattern_in_3th_data

for prefix in $(echo "_COVgt25,_COVgt50,_FRACgt20_COVgt25,_QVgt20_COVgt10" | tr "," "\n"); do
    m6mA_3th_bed="3D7-6mA-9Cell${prefix}.bed"
    validate_6mA_sliding_window_distribution_in_3th_data
    validate_6mA_motif_flank_pattern_in_3th_data
done
