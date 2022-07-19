#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_whole_genome_gene_and_mRNA() {
    grep -v "^#" $anno_gtf | awk '$3=="gene"{print}' | sort -k1,1 -k4,4n | awk '{split($9,a,";");split(a[1],b,"=");s=$4-1;print $1"\t"s"\t"$5"\t"b[2]"\t"$7}' >$gene_bed
    grep -v "^#" $anno_gtf | awk '$3=="mRNA"{print}' | awk '{split($9,a,";");split(a[1],b,"=");split(b[2],c,".");s=$4-1;print $1"\t"s"\t"$5"\t"c[1]"\t"".""\t"$7}' | sort -k1,1 -k2,2n | uniq >$mRNA_bed

}

get_top500_peak_and_submit500bp() {
    # 获得最显著的top500peak及其顶点周围500bp区域
    $code_dir/s8_sort_peak_and_submit.r "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output_jiang/3D7-T3_peaks.xls" "/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
    head -n 500 $macs2_all_peaks_one_submit_sorted | sort -k1,1 -k2,2n | awk '{print $1"\t"$2-1"\t"$3}' >$macs2_top500_peak
    head -n 500 $macs2_all_peaks_one_submit_sorted | awk '{print $1"\t"$4-1-250"\t"$4+250}' | sort -k1,1 -k2,2n | $bedtools intersect -a stdin -b $macs2_top500_peak >$macs2_top500_peak_submit500

}

get_6mA_closest_genes_nomatter_distance() {
    # 获得peak最近的基因，不限制距离
    $bedtools closest -a $macs2_all_peaks -b $gene_bed -D a >$closet_gene_info
}

get_6mA_closest_genes_peak_all_in_exons() {
    # 要求整个peak都在基因的exon里
    closet_gene_dir=$m6mA_closest_gene_global_dir/all_in_exons
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_all_peaks -b $exon_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_peak_ATG_TAA_2kb() {
    ATG_TAA_2kb_bed=$(dirname $mRNA_bed)"/ATG_TAA_2kb.bed"
    $bedtools slop -i $mRNA_bed -g $genome_size -l 2000 -r 2000 -s >$ATG_TAA_2kb_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_TAA_2kb
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_all_peaks -b $ATG_TAA_2kb_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_peak_ATG_TAA_2kb_just_flank() {
    ATG_TAA_2kb_just_flank_bed=$(dirname $mRNA_bed)"/ATG_TAA_2kb_just_flank.bed"
    $bedtools flank -i $mRNA_bed -g $genome_size -l 2000 -r 2000 -s | $bedtools subtract -a stdin -b $exon_bed >$ATG_TAA_2kb_just_flank_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_TAA_2kb_just_flank
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_all_peaks -b $ATG_TAA_2kb_just_flank_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_peak_ATG_2kb() {
    # ATG_2kb是gene body+ATG前2kb
    ATG_2kb_bed=$(dirname $mRNA_bed)"/ATG_2kb.bed"
    $bedtools slop -i $mRNA_bed -g $genomes_ize -l 2000 -r 0 -s >$ATG_2kb_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_2kb
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_all_peaks -b $ATG_2kb_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_peak_ATG_2kb_just_flank() {
    # 只是ATG前2kb
    ATG_2kb_just_flank_bed=$(dirname $mRNA_bed)"/ATG_2kb_just_flank.bed"
    $bedtools flank -i $mRNA_bed -g $genome_size -l 2000 -r 0 -s | $bedtools subtract -a stdin -b $exon_bed >$ATG_2kb_just_flank_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_2kb_just_flank
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_all_peaks -b $ATG_2kb_just_flank_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_top500_peak_all_in_exons() {
    closet_gene_dir=$m6mA_closest_gene_global_dir/all_in_exons_top500peak
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_top500_peak -b $exon_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_top500_peak_ATG_2kb_just_flank() {
    ATG_2kb_just_flank_bed=$(dirname $mRNA_bed)"/ATG_2kb_just_flank.bed"
    # $bedtools flank -i $mRNA_bed -g $genome_size -l 2000 -r 0 -s | $bedtools subtract -a stdin -b $exon_bed >$ATG_2kb_just_flank_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_2kb_just_flank_top500peak
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_top500_peak -b $ATG_2kb_just_flank_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_top500_submit500_peak_all_in_exons() {
    closet_gene_dir=$m6mA_closest_gene_global_dir/all_in_exons_top500peak_submit500
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_top500_peak_submit500 -b $exon_bed -f 1 -wo >$closet_gene_info
}

get_6mA_closest_genes_top500_submit500_peak_ATG_2kb_just_flank() {
    ATG_2kb_just_flank_bed=$(dirname $mRNA_bed)"/ATG_2kb_just_flank.bed"
    # $bedtools flank -i $mRNA_bed -g $genome_size -l 2000 -r 0 -s | $bedtools subtract -a stdin -b $exon_bed >$ATG_2kb_just_flank_bed
    closet_gene_dir=$m6mA_closest_gene_global_dir/ATG_2kb_just_flank_top500peak_submit500
    mkdir -p $closet_gene_dir
    closet_gene_info=$closet_gene_dir/6mA_closet_genes_info
    $bedtools intersect -a $macs2_top500_peak_submit500 -b $ATG_2kb_just_flank_bed -f 1 -wo >$closet_gene_info
}

get_interested_genes_genome_index() {
    gene_list=$1
    grep -f $gene_list $gene_bed >$(echo $gene_list | sed 's/_list/.bed/g')
}

get_overlap_items_count() {
    interested_gene=$1
    closet_gene_info=$2
    overlap_items=$(echo $interested_gene | sed 's/genes_list/overlap_all_peak_closet_genes/g')

    awk '{print $7}' $closet_gene_info | sort | uniq >file2.inter
    sort $interested_gene file2.inter | uniq -d >$overlap_items
    overlap_items_count=$(wc -l $overlap_items | awk '{print $1}')
    echo $overlap_items_count

    rm -rf file2.inter
}

f() {
    for pseudo in $(echo "include_pseudo exclude_pseudo" | tr " " "\n"); do
        for measure in $(echo $peak_gene_region_combination | tr " " "\n"); do
            closet_gene_info=$m6mA_closest_gene_global_dir/$measure/6mA_closet_genes_info
            peak_related_gene_count=$(awk '{print $7}' $closet_gene_info | sort | uniq | wc -l)
            let whole_genome_gene_count_except_peak_related=whole_genome_gene_count-peak_related_gene_count

            for gene_group in $(echo $gene_group_type | tr " " "\n"); do
                pseudo_gene_global_dir=$specific_gene_global_dir/$pseudo
                gene_group_count=$(wc -l $pseudo_gene_global_dir/${gene_group}/${gene_group}_genes_list | awk '{print $1}')
                get_interested_genes_genome_index $pseudo_gene_global_dir/${gene_group}/${gene_group}_genes_list
                overlap_count=$(get_overlap_items_count $pseudo_gene_global_dir/${gene_group}/${gene_group}_genes_list $closet_gene_info)
                p="$($code_dir/s8_get_enrichment_p_value.r $overlap_count $peak_related_gene_count $whole_genome_gene_count_except_peak_related $gene_group_count)"
                echo $pseudo" "$measure" "$gene_group" "$overlap_count" "$peak_related_gene_count" "$whole_genome_gene_count_except_peak_related" "$gene_group_count" "$p
            done
        done
    done

}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
anno_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno"
specific_gene_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups"

anno_gtf=$anno_dir/"PlasmoDB-36_Pfalciparum3D7.gff"
genome_size=$anno_dir/../genome/"chrom.length"
gene_bed=$anno_dir/"genes.bed"
exon_bed=$anno_dir/"exons.bed"
mRNA_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/mRNA.bed"
macs2_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"

m6mA_closest_gene_global_dir="$macs2_output_dir/6mA_closet_genes"

# macs2_all_peaks=$macs2_output_dir/"peak.bed"
macs2_all_peaks="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"
macs2_all_peaks_one_submit_sorted="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/sorted_peak.txt"
macs2_top500_peak="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/top500_peak.bed"
macs2_top500_peak_submit500="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/top500_peak_submit500.bed"

significant_enriched_gene_sets_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/significant_enriched_gene_sets"

peak_gene_region_combination="all_in_exons all_in_exons_top500peak all_in_exons_top500peak_submit500 ATG_2kb_just_flank ATG_2kb_just_flank_top500peak ATG_2kb_just_flank_top500peak_submit500 "
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
gene_group_type="BER HDR MMR  NER invasion rhoptry RIF STEVOR SURF VAR PHIST MC-2TM DnaJ DR"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# mkdir -p $closet_gene_dir
get_whole_genome_gene_and_mRNA

get_6mA_cloest_genes_peak_all_in_exons
get_6mA_closest_genes_peak_ATG_TAA_2kb
get_6mA_closest_genes_peak_ATG_2kb
get_6mA_closest_genes_peak_ATG_TAA_2kb_just_flank
get_6mA_closest_genes_peak_ATG_2kb_just_flank
get_6mA_closest_genes_top500_peak_all_in_exons
get_6mA_closest_genes_top500_peak_ATG_2kb_just_flank
get_6mA_closest_genes_top500_submit500_peak_all_in_exons
get_6mA_closest_genes_top500_submit500_peak_ATG_2kb_just_flank

whole_genome_gene_count=$(wc -l $gene_bed | awk '{print $1}')
f >$significant_enriched_gene_sets_dir/enriched_gene_sets_info.txt
$code_dir/s8_plot_gene_sets_enriched.r $significant_enriched_gene_sets_dir/enriched_gene_sets_info.txt
