#!/bin/bash
#SBATCH -n 12
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_peak_proportion() {
    flank_length=$1

    cd $whole_genome_bed_dir

    peak_control_length_in_gene=peak_control_length_in_flank${flank_length}_wholePeak
    rm -rf $peak_control_length_in_gene

    cat gene_flank${flank_length}.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc -wb | wc -l)
        value=""

        gene_name=$(echo $line | awk '{print $4}')
        gene_flank_length=$(echo $line | awk '{s+=$3-$2}END{print s}')

        gene_flank_length_in_peak=0
        if [ $peak_exist_in_gene_status -gt 0 ]; then
            gene_flank_length_in_peak=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | awk '{s+=$3-$2}END{print s}')
        fi
        let gene_flank_length_in_control=gene_flank_length-gene_flank_length_in_peak
        value+="$(echo $gene_name $gene_flank_length_in_peak $gene_flank_length_in_control)"

        echo $value >>$peak_control_length_in_gene
    done
}

export -f get_peak_proportion
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"

export whole_genome_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/whole_Genome"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export mutation_density_dir=$global_output_dir/${popu_symbol}/wholeGenome_peak25Perc/3D7
export variant_dir=$global_output_dir/$popu_symbol/variant_3D7 ## 变异信息是以3D7为祖先型的
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"

export gene_bed=$whole_genome_bed_dir/genes.bed
export gene_flank2kb_bed=$whole_genome_bed_dir/gene_flank2kb.bed
export gene_flank1kb_bed=$whole_genome_bed_dir/gene_flank1kb.bed

export peak_loc=$m6mA_bam_dir/peak.bed

export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex

$parallel -j 3 get_peak_proportion ::: "0kb" "1kb" "2kb"
