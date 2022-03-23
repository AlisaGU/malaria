#!/bin/bash
#SBATCH -n 2
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
dedup() {
    input=$1
    output=$(echo $input | sed 's/.bam/.sort.bam/g')
    sample=$(basename $input | sed 's/.bam//g')

    $sambamba markdup -t 8 --tmpdir=$bam_dir/temp $input $output 1>$bam_dir/${sample}".dedup.log" 2>&1
}

get_effective_genome_size() {
    $samtools depth -o $input_sort_bam_depth $input_sort_bam
    effective_genome_size=$(wc -l $input_sort_bam_depth)
    rm -rf $input_sort_bam_depth
    echo $effective_genome_size
}

call_peak() {
    $macs2 callpeak -t $treatment_sort_bam -c $input_sort_bam -f BAMPE -g 23008401 -n 3D7_2rd --outdir $bam_dir/macs2_output --bdg --trackline 2>$bam_dir/macs2_output/macs2.log
}
export -f dedup
# VARIABLE NAMING TODO:
export bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd"

treatment_sort_bam=$bam_dir/"3D7-T3_ChIP.sort.bam"
input_sort_bam=$bam_dir"/3D7-T3_Input.sort.bam"
input_sort_bam_depth=$(echo $input_sort_bam | sed 's/bam/depth/g')
export sambamba="/home/gushanshan/anaconda3/envs/raw_read_to_snp_calling/bin/sambamba"
export parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
samtools="/picb/evolgen/users/gushanshan/GenomeAnnotation/samtools/samtools-1.10/samtools_install/bin/samtools"
macs2="/home/gushanshan/anaconda3/envs/vscode_r/bin/macs2"
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
all_input=$(ls $bam_dir/*.bam)

# echo $all_input|tr " " "\n"|$parallel -j 2 dedup

# effective_genome_size=`get_effective_genome_size` # pf 3D7 effective genome size = 23008401 bp

call_peak

$code_dir/s2_1_peak_process.sh
