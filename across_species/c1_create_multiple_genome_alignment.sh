#!/bin/bash
#SBATCH -n 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
pairwise_alignment() {
    query_abbre=$1
    query_id=$2
    query_full_name=$1
    query_genome_sequence=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/genomes_and_anno/${query_abbre}/genome/PlasmoDB-58_${query_id}_Genome.fasta
    mkdir -p $working_dir/$query_abbre
    cd $working_dir/$query_abbre
    ${last_train} -P30 --revsym -D1e9 --sample-number=5000 $working_dir/${reference_abbre} ${query_genome_sequence} >my.train
    $lastal -P30 -D1e9 -m100 -p my.train $working_dir/${reference_abbre} ${query_genome_sequence} | $last_split -fMAF+ >many_to_one_maf

    $last_split -r many_to_one_maf | $last_postmask >one_to_one_maf
    /picb/evolgen/users/gushanshan/projects/Andrias/code/1_create_genome_alignment/LAST/maf.rename.species.S.pl one_to_one_maf $reference_full_name $query_full_name ${reference_abbre}.${query_abbre}.maf >one_to_one_maf.stat

}

pairwise_alignment_sbatch() {
    query_abbres="P_billcollinsi P_blacklocki P_adleri P_gaboni"
    query_ids="PbillcollinsiG01 PblacklockiG01 PadleriG01 PgaboniG01"

    for order in $(seq 1 4); do
        query_abbre=$(echo $query_abbres | tr " " "\n" | awk -v order=$order 'NR==order')
        query_id=$(echo $query_ids | tr " " "\n" | awk -v order=$order 'NR==order')
        sed -e "s/query_abbre_template/$query_abbre/g" -e "s/query_id_template/${query_id}/g" $code_dir/s1_LAST.sh >$code_dir/s1_LAST.sh_${query_abbre}
        sbatch $code_dir/s1_LAST.sh_${query_abbre}
    done
}

summary() {
    maffile=$1
    ref_len_in_maf=$(grep "^s P_falciparum" $maffile | awk '{print $4}' | awk '{s+=$1}END{print s}')
    ref_whole_genome_len="23332839"
    ratio=$(echo "scale=10;${ref_len_in_maf}/${ref_whole_genome_len}" | bc)
    ratio1=$(printf "%0.6f" $ratio)
    echo "P_falciparum region in maf: $ref_len_in_maf, accounting for $ratio1 in P_falciparum whole genome"
}

split_by_chromosome() {
    ChrName=$1
    chrAlignStart=$(grep -n "s P_falciparum.${ChrName} " $multiple_alignment_file | head -1 | cut -d ":" -f 1)
    chrAlignStart=$(expr $chrAlignStart - 1)
    chrAlignEnd=$(grep -n "s P_falciparum.${ChrName} " $multiple_alignment_file | tail -1 | cut -d ":" -f 1)
    chrAlignEnd=$(expr $chrAlignEnd + 2)
    sed -n "${chrAlignStart},${chrAlignEnd}p" $multiple_alignment_file >./${ChrName}.$(basename $multiple_alignment_file)
}
export -f pairwise_alignment split_by_chromosome
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/across_species"
export global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species"
export working_dir=$global_dir/alignment
export reference_genome=$global_dir/genomes_and_anno/P_falciparum/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta
export reference_abbre="P_falciparum"
export reference_full_name="P_falciparum"

export lastdb="/picb/evolgen/users/gushanshan/software/LAST/bin/lastdb"
export last_train="/picb/evolgen/users/gushanshan/software/LAST/bin/last-train"
export last_split="/picb/evolgen/users/gushanshan/software/LAST/bin/last-split"
export last_postmask="/picb/evolgen/users/gushanshan/software/LAST/bin/last-postmask"
export lastal="/picb/evolgen/users/gushanshan/software/LAST/bin/lastal"
# export multiple_alignment_file=$global_dir/maf2vcf/P_falciparum.P_reichenowi.P_praefalciparum.maf
export multiple_alignment_file=$global_dir/maf2vcf/P_falciparum.P_billcollinsi.P_reichenowi.maf

multiz="/home/gushanshan/anaconda3/bin/multiz"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
# mkdir -p $working_dir
# cd $working_dir

# $lastdb -P20 -uMAM8 ${reference_abbre} ${reference_genome}

pairwise_alignment_sbatch
# $parallel -j 2 --xapply pairwise_alignment ::: "P_praefalciparum" "P_reichenowi" ::: "PpraefalciparumG01" "PreichenowiG01"

# $multiz P_praefalciparum/P_falciparum.P_praefalciparum.maf P_reichenowi/P_falciparum.P_reichenowi.maf 1 >$multiple_alignment_file
$multiz P_reichenowi/P_falciparum.P_reichenowi.maf P_billcollinsi/P_falciparum.P_billcollinsi.maf 1 >$multiple_alignment_file

cd $global_dir/maf2vcf
$parallel -j 4 split_by_chromosome ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3
