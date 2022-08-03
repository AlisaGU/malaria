#!/bin/bash
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
query_abbre=query_abbre_template
query_id=query_id_template

export global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species"
export working_dir=$global_dir/alignment

query_full_name=$query_abbre
export last_train="/picb/evolgen/users/gushanshan/software/LAST/bin/last-train"
export last_split="/picb/evolgen/users/gushanshan/software/LAST/bin/last-split"
export last_postmask="/picb/evolgen/users/gushanshan/software/LAST/bin/last-postmask"
export lastal="/picb/evolgen/users/gushanshan/software/LAST/bin/lastal"
export reference_abbre="P_falciparum"
export reference_full_name="P_falciparum"

query_genome_sequence=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/genomes_and_anno/${query_abbre}/genome/PlasmoDB-58_${query_id}_Genome.fasta

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex

mkdir -p $working_dir/$query_abbre
cd $working_dir/$query_abbre
${last_train} -P30 --revsym -D1e9 --sample-number=5000 $working_dir/${reference_abbre} ${query_genome_sequence} >my.train
$lastal -P30 -D1e9 -m100 -p my.train $working_dir/${reference_abbre} ${query_genome_sequence} | $last_split -fMAF+ >many_to_one_maf

$last_split -r many_to_one_maf | $last_postmask >one_to_one_maf
/picb/evolgen/users/gushanshan/projects/Andrias/code/1_create_genome_alignment/LAST/maf.rename.species.S.pl one_to_one_maf $reference_full_name $query_full_name ${reference_abbre}.${query_abbre}.maf >one_to_one_maf.stat
