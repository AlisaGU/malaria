#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf
grep -v "P_falciparum" P_falciparum.P_billcollinsi.P_reichenowi.first_seq_in_block | awk '{print $1" "$2" "$3" "$4" "$5" "$6}' | while read line; do
    line_number=$(grep -n "$line" P_falciparum.P_billcollinsi.P_reichenowi.maf | awk -F":" '{print $1}')
    let lower=line_number-10
    let upper=line_number+10

    last_block_end=$(awk -v line_number=$line_number -v lower=$lower 'NR>=lower && NR<=line_number {print $1" "$2" "$3" "$4" "$5" "$6}' P_falciparum.P_billcollinsi.P_reichenowi.maf | grep "s P_falciparum" | tail -n 1 | awk '{print $3+$4}')
    first_block_start=$(awk -v line_number=$line_number -v upper=$upper 'NR>=line_number && NR<=upper {print $1" "$2" "$3" "$4" "$5" "$6}' P_falciparum.P_billcollinsi.P_reichenowi.maf | grep "s P_falciparum" | head -n 1 | awk '{print $3}')
    if [ $last_block_end -eq $first_block_start ]; then
        echo "$line:    continuous"
    else
        echo "$line:    NO"
    fi
done

/picb/evolgen/users/gushanshan/software/ucsc_pairwiseGenomeAlignmrent/data/bin/mafFilter Pf3D7_02_v3.P_falciparum.P_billcollinsi.P_reichenowi.maf -minRow=2 -maxRow=3 -reject=Pf3D7_02_v3.reject.RMsingle_seq_block.txt >Pf3D7_02_v3.P_falciparum.P_billcollinsi.P_reichenowi.RMsingle_seq_block.maf

/picb/evolgen/users/gushanshan/software/ucsc_pairwiseGenomeAlignmrent/data/bin/mafFilter Pf3D7_02_v3.P_falciparum.P_billcollinsi.P_reichenowi.maf -needComp="P_falciparum" -minRow=1 -reject=Pf3D7_02_v3.reject.noP_falciparum.txt >Pf3D7_02_v3.P_falciparum.P_billcollinsi.P_reichenowi.includingP_falciparum.maf
