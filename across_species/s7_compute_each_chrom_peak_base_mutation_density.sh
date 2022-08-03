#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $10+$14+$15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $12+$20+$21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $13+$18+$19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $11+$16+$17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $26+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $27+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $28+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${chrom}.species3_polysite.bed -b $region_bed | awk '{print $29+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        :
    else
        chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

compute_specific_base_peak_base_mutation_density() {
    specific_base=$1
    var_type=$2
    mean_variant_count=$mutation_density_dir/${specific_base}_${var_type}_mean_variant_count
    rm -rf $mean_variant_count
    for chrom in $(echo $autosomes | tr " " "\n"); do
        region_bed=$two_outgroup_consistent_dir/${chrom}.two_consistent_region.peak.bed
        value=""
        comple_base=$(get_complementary_base $specific_base)
        for base in $(echo $comple_base | tr "_" "\n"); do
            base_region_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')
            base_var_count=$(get_average_variant_count_for_chrom_base_count $chrom $region_bed "${base}_${var_type}")
            value+=$(echo ${base_var_count} ${base_region_len})
            value+=" "
        done
        echo $value >>$mean_variant_count
    done
}

get_complementary_base() {
    base=$1
    if [ $base == "G" ]; then
        echo "G_C"
    elif [ $base == "A" ]; then
        echo "A_T"
    elif [ $base == "T" ]; then
        echo "T_A"
    elif [ $base == "C" ]; then
        echo "C_G"
    elif [ $base == "A_T_C" ]; then
        echo "A_T_C_G"
    elif [ $base == "A_T" ]; then
        echo "A_T"
    elif [ $base == "A_C" ]; then
        echo "A_T_C_G"
    elif [ $base == "G_C" ]; then
        echo "G_C"
    fi
}

# VARIABLE NAMING TODO:
variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf"
mutation_density_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/peak_base_mutation_density"
two_outgroup_consistent_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi/peak_nopeak"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $mutation_density_dir
for base_type in $(echo "A T C G" | tr " " "\n"); do
    for vartype in $(echo "snps single_del" | tr " " "\n"); do
        compute_specific_base_peak_base_mutation_density $base_type $vartype
    done
done
