#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_each_site_base_propor() {
    cd $motif_each_dir
    for site in $(seq 1 6); do
        cd "site"$site
        for peak_type in $(echo "peak"); do
            each_base_appear_count=site${site}.${peak_type}.base_count
            rm -rf $each_base_appear_count
            for chrom in $(echo $autosomes | tr " " "\n"); do
                this_site_genome_pos=${chrom}.${peak_type}.intermediate
                rm -rf $this_site_genome_pos
                sort -k1,1 -k2,2n ${chrom}.${peak_type}* >$this_site_genome_pos
                $seqkit subseq --quiet --bed $this_site_genome_pos $genome | $seqkit fx2tab -i -H -n -C "A" -C "T" -C "C" -C "G" | grep -v "^#" | awk '{a+=$2;t+=$3;c+=$4;g+=$5} END{print a"\t"t"\t"c"\t"g"\t"a+t+c+g}' >>$each_base_appear_count
                rm -rf $this_site_genome_pos
            done
        done
        cd ..
    done

    cd flank_nomatter_strand.remove_close_overlap
    for flank in $(seq 1 10); do
        cd "flank_"${flank}
        for peak_type in $(echo "peak"); do
            each_base_appear_count=flank${flank}.${peak_type}.base_count
            rm -rf $each_base_appear_count
            for chrom in $(echo $autosomes | tr " " "\n"); do
                $seqkit subseq --quiet --bed ${chrom}.${peak_type}.bed $genome | $seqkit fx2tab -i -H -n -C "A" -C "T" -C "C" -C "G" | grep -v "^#" | awk '{a+=$2;t+=$3;c+=$4;g+=$5} END{print a"\t"t"\t"c"\t"g"\t"a+t+c+g}' >>$each_base_appear_count
            done
        done
        cd ..
    done
}
# VARIABLE NAMING TODO:
motif=motif_template
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent

motif_each_dir=""
# if [ "$motif_seq" = "GAWGAW" ]; then
if [ "$motif" = "GAWGAW" ]; then
    motif_each_dir=${macs2_two_outgroup_consistent_dir}/single_motif_pattern/${motif}
else
    motif_each_dir=${macs2_two_outgroup_consistent_dir}/single_motif_pattern/all_motifs/${motif}
fi

genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
compute_each_site_base_propor
