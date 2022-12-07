#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
locate_motif_in_peak_and_control_two_outgroup_consistent_remove_overlap_sites_between_different_motifs() {
    mkdir -p $motif_bed_dir
    cd $motif_bed_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_peak_seq_dir/${autosome}.peak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6}' >${autosome}.motif.peak.bed
    done

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_nopeak_seq_dir/${autosome}.nopeak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6}' >${autosome}.motif.nopeak.bed
    done
    cd -
}

find_motif_each_site_index_and_remove_overlap_site_remove_overlap_sites_between_different_motifs() {
    cd $motif_bed_dir

    for region_type in $(echo "peak nopeak" | tr " " "\n"); do
        rm -rf $(find ../site* | grep .${region_type}.)
        cat *.${region_type}.* | awk -v region=${region_type} '$4=="+"{
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site1/"region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site2/"region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site3/"region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site4/"region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site5/"region".bed";
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site6/"region".bed";

                }$4=="-"{
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site1/"region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site2/"region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site3/"region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site4/"region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site5/"region".bed";
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site6/"region".bed";

                }'

        cat $(find ../site* | grep ${region_type}) | awk '{print $1"\t"$2"\t"$3}' | sort -k 1,1 -k2,2n | uniq -d >${region_type}.allSiteDupli.bed

        rm -rf $(find ../site*/*.${region_type}.plus* 2>/dev/null)
        rm -rf $(find ../site*/*.${region_type}.minus* 2>/dev/null)
        for site_index in $(echo $sites | tr " " "\n"); do
            $bedtools subtract -a ../site${site_index}/${region_type}.bed -b ${region_type}.allSiteDupli.bed |
                awk -v site="../site${site_index}" -v region=${region_type} '$4=="+"{
                print $1"\t"$2"\t"$3 >> site"/"$1"."region".plus.bed"
            }$4=="-"{
                print $1"\t"$2"\t"$3 >> site"/"$1"."region".minus.bed"
            }'
        done
    done
}

locate_motif_in_peak_and_control_two_outgroup_consistent_only_retain_nonoverlap_motifs() {
    mkdir -p $motif_bed_dir
    cd $motif_bed_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cmd="cat $core_peak_seq_dir/${autosome}.peak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed|sort -k1,1 -k2,2n"

        $code_dir/s7_retain_nonoverlap_motifs.r "$cmd" | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6}' >${autosome}.motif.peak.bed
    done

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cmd="cat $core_nopeak_seq_dir/${autosome}.nopeak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed|sort -k1,1 -k2,2n"
        $code_dir/s7_retain_nonoverlap_motifs.r "$cmd" | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6}' >${autosome}.motif.nopeak.bed
    done
    cd -
}

find_motif_each_site_index_and_remove_overlap_site_only_retain_nonoverlap_motifs() {
    cd $motif_bed_dir

    for region_type in $(echo "peak nopeak" | tr " " "\n"); do
        rm -rf $(find ../site* 2>/dev/null | grep .${region_type}.)
        cat *.${region_type}.* | awk -v region=${region_type} '$4=="+"{
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site1/"region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site2/"region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site3/"region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site4/"region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site5/"region".bed";
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site6/"region".bed";

                }$4=="-"{
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site1/"region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site2/"region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site3/"region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site4/"region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site5/"region".bed";
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site6/"region".bed";

                }'

        cat $(find ../site* | grep ${region_type}) | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n | uniq -d >${region_type}.allSiteDupli.bed

        rm -rf $(find ../site*/*.${region_type}.plus* 2>/dev/null)
        rm -rf $(find ../site*/*.${region_type}.minus* 2>/dev/null)
        for site_index in $(echo $sites | tr " " "\n"); do
            $bedtools subtract -a ../site${site_index}/${region_type}.bed -b ${region_type}.allSiteDupli.bed |
                awk -v site="../site${site_index}" -v region=${region_type} '$4=="+"{
                print $1"\t"$2"\t"$3 >> site"/"$1"."region".plus.bed"
            }$4=="-"{
                print $1"\t"$2"\t"$3 >> site"/"$1"."region".minus.bed"
            }'
        done
    done
}

# VARIABLE NAMING TODO:
# motif=motif_template
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
motif_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs"
motif_seq=motif_seq_template
motif=$motif_dir/$motif_seq

macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
single_motif_output_dir=$single_motif_global_output_dir/$motif_seq
motif_bed_dir=$single_motif_output_dir/motif_bed
sites="1 2 3 4 5 6"

core_peak_seq_dir=$macs2_two_outgroup_consistent_dir/peak_seq
core_nopeak_seq_dir=$macs2_two_outgroup_consistent_dir/nopeak_seq

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
mkdir -p $single_motif_output_dir $single_motif_output_dir/site1 $single_motif_output_dir/site2 $single_motif_output_dir/site3 $single_motif_output_dir/site4 $single_motif_output_dir/site5 $single_motif_output_dir/site6
locate_motif_in_peak_and_control_two_outgroup_consistent_only_retain_nonoverlap_motifs
find_motif_each_site_index_and_remove_overlap_site_only_retain_nonoverlap_motifs
