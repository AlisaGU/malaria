#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
find_flank_region() {
    mkdir -p flank_bed
    for peak_type in $(echo "peak nopeak" | tr " " "\n"); do
        for chrom in $(echo $autosomes | tr " " "\n"); do
            cat motif_bed/${chrom}.motif.${peak_type}.bed | grep -E -v "GAMSAA|GAHGAW" |
                sort -k1,1 -k2,2n | $bedtools merge -i stdin >motif_bed/${chrom}.motif.${peak_type}.intermediate.bed

            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                # 去掉flank region与motif区域的重叠,要求位点在peak里
                awk -v flank_site=${flank_site} '{print $1"\t"$2-flank_site"\t"$2-flank_site+1"\n"$1"\t"$3+flank_site-1"\t"$3+flank_site}' motif_bed/${chrom}.motif.${peak_type}.intermediate.bed |
                    $bedtools subtract -a stdin -b motif_bed/${chrom}.motif.${peak_type}.intermediate.bed |
                    $bedtools intersect -a stdin -b $macs2_two_outgroup_consistent_dir/${chrom}.${peak_type}.bed >flank_bed/${chrom}.motif.${peak_type}.flank${flank_site}.intermediate.bed
            done

            # 获得flank1-6中重复出现的位点
            cat $(find flank_bed -name ${chrom}.motif.${peak_type}.flank[1-6].intermediate.bed) |
                sort -k1,1 -k2,2n |
                uniq -d >flank_bed/dupli_among_1_to_6

            # 在flank1-6中重复出现的位点要删掉
            for flank_site in $(seq 1 6); do
                $bedtools subtract -a flank_bed/${chrom}.motif.${peak_type}.flank${flank_site}.intermediate.bed -b flank_bed/dupli_among_1_to_6 >flank_bed/${chrom}.motif.${peak_type}.flank${flank_site}.bed
            done

            # 对于flank7-10，删掉在更近距离中出现的位点
            for flank_site in $(echo "7 8 9 10" | tr " " "\n"); do
                let before_count=flank_site-1
                cat $(find flank_bed -name ${chrom}.motif.${peak_type}.flank[1-${before_count}].intermediate.bed) |
                    sort -k1,1 -k2,2n |
                    uniq >flank_bed/site_among_1_to_before

                $bedtools subtract -a flank_bed/${chrom}.motif.${peak_type}.flank${flank_site}.intermediate.bed -b flank_bed/site_among_1_to_before >flank_bed/${chrom}.motif.${peak_type}.flank${flank_site}.bed
            done
            rm -rf $(find ./ -name *intermediate*) flank_bed/dupli_among_1_to_6 flank_bed/site_among_1_to_6 flank_bed/site_among_1_to_before
        done
    done
    cd flank_bed
    for flank_site in $(seq 1 10); do
        mkdir -p flank_${flank_site}
        mv *.flank${flank_site}.bed flank_${flank_site}/
        cd flank_${flank_site}
        for filename in $(/usr/bin/ls); do
            new_filename=$(echo $filename | sed -e "s/.motif//g" -e "s/.flank${flank_site}//g")
            mv $filename $new_filename
        done
        cd -
    done
}

find_nonoverlap_motif() {
    for peak_type in $(echo "peak nopeak" | tr " " "\n"); do
        for chrom in $(echo $autosomes | tr " " "\n"); do
            # 非简并式
            cmd="grep -E -v \"GAMSAA|GAHGAW\" motif_bed/${chrom}.motif.${peak_type}.bed | sort -k1,1 -k2,2n | $bedtools merge -i stdin | awk '{print \$0\"\t\"NR}' | $bedtools intersect -a motif_bed/${chrom}.motif.${peak_type}.bed -b stdin -wb | grep -E -v \"GAMSAA|GAHGAW\""
            $code_dir/s7_seperate_continuous_motif_nondegenerate.R "$(pwd)" "$cmd" "motif_bed/${chrom}.motif.${peak_type}.nonDegenerate.nonoverlap.bed"
            # 简并式
            grep -E "GAMSAA|GAHGAW" motif_bed/${chrom}.motif.${peak_type}.bed |
                awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq -d |
                awk '{print $0"\t""common"}' >motif_bed/${chrom}.motif.${peak_type}.degenerate.completely_overlap.bed

            grep -E "GAMSAA|GAHGAW" motif_bed/${chrom}.motif.${peak_type}.bed |
                $bedtools subtract -f 1 -a stdin -b motif_bed/${chrom}.motif.${peak_type}.degenerate.completely_overlap.bed >motif_bed/${chrom}.motif.${peak_type}.degenerate.unique.bed

            cat motif_bed/${chrom}.motif.${peak_type}.degenerate.completely_overlap.bed motif_bed/${chrom}.motif.${peak_type}.degenerate.unique.bed |
                sort -k1,1 -k2,2n >motif_bed/${chrom}.motif.${peak_type}.overlap_unique.bed

            cmd="$bedtools merge -i motif_bed/${chrom}.motif.${peak_type}.overlap_unique.bed | awk '{print \$0\"\t\"NR}' | $bedtools intersect -a motif_bed/${chrom}.motif.${peak_type}.overlap_unique.bed -b stdin -wb"
            $code_dir/s7_seperate_continuous_motif_degenerate.R "$(pwd)" "$cmd" "motif_bed/${chrom}.motif.${peak_type}.degenerate.nonoverlap.bed"
        done
    done
}

find_specific_site() {
    for peak_type in $(echo "peak nopeak" | tr " " "\n"); do
        for motif in $(echo $motifs | tr " " "\n"); do
            rm -rf $(find ../${motif}/ | grep ${peak_type}.bed)
            for chrom in $(echo $autosomes | tr " " "\n"); do

                mkdir -p ../${motif}/site1 ../${motif}/site2 ../${motif}/site3 ../${motif}/site4 ../${motif}/site5 ../${motif}/site6
                if [ "$motif" = "GAMSAA" -o "$motif" = "GAHGAW"]; then
                    grep -E "${motif}|common" ${chrom}.motif.${peak_type}.degenerate.nonoverlap.bed | awk -v motif=${motif} -v region=${peak_type} '$4=="+"{
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../"motif"/site1/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../"motif"/site2/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../"motif"/site3/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../"motif"/site4/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../"motif"/site5/"$1"."region".bed";
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../"motif"/site6/"$1"."region".bed";

                }$4=="-"{
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../"motif"/site1/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../"motif"/site2/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../"motif"/site3/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../"motif"/site4/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../"motif"/site5/"$1"."region".bed";
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../"motif"/site6/"$1"."region".bed";
                }'
                else
                    grep ${motif} ${chrom}.motif.${peak_type}.nonDegenerate.nonoverlap.bed | awk -v motif=${motif} -v region=${peak_type} '$4=="+"{
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../"motif"/site1/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../"motif"/site2/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../"motif"/site3/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../"motif"/site4/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../"motif"/site5/"$1"."region".bed";
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../"motif"/site6/"$1"."region".bed";

                }$4=="-"{
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../"motif"/site1/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../"motif"/site2/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../"motif"/site3/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../"motif"/site4/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../"motif"/site5/"$1"."region".bed";
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../"motif"/site6/"$1"."region".bed";
                }'
                fi
            done
        done
    done
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
genome_len="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
flank_sites=$(seq 1 10)
motifs="GAMSAA GAHGAW GAAGAA GAACAA GACGAA GACCAA GAAGAT GATGAA GATGAT GACGAT"
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -e

# 1. 计算简并式motif在基因组的位置 # TODO: CHECK DONE
sed "s/motif_seq_template/all_motifs/g" $code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh >$code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh_all_motifs
sbatch $code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh_all_motifs

# 2. 合并所有motif，寻找flank区间 # TODO: CHECK DONE
cd $single_motif_global_output_dir/all_motifs
find_flank_region

# 3.在所有motif水平上寻找nonoverlap的motif，先定两头 # TODO: CHECK DONE
cd $single_motif_global_output_dir/all_motifs
find_nonoverlap_motif

# 4. 为每种motif定义1，2，3，4，5，6号位 # TODO: CHECK DONE
cd $single_motif_global_output_dir/all_motifs/motif_bed
find_specific_site

# 5. 创建和之前相同的文件夹模式
cd $single_motif_global_output_dir/all_motifs
for motif in $(echo "GAMSAA GAHGAW GAAGAA GAACAA GACGAA GACCAA GAAGAT GATGAA GATGAT GACGAT" | tr " " "\n"); do
    cd $motif
    ln -s $single_motif_global_output_dir/all_motifs/flank_bed flank_nomatter_strand.remove_close_overlap
    cd -
done

# 6. 计算motif位点的突变速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    for motif_seq in $(echo $motifs | tr " " "\n"); do
        sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh | sed "s/motif_seq_template/${motif_seq}/g" >$code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh_${popu_symbol}_${motif_seq}
        sbatch $code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh_${popu_symbol}_${motif_seq}
    done
done

# 7. 计算flank位点的突变速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    for motif_seq in $(echo $motifs | tr " " "\n"); do
        sed -e "s/motif_seq_template/$motif_seq/g" -e "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_mutation_load_for_each_site_of_flank_asWhole.sh >$code_dir/s7_compute_mutation_load_for_each_site_of_flank_asWhole.sh_${popu_symbol}_$motif_seq
        sbatch $code_dir/s7_compute_mutation_load_for_each_site_of_flank_asWhole.sh_${popu_symbol}_$motif_seq
    done
done

# 8. 计算A、T、C、G四种碱基在信号区的基底速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh >$code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh_$popu_symbol
    sbatch $code_dir/s7_compute_each_chrom_peak_base_mutation_density.sh_$popu_symbol
done

# 9. 计算各位点碱基比例
for motif_seq in $(echo $motifs | tr " " "\n"); do
    sed "s/motif_template/$motif_seq/g" $code_dir/s7_compute_each_site_base_propor_asWhole.sh >$code_dir/s7_compute_each_site_base_propor_asWhole.sh_$motif_seq
    sbatch $code_dir/s7_compute_each_site_base_propor_asWhole.sh_$motif_seq
done
