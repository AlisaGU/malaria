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
            cat motif_bed/${chrom}.motif.${peak_type}.bed |
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
            cmd="sort -k1,1 -k2,2n motif_bed/${chrom}.motif.${peak_type}.bed  | $bedtools merge -i stdin | awk '{print \$0\"\t\"NR}' | $bedtools intersect -a motif_bed/${chrom}.motif.${peak_type}.bed -b stdin -wb "
            $code_dir/s7_seperate_continuous_motif_nondegenerate.R "$(pwd)" "$cmd" "motif_bed/${chrom}.motif.${peak_type}.nonoverlap.bed"
        done
    done
}

find_specific_site() {
    for peak_type in $(echo "peak nopeak" | tr " " "\n"); do
        rm -rf $(find ../site* | grep ${peak_type}.bed)
        for chrom in $(echo $autosomes | tr " " "\n"); do
            mkdir -p ../site1 ../site2 ../site3 ../site4 ../site5 ../site6

            awk -v region=${peak_type} '$4=="+"{
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site1/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site2/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site3/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site4/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site5/"$1"."region".bed";
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site6/"$1"."region".bed";

                }$4=="-"{
                    print $1"\t"$3-1"\t"$3"\t"$4 >> "../site1/"$1"."region".bed";
                    print $1"\t"$3-2"\t"$3-1"\t"$4 >> "../site2/"$1"."region".bed";
                    print $1"\t"$3-3"\t"$3-2"\t"$4 >> "../site3/"$1"."region".bed";
                    print $1"\t"$3-4"\t"$3-3"\t"$4 >> "../site4/"$1"."region".bed";
                    print $1"\t"$3-5"\t"$3-4"\t"$4 >> "../site5/"$1"."region".bed";
                    print $1"\t"$3-6"\t"$3-5"\t"$4 >> "../site6/"$1"."region".bed";
                }' ${chrom}.motif.${peak_type}.nonoverlap.bed
        done
    done
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
flank_sites=$(seq 1 10)
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -e

# 1. 计算简并式motif在基因组的位置 # TODO: CHECK DONE
sed "s/motif_seq_template/GAWGAW/g" $code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh >$code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh_GAWGAW
sbatch $code_dir/s7_each_site_genome_index_for_motifs_asWhole.sh_GAWGAW

# 2. 合并所有motif，寻找flank区间 # TODO: CHECK DONE
cd $single_motif_global_output_dir/GAWGAW
find_flank_region

# 3.在所有motif水平上寻找nonoverlap的motif，先定两头 # TODO: CHECK DONE
cd $single_motif_global_output_dir/GAWGAW
find_nonoverlap_motif

# 4. 为每种motif定义1，2，3，4，5，6号位 # TODO: CHECK DONE
cd $single_motif_global_output_dir/GAWGAW/motif_bed
find_specific_site

# 5. 创建和之前相同的文件夹模式
cd $single_motif_global_output_dir/GAWGAW
ln -s $single_motif_global_output_dir/GAWGAW/flank_bed flank_nomatter_strand.remove_close_overlap

# 6. 计算motif位点的突变速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    for motif_seq in $(echo "GAWGAW" | tr " " "\n"); do
        sed "s/popu_symbol_template/${popu_symbol}/g" $code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh | sed "s/motif_seq_template/${motif_seq}/g" >$code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh_${popu_symbol}_${motif_seq}
        sbatch $code_dir/s7_compute_mutation_load_for_each_site_of_motif_asWhole.sh_${popu_symbol}_${motif_seq}
    done
done

# 7. 计算flank位点的突变速率
for popu_symbol in $(echo "ESEA_WSEA_OCE_SAM_SAS" | tr " " "\n"); do
    for motif_seq in $(echo "GAWGAW" | tr " " "\n"); do
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
for motif_seq in $(echo "GAWGAW" | tr " " "\n"); do
    sed "s/motif_template/$motif_seq/g" $code_dir/s7_compute_each_site_base_propor_asWhole.sh >$code_dir/s7_compute_each_site_base_propor_asWhole.sh_$motif_seq
    sbatch $code_dir/s7_compute_each_site_base_propor_asWhole.sh_$motif_seq
done
