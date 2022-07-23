#!/bin/bash
#SBATCH -n 28
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_motif_and_flank_each_site_genome_position_each_chrom() {
    peak_type=$1
    chrom=$2

    mkdir -p $working_dir/motif_and_flank_each_site_genome_position
    cd $working_dir/motif_and_flank_each_site_genome_position

    motif_and_flank_each_site=${chrom}.${peak_type}.bed
    rm -rf $motif_and_flank_each_site
    cat $motif_bed_dir/${chrom}.motif.${peak_type}.nonoverlap.bed | while read line; do
        # head test | while read line; do
        one_motif_and_flank_each_site_genome_position=""
        # 左侧flank
        for flank_site in $(echo "10 9 8 7 6 5 4 3 2 1" | tr " " "\n"); do

            if [ $flank_site -eq 10 ]; then
                value=$(echo $"$line" | awk -v flank_site=${flank_site} '$4=="+"{
                        print $1"\t"$2-flank_site"\t"$2-flank_site+1
                    }$4=="-"{
                        print $1"\t"$3+flank_site-1"\t"$3+flank_site
                    }' | $bedtools intersect -loj -a stdin -b $whole_motif_and_flank_dir/flank_nomatter_strand.remove_close_overlap/flank_${flank_site}/${chrom}.${peak_type}.bed | uniq | awk '{print $1"\t"$5"\t"$6}')
                one_motif_and_flank_each_site_genome_position+="$value"
            else
                value=$(echo $"$line" | awk -v flank_site=${flank_site} '$4=="+"{
                        print $1"\t"$2-flank_site"\t"$2-flank_site+1
                    }$4=="-"{
                        print $1"\t"$3+flank_site-1"\t"$3+flank_site
                    }' | $bedtools intersect -loj -a stdin -b $whole_motif_and_flank_dir/flank_nomatter_strand.remove_close_overlap/flank_${flank_site}/${chrom}.${peak_type}.bed | uniq | awk '{print $5"\t"$6}')
                one_motif_and_flank_each_site_genome_position+="\t$value"
            fi
        done
        # motif位点
        for motif_site in $(seq 1 6); do
            value=$(echo $"$line" | awk -v motif_site=${motif_site} '$4=="+"{
                        print $1"\t"$2+motif_site-1"\t"$2+motif_site
                    }$4=="-"{
                        print $1"\t"$3-motif_site"\t"$3-motif_site+1
                    }' | $bedtools intersect -loj -a stdin -b $whole_motif_and_flank_dir/site${motif_site}/${chrom}.${peak_type}.bed | uniq | awk '{print $5"\t"$6}')
            one_motif_and_flank_each_site_genome_position+="\t$value"
        done
        # 右侧flank
        for flank_site in $(seq 1 10); do
            value=$(echo $"$line" | awk -v flank_site=${flank_site} '$4=="-"{
                        print $1"\t"$2-flank_site"\t"$2-flank_site+1
                    }$4=="+"{
                        print $1"\t"$3+flank_site-1"\t"$3+flank_site
                    }' | $bedtools intersect -loj -a stdin -b $whole_motif_and_flank_dir/flank_nomatter_strand.remove_close_overlap/flank_${flank_site}/${chrom}.${peak_type}.bed | uniq | awk '{print $5"\t"$6}')
            one_motif_and_flank_each_site_genome_position+="\t$value"
        done
        echo -e $"$one_motif_and_flank_each_site_genome_position" >>$motif_and_flank_each_site
    done
}

get_6mA_status_for_motif_and_flank_site() {
    peak_type=$1
    chrom=$2

    mkdir -p $working_dir/m6mA_status_motif_and_flank_site
    cd $working_dir/m6mA_status_motif_and_flank_site

    pos_dir=$working_dir/motif_and_flank_each_site_genome_position
    m6mA_status=$chrom.$peak_type.txt
    rm -rf $m6mA_status
    # cat test.bed | while read line; do
    cat $pos_dir/${chrom}.${peak_type}.bed | while read line; do
        m6mA_status_for_one_motif=""

        for i in $(seq 1 26); do
            let column_1=2*$i
            let column_2=2*$i+1
            suitable_pos_status=$(echo $"$line" | awk -v c1=$column_1 '{print $c1}')

            if [ $i == 1 ]; then
                if [ $suitable_pos_status -ne -1 ]; then
                    value=$(echo $"$line" | awk -v c1=$column_1 -v c2=$column_2 '{print $1"\t"$c1"\t"$c2}' | $bedtools intersect -a stdin -b $m6mA_3th_bed | wc -l)
                    if [ $value -eq 0 ] 2>/dev/null; then
                        m6mA_status_for_one_motif+="0"
                    else
                        m6mA_status_for_one_motif+="1"
                    fi
                else
                    m6mA_status_for_one_motif+="-1"
                fi
            else
                if [ $suitable_pos_status -ne -1 ]; then
                    value=$(echo $"$line" | awk -v c1=$column_1 -v c2=$column_2 '{print $1"\t"$c1"\t"$c2}' | $bedtools intersect -a stdin -b $m6mA_3th_bed | wc -l)
                    if [ $value -eq 0 ] 2>/dev/null; then
                        m6mA_status_for_one_motif+=" 0"
                    else
                        m6mA_status_for_one_motif+=" 1"
                    fi
                else
                    m6mA_status_for_one_motif+=" -1"
                fi
            fi
        done
        echo $"$m6mA_status_for_one_motif" >>$m6mA_status
    done
}
export -f get_motif_and_flank_each_site_genome_position_each_chrom get_6mA_status_for_motif_and_flank_site
# VARIABLE NAMING TODO:
export bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW"
export motif_bed_dir=$bed_dir/motif_bed
export whole_motif_and_flank_dir=$bed_dir
export global_working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/3th/public"
export working_dir=$global_working_dir/motif_each_site
export m6mA_3th_bed=$global_working_dir/3D7-6mA-9Cell_FRACgt20_COVgt25.bed
export parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $working_dir

time $parallel -j 28 get_motif_and_flank_each_site_genome_position_each_chrom ::: peak ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3

time $parallel -j 28 get_6mA_status_for_motif_and_flank_site ::: peak ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3
