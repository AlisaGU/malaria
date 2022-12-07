#!/bin/bash
#SBATCH -n 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
summary() {
    file=$(ls | grep ${chrom}.window.${prefix}.bed-[0]*$order$)
    rm -rf ${chrom}.window.${prefix}.summary.${order}
    cat $file | while read line; do
        ## 计算窗口里有多少碱基位于motif里
        base_count_in_motif=$(echo "$line" | $bedtools intersect -a stdin -b $motif_genome_loc/${chrom}.genome.motif.bed | sort -k1,1 -k2,2n | $bedtools merge -i stdin | awk '{s+=$3-$2}END{print s}')
        if [[ $base_count_in_motif -gt 0 ]]; then
            :
        else
            base_count_in_motif=0
        fi

        ## 计算窗口里有几个motif
        motif_count_in_window=$(echo "$line" | $bedtools intersect -a stdin -b $motif_genome_loc/${chrom}.genome.motif.bed | wc -l)
        if [[ $motif_count_in_window -gt 0 ]]; then
            :
        else
            motif_count_in_window=0
        fi

        ## 计算窗口里的max FE、mean FE、mean FE
        let window_start=$(echo $line | awk '{print $2}')+1
        window_end=$(echo $line | awk '{print $3}')
        let window_start_minus_1=window_start-1
        let window_end_minus_1=window_end-1
        seq $window_start_minus_1 $window_end_minus_1 >${chrom}.${prefix}.${order}.a
        seq $window_start $window_end >${chrom}.${prefix}.${order}.b

        let specific_window_size=window_end-window_start_minus_1
        seq $specific_window_size | sed "c ${chrom}" | paste - ${chrom}.${prefix}.${order}.a ${chrom}.${prefix}.${order}.b >${chrom}.${prefix}.${order}.test
        $bedtools intersect -a $bdg_dir/${chrom}.treat.bdg -b ${chrom}.${prefix}.${order}.test >${chrom}.${prefix}.${order}.treat
        $bedtools intersect -a $bdg_dir/${chrom}.cont.bdg -b ${chrom}.${prefix}.${order}.test >${chrom}.${prefix}.${order}.cont
        FE_statistics=$(paste ${chrom}.${prefix}.${order}.treat ${chrom}.${prefix}.${order}.cont | awk '{print $4/$8}' | $datamash max 1 mean 1 median 1)

        ##总结
        echo -e $chrom"\t"$window_start_minus_1"\t"$window_end"\t"$base_count_in_motif"\t"$motif_count_in_window"\t"$FE_statistics >>${chrom}.window.${prefix}.summary.${order}
    done

    rm -rf ${chrom}.${prefix}.${order}.*
}
# VARIABLE NAMING TODO:
order=order_moban
chrom=chrom_moban
window_size=window_size_moban ## 单位：bp
step_size=step_size_moban     ## 单位：bp
output_dir=output_dir_moban
prefix=prefix_moban
motif_genome_loc=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/motif_GAWGAW_genome_loc
treat_bdg=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7-T3_treat_pileup.bdg
cont_bdg=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7-T3_control_lambda.bdg
bdg_dir=$output_dir/../../whole_genome_bdg

seqkit=/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
datamash=/picb/evolgen/users/gushanshan/software/datamash/bin/datamash
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
cd $output_dir
summary
