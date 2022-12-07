#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

get_flank_site_by_strand_closest() {
    for peak_type in $(echo "nopeak peak" | tr " " "\n"); do
        for flank_type in $(echo $flank_types | tr " " "\n"); do
            rm -rf $(find $single_motif_output_dir/flank_${flank_type} -name ${peak_type}.allChroms.bed 2>/dev/null)
            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                mkdir -p $flank_dir
            done

            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                intermediate=$flank_dir/${peak_type}.allChroms.bed.intermediate
                if [ $flank_type == "upstream" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate

                    sort -k1,1 -k2,2n $intermediate | uniq | $bedtools intersect -a stdin -b $macs2_two_outgroup_consistent_dir/core_consistent_${peak_type}.bed | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | $bedtools closest -b $single_motif_global_output_dir/${motif_seq}.${peak_type}.bed -a stdin -d -D "b" -id | awk -v out=$flank_dir -v type=${peak_type} '{a=sqrt($NF^2);print $1"\t"$2"\t"$3>> out"/../flank_"a"/"type".allChroms.bed"}'

                elif [ $flank_type == "downstream" ]; then
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate

                    sort -k1,1 -k2,2n $intermediate | uniq | $bedtools intersect -a stdin -b $macs2_two_outgroup_consistent_dir/core_consistent_${peak_type}.bed | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | $bedtools closest -b $single_motif_global_output_dir/${motif_seq}.${peak_type}.bed -a stdin -d -D "b" -iu | awk -v out=$flank_dir -v type=${peak_type} '{a=sqrt($NF^2);print $1"\t"$2"\t"$3>> out"/../flank_"a"/"type".allChroms.bed"}'

                elif [ $flank_type == "nomatter_strand" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate

                    sort -k1,1 -k2,2n $intermediate | uniq | $bedtools intersect -a stdin -b $macs2_two_outgroup_consistent_dir/core_consistent_${peak_type}.bed | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | $bedtools closest -b $single_motif_global_output_dir/${motif_seq}.${peak_type}.bed -a stdin -d -D "b" | awk -v out=$flank_dir -v type=${peak_type} '{a=sqrt($NF^2);print $1"\t"$2"\t"$3>> out"/../flank_"a"/"type".allChroms.bed"}'
                fi
                rm -rf $intermediate
            done

            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                cd $flank_dir
                sort -k1,1 -k2,2n ${peak_type}.allChroms.bed | uniq | awk -v type=${peak_type} '{print >$1"."type".bed"}'
            done
        done
    done
}

get_flank_site_by_strand_remove_overlap() {
    for peak_type in $(echo "nopeak peak" | tr " " "\n"); do
        for flank_type in $(echo "nomatter_strand" | tr " " "\n"); do
            dupli_all_related_site=$single_motif_output_dir/flank_${flank_type}/dupli_among_any_sites_${peak_type}
            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                mkdir -p $flank_dir
                intermediate=$flank_dir/${peak_type}.allChroms.bed.intermediate
                if [ $flank_type == "upstream" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate

                elif [ $flank_type == "downstream" ]; then
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate
                elif [ $flank_type == "nomatter_strand" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate
                fi
            done
            cat $(find $single_motif_output_dir/flank_${flank_type}/ -name *.intermediate) | sort -k1,1 -k2,2n | uniq -d >$dupli_all_related_site
            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                cd $flank_dir
                sort -k1,1 -k2,2n ${peak_type}.allChroms.bed.intermediate | $bedtools subtract -a stdin -b $dupli_all_related_site | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | uniq | awk -v type=${peak_type} '{print >$1"."type".bed"}'
                rm -rf *.intermediate
            done
            rm -rf $dupli_all_related_site
        done
    done
}

get_flank_site_by_strand_remove_close_overlap() {
    for peak_type in $(echo "nopeak peak" | tr " " "\n"); do
        for flank_type in $(echo "nomatter_strand.remove_close_overlap" | tr " " "\n"); do
            dupli_all_related_site=$single_motif_output_dir/flank_${flank_type}/dupli_among_any_sites_${peak_type}
            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                mkdir -p $flank_dir
                intermediate=$flank_dir/${peak_type}.allChroms.bed.intermediate
                if [ $flank_type == "upstream" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate

                elif [ $flank_type == "downstream" ]; then
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate
                elif [ $flank_type == "nomatter_strand.remove_close_overlap" ]; then
                    cat $site1_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >$intermediate
                    cat $site1_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.plus.* | awk -v flank_site=$flank_site '{print $1"\t"$2+flank_site"\t"$3+flank_site}' >>$intermediate
                    cat $site6_dir/*.${peak_type}.minus.* | awk -v flank_site=$flank_site '{print $1"\t"$2-flank_site"\t"$3-flank_site}' >>$intermediate
                fi
            done
            cat $(find $single_motif_output_dir/flank_${flank_type}/flank_[1-6] -name *.intermediate) | sort -k1,1 -k2,2n | uniq -d >$dupli_all_related_site
            for flank_site in $(echo $flank_sites | tr " " "\n"); do
                flank_dir=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}
                cd $flank_dir
                # sort -k1,1 -k2,2n ${peak_type}.allChroms.bed.intermediate | $bedtools subtract -a stdin -b $dupli_all_related_site | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | uniq | awk -v type=${peak_type} '{print >$1"."type".bed"}'
                sort -k1,1 -k2,2n ${peak_type}.allChroms.bed.intermediate | $bedtools subtract -a stdin -b $dupli_all_related_site | $bedtools subtract -a stdin -b $single_motif_global_output_dir/all_motifs_position.${peak_type}.bed | uniq | $bedtools intersect -a stdin -b $core_consistent_dir/core_consistent_${peak_type}.bed | awk -v type=${peak_type} '{print >$1"."type".bed"}'

                rm -rf *.intermediate
            done
            rm -rf $dupli_all_related_site
        done
    done
}

get_average_variant_count_for_each_motif_global_variants() {
    flank_type=$1
    flank_site=$2
    get_average_variant_count_for_global_variants "allvar"
    get_average_variant_count_for_global_variants "snps"
    get_average_variant_count_for_global_variants "indels"
    get_average_variant_count_for_global_variants "single_del"

}

get_average_variant_count_for_global_variants() {
    var_type=$1
    specific_mutation_density_dir=$mutation_density_dir/flank_${flank_type}/flank_${flank_site}
    mean_variant_count=$specific_mutation_density_dir/flank_${flank_site}.${var_type}_mean_variant_count
    mkdir -p $specific_mutation_density_dir
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        peak_file=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}/${chrom}.peak.*
        base_count_on_this_site_peak=$(cat $peak_file | wc -l)
        chrom_var_count_peak="$(get_average_variant_count_for_chrom_base_count $chrom $peak_file $var_type)"

        nopeak_file=$single_motif_output_dir/flank_${flank_type}/flank_${flank_site}/${chrom}.nopeak.*
        base_count_on_this_site_nopeak=$(cat $nopeak_file | wc -l)
        chrom_var_count_nopeak="$(get_average_variant_count_for_chrom_base_count $chrom $nopeak_file $var_type)"

        echo $chrom_var_count_peak $base_count_on_this_site_peak $chrom_var_count_nopeak $base_count_on_this_site_nopeak >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    varType=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_mutation_count=""
    if [ $varType == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_insert" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_multi_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "A_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "T_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "C_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "G_indel" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    elif [ $varType == "single_del" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $26+$27+$28+$29}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        $chrom_mutation_count=0
    fi

    echo $chrom_mutation_count
}

# VARIABLE NAMING TODO:
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
motif_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs"
motif_seq=motif_seq_template
popu_symbol=popu_symbol_template

# motif_seq="GAHGAW"
motif=$motif_dir/$motif_seq

global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
variant_dir=$global_output_dir/$popu_symbol/variant_two_group
mutation_density_dir=$global_output_dir"/${popu_symbol}/motif/mutation_dentisy_two_outgroup_consistent/each_site/$motif_seq"

macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
single_motif_output_dir=$single_motif_global_output_dir/$motif_seq
motif_bed_dir=$single_motif_output_dir/motif_bed
flank_sites=$(seq 1 12)
core_consistent_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent"

flank_types="upstream downstream nomatter_strand"
site1_dir=$single_motif_output_dir/site1
site6_dir=$single_motif_output_dir/site6

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
# get_flank_site_by_strand_closest
# for flank_type in $(echo $flank_types | tr " " "\n"); do
#     for flank_site in $(echo $flank_sites | tr " " "\n"); do
#         get_average_variant_count_for_each_motif_global_variants ${flank_type} ${flank_site}
#     done
# done

# get_flank_site_by_strand_remove_overlap
# for flank_type in $(echo "nomatter_strand" | tr " " "\n"); do
#     for flank_site in $(echo $flank_sites | tr " " "\n"); do
#         get_average_variant_count_for_each_motif_global_variants ${flank_type} ${flank_site}
#     done
# done

get_flank_site_by_strand_remove_close_overlap
for flank_type in $(echo "nomatter_strand.remove_close_overlap" | tr " " "\n"); do
    for flank_site in $(echo $flank_sites | tr " " "\n"); do
        get_average_variant_count_for_each_motif_global_variants ${flank_type} ${flank_site}
    done
done
