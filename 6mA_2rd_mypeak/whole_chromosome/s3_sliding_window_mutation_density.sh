#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
# FUNCTIONS TODO:
get_average_variant_count() {
    var_type=$1

    for por in $(echo $relative_por | tr " " "\n"); do
        mutation_density_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$por/mutation_density
        mkdir -p $mutation_density_dir
        mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
        rm -rf $mean_variant_count

        for chrom in $(echo $autosomes | tr " " "\n"); do
            chrom_var_count=""

            regions="$chrom_bed_dir/${chrom}.nopeak.bed"

            for w in $(seq 1 $window_number); do
                regions+=" $chrom_bed_dir/${sliding_window_type}/relative.${por}/${chrom}/w$w.bed"
            done
            for region_bed in $(echo $regions | tr " " "\n"); do
                if [[ $(basename $region_bed) =~ "nopeak" ]] || [ $sliding_window_type == "sliding_window" ]; then
                    chrom_var_count+=" $(get_average_variant_count_for_chrom $chrom $region_bed)"
                else
                    chrom_var_count+=" $(get_average_variant_count_for_chrom_mean_window $chrom $region_bed)"
                fi
            done
            echo $chrom_var_count >>$mean_variant_count
        done
    done

}

get_average_variant_count_3control_N2N_indel() {
    var_type=$1

    for por in $(echo $relative_por | tr " " "\n"); do
        mutation_density_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$por/mutation_density
        mkdir -p $mutation_density_dir
        mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_3control
        rm -rf $mean_variant_count

        for chrom in $(echo $autosomes | tr " " "\n"); do
            chrom_var_count=""
            if [ $var_type == "insertion" ] || [ $var_type == "deletion" ]; then
                chrom_var_count+=" $(get_average_variant_count_for_chrom $chrom $chrom_bed_dir/${chrom}.nopeak.bed)"
            else
                chrom_var_count="$(get_average_variant_count_for_chrom_base_count $chrom $chrom_bed_dir/${chrom}.nopeak.bed $var_type)"
            fi

            if [ $var_type == "insertion" ] || [ $var_type == "deletion" ]; then
                chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif $chrom $chrom_bed_dir/${chrom}.nopeak.bed)"
            else
                chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_base_count $chrom $chrom_bed_dir/${chrom}.nopeak.bed $var_type)"
            fi

            if [ $var_type == "insertion" ] || [ $var_type == "deletion" ]; then
                chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif $chrom $chrom_bed_dir/$chrom.remove_margin_500_nopeak.bed)"
            else
                chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif_base_count $chrom $chrom_bed_dir/$chrom.remove_margin_500_nopeak.bed $var_type)"
            fi

            regions=""
            for w in $(seq 1 $window_number); do
                regions+=" $chrom_bed_dir/${sliding_window_type}/relative.${por}/${chrom}/w$w.bed"
            done
            for region_bed in $(echo $regions | tr " " "\n"); do
                if [ $var_type == "insertion" ] || [ $var_type == "deletion" ]; then
                    chrom_var_count+=" $(get_average_variant_count_for_chrom $chrom $region_bed)"
                else
                    chrom_var_count+=" $(get_average_variant_count_for_chrom_base_count $chrom $region_bed $var_type)"
                fi
            done
            echo $chrom_var_count >>$mean_variant_count
        done
    done

}

get_average_variant_count_only_control() {
    var_type=$1

    for por in $(echo $relative_por | tr " " "\n"); do
        mutation_density_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$por/mutation_density
        mkdir -p $mutation_density_dir
        mean_variant_count=$mutation_density_dir/${var_type}_only_control_mean_variant_count
        rm -rf $mean_variant_count

        for chrom in $(echo $autosomes | tr " " "\n"); do
            chrom_var_count="$(get_average_variant_count_for_chrom $chrom $chrom_bed_dir/${chrom}.nopeak.bed)"
            chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif $chrom $chrom_bed_dir/${chrom}.nopeak.bed)"
            chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif $chrom $chrom_bed_dir/$chrom.remove_margin_500_nopeak.bed)"

            echo $chrom_var_count >>$mean_variant_count
        done
    done

}

get_average_variant_count_only_peak() {
    var_type=$1

    for por in $(echo $relative_por | tr " " "\n"); do
        mutation_density_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$por/mutation_density
        mkdir -p $mutation_density_dir
        mean_variant_count=$mutation_density_dir/${var_type}_only_peak_mean_variant_count
        rm -rf $mean_variant_count

        for chrom in $(echo $autosomes | tr " " "\n"); do
            chrom_var_count="$(get_average_variant_count_for_chrom $chrom $chrom_bed_dir/${chrom}.peak.bed)"
            chrom_var_count+=" $(get_average_variant_count_for_chrom_only_motif $chrom $chrom_bed_dir/${chrom}.peak.bed)"

            echo $chrom_var_count >>$mean_variant_count
        done
    done

}

get_average_variant_count_for_chrom_only_motif() {
    chrom=$1
    region_bed=$2
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    motif_dir=""
    if [ $(basename $region_bed | grep "nopeak" | wc -l) -eq 1 ]; then
        motif_dir=$(echo $nopeak_motif_bed_dir)
    else
        motif_dir=$(echo $peak_motif_dir)
    fi

    awk '{print $0"\t""motif_"NR}' ${motif_dir}/${chrom}.* | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' | $bedtools intersect -wa -a stdin -b $region_bed >$motif_dir/${chrom}.intermediate

    chrom_region_len=$(awk '{print $3-$2}' $motif_dir/${chrom}.intermediate | awk '{s+=$1} END {print s}')
    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1

    rm -rf $motif_dir/${chrom}.intermediate
}

get_average_variant_count_for_chrom() {
    chrom=$1
    region_bed=$2
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

get_average_variant_count_for_chrom_base_count() {
    chrom=$1
    region_bed=$2
    var_type=$3

    # chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    ori_base=$(echo $var_type | awk -F"2" '{print $1}')
    chrom_region_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${ori_base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

get_average_variant_count_for_chrom_only_motif_base_count() {
    chrom=$1
    region_bed=$2
    var_type=$3

    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    ori_base=$(echo $var_type | awk -F"2" '{print $1}')

    motif_dir=""
    if [ $(basename $region_bed | grep "nopeak" | wc -l) -eq 1 ]; then
        motif_dir=$(echo $nopeak_motif_bed_dir)
    else
        motif_dir=$(echo $peak_motif_dir)
    fi

    awk '{print $0"\t""motif_"NR}' ${motif_dir}/${chrom}.* | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' | $bedtools intersect -wa -a stdin -b $region_bed >$motif_dir/${chrom}.intermediate

    chrom_region_len=$($seqkit subseq --quiet --bed $motif_dir/${chrom}.intermediate $genome | $seqkit fx2tab -i -H -n -C "${ori_base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

    # chrom_region_len=$(awk '{print $3-$2}' $motif_dir/${chrom}.intermediate | awk '{s+=$1} END {print s}')
    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $motif_dir/${chrom}.intermediate | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1

    rm -rf $motif_dir/${chrom}.intermediate
}
get_average_variant_count_for_chrom_mean_window() {
    chrom=$1
    region_bed=$2
    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    all_mean=$(echo $region_bed | sed 's/bed/mean.txt/g')
    rm -rf $all_mean

    while read -r window_bed; do
        b=$(get_one_window_mutation_density_chrom "$window_bed")
        echo $b >>$all_mean
    done <"$region_bed"

    mean_window_mutation_density=$($code_dir/ss3_1_average.r $all_mean)
    echo $mean_window_mutation_density
}

get_one_window_mutation_density_chrom() {
    window_bed=$1
    chrom_region_len=$(echo $window_bed | awk '{print $3-$2}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $window_bed | sed 's/ /\t/g' | $bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b - | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

compute_base_proporation() {
    base=$1
    base_proporation=$base_proporation_dir/${base}_mean_proporation
    rm -rf $base_proporation

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_base_propor=""

        regions="$chrom_bed_dir/${chrom}.nopeak.bed"

        for w in $(seq 1 $window_number); do
            regions+=" $chrom_bed_dir/${sliding_window_type}/relative.${relative_por}/${chrom}/w$w.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            whole_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
            base_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

            base_proportion=$(echo "scale=10;${base_len}/${whole_region_len}" | bc)
            base_proportion1=$(printf "%0.6f" $base_proportion)
            chrom_base_propor+=" $(echo $base_proportion1)"
        done

        echo $chrom_base_propor >>$base_proporation
    done
}

compute_base_proporation_abso_10() {
    base=$1
    base_proporation=$base_proporation_dir/${base}_mean_proporation
    rm -rf $base_proporation

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_base_propor=""

        regions="$chrom_bed_dir/${chrom}.nopeak.bed"
        regions+=" $chrom_bed_dir/absolute.10/${chrom}.peak.absolute.10.bed"

        # for w in $(seq 1 $window_number); do
        #     regions+=" $chrom_bed_dir/${sliding_window_type}/relative.${relative_por}/${chrom}/w$w.bed"
        # done

        for region_bed in $(echo $regions | tr " " "\n"); do
            whole_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
            base_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

            base_proportion=$(echo "scale=10;${base_len}/${whole_region_len}" | bc)
            base_proportion1=$(printf "%0.6f" $base_proportion)
            chrom_base_propor+=" $(echo $base_proportion1)"
        done

        echo $chrom_base_propor >>$base_proporation
    done
}
# VARIABLE NAMING TODO:
popu_symbol=$1
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd/whole_chromosome"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
variant_dir=$global_output_dir/$popu_symbol/variant
chrom_bed_dir=$macs2_out_dir/chrom

relative_por="0.05"
window_number="10"
motif_id="2"
# sliding_window_type="sliding_window_noCollapseWindow"
sliding_window_type="sliding_window"

# base_proporation_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/relative$relative_por/base_proporation
base_proporation_dir=$global_output_dir/$popu_symbol/${sliding_window_type}/absolute10/base_proporation

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
nopeak_motif_bed_dir=$macs2_two_outgroup_consistent_dir/motif_pattern_${motif_id}_nopeak/motif_loc
peak_motif_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_${motif_id}/motif_loc"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
mkdir -p $base_proporation_dir
# mkdir -p $mutation_density_dir
# get_average_variant_count "allvar"
# get_average_variant_count "snps"
# get_average_variant_count "indels"
# get_average_variant_count "A2T"
# get_average_variant_count "A2G"
# get_average_variant_count "A2C"
# get_average_variant_count "T2A"
# get_average_variant_count "T2G"
# get_average_variant_count "T2C"
# get_average_variant_count "C2T"
# get_average_variant_count "C2G"
# get_average_variant_count "C2A"
# get_average_variant_count "G2A"
# get_average_variant_count "G2C"
# get_average_variant_count "G2T"

# compute_base_proporation "A"
# compute_base_proporation "T"
# compute_base_proporation "C"
# compute_base_proporation "G"

# compute_base_proporation_abso_10 "A"
# compute_base_proporation_abso_10 "T"
# compute_base_proporation_abso_10 "C"
# compute_base_proporation_abso_10 "G"

# get_average_variant_count_only_control "snps"
# get_average_variant_count_only_peak "snps"

get_average_variant_count_3control_N2N_indel "A2T"
get_average_variant_count_3control_N2N_indel "A2G"
get_average_variant_count_3control_N2N_indel "A2C"
get_average_variant_count_3control_N2N_indel "T2A"
get_average_variant_count_3control_N2N_indel "T2G"
get_average_variant_count_3control_N2N_indel "T2C"
get_average_variant_count_3control_N2N_indel "C2T"
get_average_variant_count_3control_N2N_indel "C2G"
get_average_variant_count_3control_N2N_indel "C2A"
get_average_variant_count_3control_N2N_indel "G2A"
get_average_variant_count_3control_N2N_indel "G2C"
get_average_variant_count_3control_N2N_indel "G2T"
get_average_variant_count_3control_N2N_indel "insertion"
get_average_variant_count_3control_N2N_indel "deletion"
