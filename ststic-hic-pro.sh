#!/bin/bash

# this script is used to extract the stastic info from the hic-pro result
## can't parallel you only to input the dir of hic-pro result and this script
## will produce a dir _stastic_ and in the sumsta, one file named final2_formal is the 
## stastic file; you can trans it into excel in the excel;

# Boyuan-Li
# 2018-5-27

function helps
{
        echo ""
        echo -e "Usage: $0 hic-pro-result"
   

}

while getopts "h" optionName
do
        case "$optionName" in
                
                h)
                        helps
                        exit 0
                        ;;
        esac
done

##########################################################################################
# $1 refore the hic-pro product dir;
##########################################################################################

mkdir _stastic_
ls -d $1/hic_results/data/sample* | parallel cp {}/{/}_allValidPairs.mergestat ./_stastic_/
ls -d $1/bowtie_results/bwt2/* | parallel -v cp {}/{/}*.bwt2pairs.pairstat ./_stastic_/
ls -d $1/hic_results/data/sample* | parallel -v cp {}/{/}*.mRSstat ./_stastic_/

cd ./_stastic_/
mkdir sumsta

ls sample*.mRSstat | parallel -v -j1 cat {} '|' egrep -w \"Valid_interaction_pairs\|Dangling_end_pairs\|Religation_pairs\|Self_Cycle_pairs\|Single-end_pairs\|Dumped_pairs\" >> ./sumsta/filter.stat
ls sample*_allValidPairs.mergestat | parallel -v -j1 cat {}  >> ./sumsta/interaction.stat
ls sample*.bwt2pairs.pairstat | parallel -v -j1 cat {} >> ./sumsta/pair.stat

cd sumsta

sed -i -e 's/^cat //g' -e 's/ |.*$//g'  filter.stat
sed -ie 's/.mRSstat//g' filter.stat
cat filter.stat | xargs  -l7 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n","experiment","Dangling_end_pairs", \
"Self_Cycle_pairs","Dumped_pairs","Valid_interaction_pairs","Religation_pairs","Single-end_pairs"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$5,$9,$13,$3,$7,$11}'>> format_filter

sed -i -e 's/^cat//g' -e 's/_allValidPairs.mergestat$//g' interaction.stat
cat interaction.stat | xargs -l7 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n", \
"experiment","valid_interaction_rmdup","cis_shortRange","cis_longRange","trans_interaction","valid_interaction","cis_interaction"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$5,$11,$13,$7,$3,$9}' >> format_interaction

# sed -i -e 's/^cat //g' -e 's/.mpairstat.*$//g' map.stat
# cat map.stat | xargs -l5 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s\n", \
# "experiment","Total_pairs_processed","Unmapped_pairs","Unique_singleton_alignments","Multiple_singleton_alignments"}; \
# {printf "%15s \t %15i \t %15i \t %15i \t %15i \n",$1,$3,$6,$9,$12; }' >> format_map2

sed -i -e 's/^cat //g' -e 's/_mm10.bwt2pairs.*$//g' pair.stat
cat pair.stat | xargs -l11 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \n","experiment","Total_pairs_processed","Unmapped_pairs","Pairs_with_singleton","Multiple_pairs_alignments","Unique_paired_alignments"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$3,$6,$18,$15,$12}' >> format_pair

join <(sort format_pair) <(sort format_filter )>> pair_filter
join <(sort pair_filter ) <(sort format_interaction ) >> pair_filter_interaction
cat pair_filter_interaction | awk 'NR==1{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s  \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >> final2_formal

cat pair_filter_interaction | awk 'NR>=2{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i\n \
%15s \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15f \t %15.5f \t %15.5f\t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18, \
$1,$2/$2,$3/$2,$4/$2,$5/$2,$6/$2,$7/$2,$8/$2,$9/$2,$10/$2,$11/$2,$12/$2,$13/$2,$14/$2,$15/$2,$16/$2,$17/$2,$18/$2}' >> final2_formal


