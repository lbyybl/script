#!/usr/bin/bash

###########################################################################################################################
## to extract 1000000 reads from the .gz files
###########################################################################################################################
ls sample0419/* sample1018/* sample1026/* sample1124/* | parallel mkdir /DATA/work/lbyybl/wangcl/hic/KPNA2/5_21_test/{//}
 
ls sample0419/* sample1018/* sample1026/* sample1124/* | parallel --dryrun head -4000000 \<\(less {}\) '>' \
/DATA/work/lbyybl/wangcl/hic/KPNA2/5_21_test/{.}

###########################################################################################################################
## to trim the linker insert the reads #batch is except sample0409; means the sample0409 is treated only;
###########################################################################################################################
trimLinker -t 12 -m 1 -k 1 -l 16 -o ./6h_1 -n 6h_1 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT s1215Rpb1-6h-1_HCKYCCCXY_L2_1.fq.gz \
 s1215Rpb1-6h-1_HCKYCCCXY_L2_2.fq.gz

cat filetotrimlinker.txt | parallel -C "\s" trimlinker -t 2 -m 1 -k 1 -l 16  -o ./trimLinker/{1//} -n {1//} \
-A CGGTGGCT -B GCCACCGT {1} {2}

ls sample0419/* | xargs -l2 | parallel --dryrun -C "\s" trimlinker -t 2 -m 1 -k 1 -l 16  -o ./trimLinker/{1//} -n {1//} \
-A CGCGATATCTTATCTGACT -B GTCAGATAAGATATCGCGT {1} {2}





###########################################################################################################################
## the stastic of first step
###########################################################################################################################
ls -d ./5-21result/hic_results/data/sample* | parallel cp {}/{/}_allValidPairs.mergestat ./stastic/
ls sample*.trim.stat | parallel -v -j1 cat {} '|' egrep \"Total\|Valid\" >> ./sumsta/trim.stat
ls sample*.mpairstat | parallel -v -j1 cat {} '|' egrep \"Total_pairs_processed\|Unmapped_pairs\|Unique_singleton_alignments\|Multiple_singleton_alignments\" >> ./sumsta/map.stat
ls sample*.mRSstat | parallel -v -j1 cat {} '|' egrep -w \"Valid_interaction_pairs\|Dangling_end_pairs\|Religation_pairs\|Self_Cycle_pairs\|Single-end_pairs\|Dumped_pairs\" >> ./sumsta/filter.stat
ls sample*_allValidPairs.mergestat | parallel -v -j1 cat {}  >> ./sumsta/interaction.stat
sed -i -e 's/^cat //g' -e 's/ |.*$//g'  filter.stat
sed -ie 's/.mRSstat//g' filter.stat
cat filter.stat | xargs  -l7 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n","experiment","Dangling_end_pairs", \
"Self_Cycle_pairs","Dumped_pairs","Valid_interaction_pairs","Religation_pairs","Single-end_pairs"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$5,$9,$13,$3,$7,$11}'>>format_filter

sed -i -e 's/^cat//g' -e 's/_allValidPairs.mergestat$//g' interaction.stat
cat interaction.stat | xargs -l7 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n", \
"experiment","valid_interaction_rmdup","cis_shortRange","cis_longRange","trans_interaction","valid_interaction","cis_interaction"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$5,$11,$13,$7,$3,$9}' >> format_interaction

sed -i -e 's/^cat //g' -e 's/.mpairstat.*$//g' map.stat
cat map.stat | xargs -l5 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s\n", \
"experiment","Total_pairs_processed","Unmapped_pairs","Unique_singleton_alignments","Multiple_singleton_alignments"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \n",$1,$3,$6,$9,$12; }' >> format_map2
#printf "%15s \t %10.5f \t %10.5f \t %10.5f \t %10.5f \n",$1,$4,$7,$10,$13}' >> format_map

sed -i -e "s/^cat //g" -e "s/\.trim.stat.*$//g" trim.stat
cat trim.stat | xargs -l3 | awk 'BEGIN{printf "%15s \t %15s \t %15s \n","experiment","Total_PETs","Valid_PETs"}; \
{printf "%15s \t %15i \t %15i \n",$1,$4,$7}' >> format_trim

sed -i -e 's/^cat //g' -e 's/_mm10.bwt2pairs.*$//g' pair.stat
cat pair.stat | xargs -l11 | awk 'BEGIN{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \n","experiment","Total_pairs_processed","Unmapped_pairs","Pairs_with_singleton","Multiple_pairs_alignments","Unique_paired_alignments"}; \
{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i\n",$1,$3,$6,$18,$15,$12}' >> format_pair


join <(sort format_map2) <(sort format_filter )>> map2_filter
join <(sort map2_filter ) <(sort format_interaction ) >> map2_filter_interaction


#cat map2_filter_interaction | awk 'NR>=2{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t \
#%15i \t  %15i \t %15i \t %15i \t %15i \t %15i\n %15s \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\t \
#%15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17, \
#$1,$2/$2,$3/$2,$4/$2,$5/$2,$6/$2,$7/$2,$8/$2,$9/$2,$10/$2,$11/$2,$12/$2,$13/$2,$14/$2,$15/$2,$16/$2,$17/$2}' | less

cat map2_filter_interaction | awk 'NR==1{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s  \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' >> final_formal

cat map2_filter_interaction | awk 'NR>=2{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i\n \
%15s \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17, \
$1,$2/$2,$3/$2,$4/$2,$5/$2,$6/$2,$7/$2,$8/$2,$9/$2,$10/$2,$11/$2,$12/$2,$13/$2,$14/$2,$15/$2,$16/$2,$17/$2}' >> final_formal

ls -d ./5-21result/bowtie_results/bwt2/* | parallel -v --dryrun cp {}/{/}*.bwt2pairs.pairstat ./stastic/
join <(sort format_pair) <(sort format_filter )>> pair_filter
join <(sort pair_filter ) <(sort format_interaction ) >> pair_filter_interaction
cat pair_filter_interaction | awk 'NR==1{printf "%15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s  \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s \t %15s\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >> final2_formal

cat pair_filter_interaction | awk 'NR>=2{printf "%15s \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i \t  %15i \t %15i \t %15i \t %15i \t %15i\n \
%15s \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15f \t %15.5f \t %15.5f\t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f \t %15.5f\n", \
$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18, \
$1,$2/$2,$3/$2,$4/$2,$5/$2,$6/$2,$7/$2,$8/$2,$9/$2,$10/$2,$11/$2,$12/$2,$13/$2,$14/$2,$15/$2,$16/$2,$17/$2,$18/$2}' >> final2_formal


###########################################################################################################################
##
###########################################################################################################################

###########################################################################################################################
##
###########################################################################################################################