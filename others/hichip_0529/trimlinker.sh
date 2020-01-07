#!/bin/bash
#cd /DATA2/pipeline/hichip/sample0529

Output=/DATA/work/lbyybl/ypjiang/hichip/sample0529_2018

trimLinker -t 12 -m 1 -k 1 -l 16 -o $Output/pol1_hichip_1 -n pol1_hichip_1 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT s-4-28-WT-Pol1-HiChIP-1_HL72KCCXY_L4_1.fq.gz s-4-28-WT-Pol1-HiChIP-1_HL72KCCXY_L4_2.fq.gz &
trimLinker -t 12 -m 1 -k 1 -l 16 -o $Output/pol1_hichip_2 -n pol1_hichip_2 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT s-4-28-WT-Pol1-HiChIP-2_HL72KCCXY_L4_1.fq.gz s-4-28-WT-Pol1-HiChIP-2_HL72KCCXY_L4_2.fq.gz &
trimLinker -t 12 -m 1 -k 1 -l 16 -o $Output/pol2_hichip_1 -n pol2_hichip_1 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT s-4-28-WT-Pol2-HiChIP-1_HL5F2CCXY_L6_1.fq.gz s-4-28-WT-Pol2-HiChIP-1_HL5F2CCXY_L6_2.fq.gz &
trimLinker -t 12 -m 1 -k 1 -l 16 -o $Output/pol2_hichip_2 -n pol2_hichip_2 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT s-4-28-WT-Pol2-HiChIP-2_HL5F2CCXY_L6_1.fq.gz s-4-28-WT-Pol2-HiChIP-2_HL5F2CCXY_L6_2.fq.gz &
