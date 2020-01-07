#!/bin/bash
perl ../scripts/perl/generateBins.pl --a mm10 --r chr1:1-195471971 --bsize 25000 --bstep 1 chr1-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr1.headers.bed-x
mv myHeaders.headers.bed chr1.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr2:1-182113224 --bsize 25000 --bstep 1 chr2-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr2.headers.bed-x
mv myHeaders.headers.bed chr2.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr3:1-160039680 --bsize 25000 --bstep 1 chr3-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr3.headers.bed-x
mv myHeaders.headers.bed chr3.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr4:1-156508116 --bsize 25000 --bstep 1 chr4-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr4.headers.bed-x
mv myHeaders.headers.bed chr4.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr5:1-151834684 --bsize 25000 --bstep 1 chr5-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr5.headers.bed-x
mv myHeaders.headers.bed chr5.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr6:1-149736546 --bsize 25000 --bstep 1 chr6-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr6.headers.bed-x
mv myHeaders.headers.bed chr6.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr7:1-145441459 --bsize 25000 --bstep 1 chr7-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr7.headers.bed-x
mv myHeaders.headers.bed chr7.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr8:1-129401213 --bsize 25000 --bstep 1 chr8-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr8.headers.bed-x
mv myHeaders.headers.bed chr8.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr9:1-124595110 --bsize 25000 --bstep 1 chr9-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr9.headers.bed-x
mv myHeaders.headers.bed chr9.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr10:1-130694993 --bsize 25000 --bstep 1 chr10-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr10.headers.bed-x
mv myHeaders.headers.bed chr10.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr11:1-122082543 --bsize 25000 --bstep 1 chr11-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr11.headers.bed-x
mv myHeaders.headers.bed chr11.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr12:1-120129022 --bsize 25000 --bstep 1 chr12-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr12.headers.bed-x
mv myHeaders.headers.bed chr12.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr13:1-120421639 --bsize 25000 --bstep 1 chr13-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr13.headers.bed-x
mv myHeaders.headers.bed chr13.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr14:1-124902244 --bsize 25000 --bstep 1 chr14-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr14.headers.bed-x
mv myHeaders.headers.bed chr14.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr15:1-104043685 --bsize 25000 --bstep 1 chr15-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr15.headers.bed-x
mv myHeaders.headers.bed chr15.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr16:1-98207768 --bsize 25000 --bstep 1 chr16-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr16.headers.bed-x
mv myHeaders.headers.bed chr16.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr17:1-94987271 --bsize 25000 --bstep 1 chr17-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr17.headers.bed-x
mv myHeaders.headers.bed chr17.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr18:1-90702639 --bsize 25000 --bstep 1 chr18-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr18.headers.bed-x
mv myHeaders.headers.bed chr18.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chr19:1-61431566 --bsize 25000 --bstep 1 chr19-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chr19.headers.bed-x
mv myHeaders.headers.bed chr19.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chrX:1-171031299 --bsize 25000 --bstep 1 chrX-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chrX.headers.bed-x
mv myHeaders.headers.bed chrX.headers.bed-y
perl ../scripts/perl/generateBins.pl --a mm10 --r chrY:1-91744698 --bsize 25000 --bstep 1 chrY-0h
sed -i '1d' myHeaders.headers.bed
cut -f4 myHeaders.headers.bed > chrY.headers.bed-x
mv myHeaders.headers.bed chrY.headers.bed-y
perl ../scripts/perl/addMatrixHeaders.pl -i chr1_0h_25kb_matrix --xhf chr1.headers.bed-x --yhf chr1.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr2_0h_25kb_matrix --xhf chr2.headers.bed-x --yhf chr2.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr3_0h_25kb_matrix --xhf chr3.headers.bed-x --yhf chr3.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr4_0h_25kb_matrix --xhf chr4.headers.bed-x --yhf chr4.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr5_0h_25kb_matrix --xhf chr5.headers.bed-x --yhf chr5.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr6_0h_25kb_matrix --xhf chr6.headers.bed-x --yhf chr6.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr7_0h_25kb_matrix --xhf chr7.headers.bed-x --yhf chr7.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr8_0h_25kb_matrix --xhf chr8.headers.bed-x --yhf chr8.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr9_0h_25kb_matrix --xhf chr9.headers.bed-x --yhf chr9.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr10_0h_25kb_matrix --xhf chr10.headers.bed-x --yhf chr10.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr11_0h_25kb_matrix --xhf chr11.headers.bed-x --yhf chr11.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr12_0h_25kb_matrix --xhf chr12.headers.bed-x --yhf chr12.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr13_0h_25kb_matrix --xhf chr13.headers.bed-x --yhf chr13.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr14_0h_25kb_matrix --xhf chr14.headers.bed-x --yhf chr14.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr15_0h_25kb_matrix --xhf chr15.headers.bed-x --yhf chr15.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr16_0h_25kb_matrix --xhf chr16.headers.bed-x --yhf chr16.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr17_0h_25kb_matrix --xhf chr17.headers.bed-x --yhf chr17.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr18_0h_25kb_matrix --xhf chr18.headers.bed-x --yhf chr18.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chr19_0h_25kb_matrix --xhf chr19.headers.bed-x --yhf chr19.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chrX_0h_25kb_matrix --xhf chrX.headers.bed-x --yhf chrX.headers.bed-x
perl ../scripts/perl/addMatrixHeaders.pl -i chrY_0h_25kb_matrix --xhf chrY.headers.bed-x --yhf chrY.headers.bed-x
