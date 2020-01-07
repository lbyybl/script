#!/bin/bash 
for j in $(seq 1 19) X Y;
do
    perl ../scripts/perl/matrix2insulation.pl -i chr${j}_0h_25kb_matrix.addedHeaders.matrix.gz --is 500000 --ids 250000 --im iqrMean --ss 125000 --yb 1.5 --bmoe 1 --nt 0.01 --ez --bg
    perl ../scripts/perl/insulation2tads.pl -i chr${j}_0h_25kb_matrix.addedHeaders--is500001--nt0.01--ids250001--ss125001--imiqrMean.insulation -b chr${j}_0h_25kb_matrix.addedHeaders--is500001--nt0.01--ids250001--ss125001--imiqrMean.insulation.boundaries
done
for k in $(seq 1 19) X Y;
do
    perl ../scripts/perl/matrix2compartment.pl -i chr${k}_0h_25kb_matrix.addedHeaders.matrix.gz --ca 0.005
done

