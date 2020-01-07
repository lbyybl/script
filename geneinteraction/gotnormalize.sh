#/bin/bash
Resolution=$1



join -1 1 -2 1 <(sort -k 1,1 -k 2,2 interaction_ko${Resolution}K.txt) <(sort -k1 region${Resolution}K.bed) | sort -k 2,2 -k 1,1 -o interactionnomalize_ko${Resolution}K
join -1 2 -2 1 interactionnomalize_ko${Resolution}K <(sort -k1 region${Resolution}K.bed) | sort -o interactionnomalize_ko${Resolution}K
join -1 1 -2 1 <(sort -k 1,1 -k 2,2 interaction_wt${Resolution}K.txt) <(sort -k1 region${Resolution}K.bed) | sort -k 2,2 -k 1,1 -o interactionnomalize_wt${Resolution}K
join -1 2 -2 1 interactionnomalize_wt${Resolution}K <(sort -k1 region${Resolution}K.bed) | sort -o interactionnomalize_wt${Resolution}K
