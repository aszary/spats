#! /bin/bash
# uses PSRSALSA
pmod -debase -onpulse '473 556' ../input/1.spCF
pspec -w -2dfs -lrfs -nfft 256 -onpulse '473 556' ../input/1.debase.gg # 2DFS and LRFS creates 
pspecDetect -v ../input/1.debase.gg # detects P_3 # P3[P0]  = 4.776026 +- 0.006533
pfold -p3fold "4.78 18" -p3fold_nritt 50 -p3fold_cpb 50 -w -onpulse '473 556' -oformat ascii ../input/1.debase.gg # P3-folded profile in file 1.debase.p3fold
