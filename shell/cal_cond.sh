#!/bin/bash
ele=al
density=0.215
allT=(1 2 3 4 6 8 12 16 25 32 50 64 100 128 200 256 384 512 750 1000)
for T in ${allT[@]}
do
cat>input<<EOF
3
$ele
$density
$T
EOF
tool.exe<input >_tmp
sigma=`grep "electrical conductivity" _tmp | awk '{print $3}' `
kappa=`grep "thermal conductivity" _tmp | awk '{print $3}' `
rm -f _tmp
echo $T $sigma $kappa
done

rm -f input