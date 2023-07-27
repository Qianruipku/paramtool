#!/bin/bash
ele=al
density=2.7
allT=(1 2 3 4 6 8 12 16 25 32 50 64 100 128 200 256 384 512 750 1000)
model="Spitzer"
# model="Lee-More"
for T in ${allT[@]}
do
cat>input<<EOF
3
$ele
$density
$T
EOF
../tool.exe<input >_tmp
sigma=`grep $model _tmp | awk '{print $2}' `
kappa=`grep $model _tmp | awk '{print $3}' `
Lorentz=`grep $model _tmp | awk '{print $4}' `
rm -f _tmp
echo $T $sigma $kappa $Lorentz
done

rm -f input