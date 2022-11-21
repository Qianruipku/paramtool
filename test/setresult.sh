#!/bin/bash
allcases=`ls | grep ^[0-9]`
for case in $allcases
do
    cd $case
    mv result.txt result.ref
    cd ../
done