#!/usr/bin/bash
flist="${1}"
echo {} > result.json
cat result.json
for file in  $flist ; do 
    echo Doing $file
    compareJSON.py --or $file result.json > tmp
    mv tmp result.json
done
mv result.json $2
