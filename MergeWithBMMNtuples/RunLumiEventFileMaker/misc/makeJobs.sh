files=`ls srcFiles/bmmg*.fls `
for f in $files: 
    do
    f=`echo $f | sed 's/://' `
    echo $f
    fname=`echo $f | awk -F'/' '{print $NF}' `
    tag=`echo $fname |sed 's/.fls//' | awk -F'_' '{print $NF}'`
    ftype=`echo $fname |sed 's/.fls//' | awk -F'_' '{print $1}'`
    echo $tag
    set -x
    ./jobFileMakerBMMG.py $f RunLumiFiles/${ftype}_mc_bkg_${tag}/ 1 2000 $tag
    set +x
    done

files=`ls srcFiles/bmm5*.fls `
for f in $files: 
    do
    f=`echo $f | sed 's/://' `
    echo $f
    fname=`echo $f | awk -F'/' '{print $NF}' `
    tag=`echo $fname |sed 's/.fls//' | awk -F'_' '{print $NF}'`
    ftype=`echo $fname |sed 's/.fls//' | awk -F'_' '{print $1}'`
    echo $tag
    set -x
    ./jobFileMakerBMM5.py $f RunLumiFiles/${ftype}_mc_bkg_${tag}/ 1 2000 $tag
    set +x
    done


