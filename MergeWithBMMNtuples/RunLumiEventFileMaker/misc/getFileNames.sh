files=`xrdfst2 ls /cms/store/user/athachay/BsToMuMuGamma/Run2Studies/BsToMuMuGammaNtuples/MC/`
for f in $files: 
    do
    f=`echo $f |sed 's/://'`
    tag=`echo $f | awk -F'/' '{print $NF}'`
    echo $tag
    set -x
    `xrdfst2 ls $f | sed 's/\/cms/root:\/\/se01.indiacms.res.in\//' > bmmg\_${tag}.fls`
    set +x
    done
