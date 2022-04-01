# ./makeTheMergerCondorJobs.py < exe > < cfg template  > <File source> <mother files list> <son file list> <destination> <tag> <njobs> <bkpDestination>

EXECUTABLE=analysis.exe
CFG_TEMPLATE=analysis2018Temp.cfg
ANALYSIS_OPT=1

declare -a SourceFiles=(\
"sourceFiles/mc_sig_BsTMMG_v2.files" \
"sourceFiles/mc_bkg_hPh_B0ToKK_v2.files" \
"sourceFiles/mc_bkg_hPh_B0ToKPi_v2.files" \
"sourceFiles/mc_bkg_hPh_B0ToPiPi_v2.files" \
"sourceFiles/mc_bkg_hPh_BsToKK_v2.files" \
"sourceFiles/mc_bkg_hPh_BsToKPi_v2.files" \
"sourceFiles/mc_bkg_hPh_BsToPiPi_v2.files" \
"sourceFiles/mc_bkg_hPMuNu_B0ToPiMuNu_v2.files" \
"sourceFiles/mc_bkg_hPMuNu_BsToKMuNu_v2.files" \
)

declare -a tagArr=(\
"mc_sig_BsToMMG_v2" \
"mc_bkg_hPh_B0ToKK_v2" \
"mc_bkg_hPh_B0ToKPi_v2" \
"mc_bkg_hPh_B0ToPiPi_v2" \
"mc_bkg_hPh_BsToKK_v2" \
"mc_bkg_hPh_BsToKPi_v2" \
"mc_bkg_hPh_BsToPiPi_v2" \
"mc_bkg_hPMuNu_B0ToPiMuNu_v2" \
"mc_bkg_hPMuNu_BsToKMuNu_v2" \
)

declare -a jobArr=(\
"BsToMMG" \
"B0ToKK" \
"B0ToKPi" \
"B0ToPiPi" \
"BsToKK" \
"BsToKPi" \
"BsToPiPi" \
"B0ToPiMuNu" \
"BsToKMuNu" \
)

for i in "${!jobArr[@]}"; do 
    echo $i : ${jobArr[$i]}
    src=${SourceFiles[$i]}
    TAG=${tagArr[$i]}
    job=${jobArr[$i]}
 #   set -x
    ./makeCondorJobForAnalysis.py \
        $EXECUTABLE \
        $src \
        $CFG_TEMPLATE \
        $ANALYSIS_OPT \
        /grid_mnt/t3storage3/athachay/bs2mumug/run2studies/analysis/CMSSW_10_6_4_patch1/src/BsMMGAnalysis/Analysis/results/MC/$TAG \
        200 \
        8 \
        $job 
#    set +x
done
