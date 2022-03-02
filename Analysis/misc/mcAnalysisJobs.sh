# ./makeTheMergerCondorJobs.py < exe > < cfg template  > <File source> <mother files list> <son file list> <destination> <tag> <njobs> <bkpDestination>
#./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
NJOBS=${1-5}
FILES_PER_JOB=${2-1}
MAXEVENTS=${3--1}
echo NJOBS : $NJOBS
echo FILES_PER_JOB : $FILES_PER_JOB
echo MAXEVENTS : $MAXEVENTS
echo ""
EXECUTABLE=mcAnalysis.exe

declare -a SourceFiles=(\
# "fileList/mc_sig_BsTMMG_v2.files" \
"fileList/mc_bkg_bdToKK.files" \
"fileList/mc_bkg_bdToKPi.files" \
"fileList/mc_bkg_bdToPiMuNu.files" \
"fileList/mc_bkg_bdToPiPi.files" \
"fileList/mc_bkg_bsToKK.files" \
"fileList/mc_bkg_bsToKMuNu.files" \
"fileList/mc_bkg_bsToKPi.files" \
"fileList/mc_bkg_bsToPiPi.files" \
)

declare -a tagArr=(\
#"mc_sig_bs2mmg" \
"mc_bkg_bdToKKg" \
"mc_bkg_bdToKPig" \
"mc_bkg_bdToPiMuNug" \
"mc_bkg_bdToPiPig" \
"mc_bkg_bsToKKg" \
"mc_bkg_bsToKMuNug" \
"mc_bkg_bsToKPig" \
"mc_bkg_bsToPiPig" \
)

declare -a AnalysisOption=(\
# 0 \
1 \
1 \
1 \
1 \
1 \
1 \
1 \
1 \
)

declare -a CfgTemplate=(\
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
"configs/analysisMC.tpl.cfg" \
)

# ./makeCondorJobForAnalysis.py mvaDataMaker.exe configs/QcdSample.tpl 1 /home/athachay/t3store3/bs2mumug/photonID/analysis/CMSSW_10_6_29/src/BsMMGAnalysis/PhotonID/test/results/MC/mc_bkg_qcd30To50EMEnriched 5 1 -1 mc_bkg_qcd30To50EMEnriched
for i in "${!tagArr[@]}"; do 
    echo $i : ${jobArr[$i]}
    src=${SourceFiles[$i]}
    TAG=${tagArr[$i]}
    ANALYSIS_OPT=${AnalysisOption[$i]}
    CFG_TEMPLATE=${CfgTemplate[$i]}
    set -x
    ./makeCondorJobForAnalysis.py \
        $EXECUTABLE \
        $src \
        $CFG_TEMPLATE \
        $ANALYSIS_OPT \
        $PWD/results/MC/$TAG \
        $NJOBS \
        $FILES_PER_JOB \
        $MAXEVENTS \
        $TAG
    set +x
done
