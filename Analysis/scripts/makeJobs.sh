#Usage
#    ./makeCondorJobForAnalysis.py <InputFileListFname> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
#           InputFileListFname : path to a file which contains the paths to files to be analysed , one per line [ xrootd adresses will work ]
#           destination : place where final files have to be copied
#

./makeCondorJobForAnalysis.py \
        charmD_mergedFileList.txt \
        /eos/cms/store/group/phys_bphys/athachay/bs2mmg/AnalysisNtuples/Data/Charmoniun/2018D/ \
        10 \
        40 \
        Charm18D
