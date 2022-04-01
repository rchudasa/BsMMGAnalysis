#	./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>
./makeCondorJobForAnalysis.py \
    analysis.exe \
    fileList/bph6A.fls \
    configs/analysis.tpl.cfg \
    1 \
    results/bParkAnalysis/bph6A  \
    1000 \
    1  \
    bph6A
