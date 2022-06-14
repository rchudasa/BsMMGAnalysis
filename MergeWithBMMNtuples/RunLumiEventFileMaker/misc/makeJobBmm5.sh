#!/usr/bin/bash
#./jogFileMakerBMMG.py <FileSource> <destination> <filesPerJob> <njobs> <tag> 

./jobFileMakerBMM5.py fileLists/bmm5_bph2A.fls  RunLumiFiles/BMM5/bph2A 10 20000 bph2A
./jobFileMakerBMM5.py fileLists/bmm5_bph3A.fls  RunLumiFiles/BMM5/bph3A 10 20000 bph3A
./jobFileMakerBMM5.py fileLists/bmm5_bph4A.fls  RunLumiFiles/BMM5/bph4A 10 20000 bph4A
./jobFileMakerBMM5.py fileLists/bmm5_bph5A.fls  RunLumiFiles/BMM5/bph5A 10 20000 bph5A
./jobFileMakerBMM5.py fileLists/bmm5_bph6A.fls  RunLumiFiles/BMM5/bph6A 10 20000 bph6A
./jobFileMakerBMM5.py fileLists/bmm5_bph6B.fls  RunLumiFiles/BMM5/bph6B 10 20000 bph6B
