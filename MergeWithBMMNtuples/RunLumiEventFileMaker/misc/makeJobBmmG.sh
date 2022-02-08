#!/usr/bin/bash
#./jogFileMakerBMMG.py <FileSource> <destination> <filesPerJob> <njobs> <tag> 

./jobFileMakerBMMG.py fileLists/bmmg_bph2A.fls  RunLumiFiles/BMMG/bph2A 10 20000 bph2A
./jobFileMakerBMMG.py fileLists/bmmg_bph3A.fls  RunLumiFiles/BMMG/bph3A 10 20000 bph3A
./jobFileMakerBMMG.py fileLists/bmmg_bph4A.fls  RunLumiFiles/BMMG/bph4A 10 20000 bph4A
./jobFileMakerBMMG.py fileLists/bmmg_bph5A.fls  RunLumiFiles/BMMG/bph5A 10 20000 bph5A
./jobFileMakerBMMG.py fileLists/bmmg_bph6A.fls  RunLumiFiles/BMMG/bph6A 10 20000 bph6A
