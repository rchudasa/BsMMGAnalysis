# Analysis 

## Helper scripts

#### bmmX_distributionStudy.cc
- Makes an ntuple with hists of branches described in *.cfg files, see _charmoniumAnalysisConfig.cfg_  for an example cfg file
- the root file will also have a root tree with all the  BMMG Candidates and corresponding branches from the merged ntuples
  - the branches are described in assignXX.cc and setupXX.cc files
- For making the "assignXX.cc" and "setupXX.cc" files use  `makeSetupAndFillCodeFromBranches.py`
  - this takes in branches.br.XX files to setup branches
- Assumes the MergedBMMX.h and MergedBMMX.C are alredy made and compiled ready.
- Usage 
  ```
  $  root -b -q -l 'bmmX_distributionStudy.cc("charmoniumAnalysisConfig.cfg")'
  $  root -b -q -l 'bmmX_distributionStudy.cc("mcRun3AnalysisConfig.cfg")'
  ```
## Making Analysis Ntuples
  Selections Implemented :
    - Triggers of Interest
    - number of PV ( > 0 )
    - golden json based run-lumi masking
    Note : other selection criterias in the *.cfg files are not for used for making the ntuple
  #### Setting up condor jobs :
    #clone repo
    cd BsMMGAnalysis
    git checkout analysis
    cd Analysis
    make  analysisNtupleMaker
    # edit the makeCondorJobForAnalysis.py to set `executable='analysisNtupleMaker'` 
    cp scripts/makeJobs.sh .
    # Edit the makeJobs. to reflect ur requirement
    ./makeJobs.sh
    # Now JobsXX folder will be made with for the condor submission scripts
    condor_submit subjobXX.sub

 After all the jobs end [ or are past the IDLE state in condor ] one can do `condor_submit subjobXX.sub` again to resubmit the failed jobs.
 #### example , for making analysis ntuples for charmonium 2018 D
 Getting Files to run on
     
     realpath /eos/cms/store/group/phys_bphys/athachay/bs2mmg/MergedNtuples/Data/Charmoniun/2018D/*.root >  charmD_mergedFileList.txt
 
 Set makeJobs.sh as 
 
     ./makeCondorJobForAnalysis.py \
        charmD_mergedFileList.txt \
        /eos/cms/store/group/phys_bphys/athachay/bs2mmg/AnalysisNtuples/Data/Charmoniun/2018D/ \
        10000 \
        40\
        Charm18D
        
  (40 is a normal number u can keep for merging .. can try our other options as well )
    
  
