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
