# Merging BMMG Ntuples with the BMM5 Ntuples

- For Datas the merging process foillows 3 steps :
  - Step 1 : Make the BMMG and BMM5 Ntuples of the smaple of interest.
  - Step 2 : Processing of BMM5 ntuples to make a list of events/lumis/runs stored in them.
  - Step 3 : Processing od BMMG ntuples to make a list of events/lumis/runs stored in them.
  - Step 4 : Submitting the condor jobs for merging the ntuples

## Workflow 
#### Processing of BMM5 / BMMG ntuples to get events/lumis/runs list.
- `BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/getEventListFromBMM5.cc` script does this job
- Need the BMM5 tree interface compiled beforehand (_use Bmm5Ntuple.h and Bmm5Ntuple.C_ ) 
- Usage 
```
$ root -b -q 'getEventListFromBMM5.cc("bmm5Selection.cfg")'
$ root -b -q 'getEventListFromBMMG.cc("bmmGSelection.cfg")'
```
this step will produce a *.txt and *.root file
- Use scripts `jobFileMakerBMM5.py ` and `jobFileMakerBMMG.py` to make jobs in condor
  ```
  # ./jogFileMakerBMMX.py <FileSource> <destination> <filesPerJob> <tag> 
  ./jobFileMakerBMMX_test.py bmm5_mc_bkg_hPMuNu_BsToKMuNu_v2.files \
          /eos/home-a/athachay/workarea/bs2mumug/run2studies/ntuplizer/RunLumiFileMaker/MC/bmm5_mc_bkg_hPMuNu_BsToKMuNu_v2\
          5 \
          100 \
          BsToKMuNu
  condor_submit subCondorBMMXBsToKMuNu.sub
  ```

#### Merging the Ntuples
- Update the `Bmm5Ntuple.h` and `BmmGNtuple.h` with relavant headers of the Ntuples of interest.
- Update the `branchList` file wit the branches of interest.
    - I usually copy the branch list from BmmXNtuple.h and edit in vi to the format described below. 
- Customize the `genericTreeBranchSelector` class to have the proper data members and functions.Use `makeTheClassTemplate.py` script to make the customizations. See usage as :
  - Fill the `branchList` file with template
     ```
     className, datatype, BranchName,arrayVar
     ```
     - classname : tag to the source tree
     - datatype  : the datatype of the branch
     - BranchName : Name of the branch in original tree
     - arrayVar   : length of the array in the event [ not mandatory for non array branches]
     - Note : vector<DTYPE> branches will be converted to array type branches automatically
  - Then execute the `makeTheClassTemplate.py`  to produce `Functoion.cc` file.
  ```
  ./makeTheClassTemplate.py
  ```
  - Copy the Data Member section from `Functoion.cc` to the class defenition in `genericTreeBranchSelector.h`
  - Copy  these functions from `Functoion.cc` to `genericTreeBranchSelector.C`
    - `void genericTreeBranchSelector::FillTheArraysFromVectors()`
    - `void genericTreeBranchSelector::setCompiledTreeBranches()`
- Copy the  file produced using `getEventListFromBMMX.cc` [if its condor submission, merge all txt outputs to a single file] to `
BsMMGAnalysis/MergeWithBMMNtuples/MergeFiles/
`
- Make the  condor merger jobs
  ```
  ./makeTheMergerCondorJobs.py 
  condor_submit subMergerJobs.sub
  ```
- Done !!
