# Merging BMMG Ntuples with the BMM5 Ntuples

For Dats the merging process foillows 3 steps :
    - Step 1 : Processing of BMM5 ntuples to make a list of events/lumis/runs stored in them.
        
## Step 1 

- `BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/getEventListFromBMM5.cc` script does this job
- Need the BMM5 tree interface compiled beforehand (_use Bmm5Ntuple.h and Bmm5Ntuple.C_ ) 
- Usage 
```
$ root -b -q 'getEventListFromBMM5.cc("bmm5Selection.cfg")'
```
    this step will produce a *.txt and *.root file

## Step 2 
-  Convert the event file list produced in step 1 to the format compatable fot `PickEvents` configuration for crab
- use RunLumiEventFileMaker/getRunLumiEventFileForPicEvent.py for doing this 
- Follow [this twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents) for pickevents documentation
