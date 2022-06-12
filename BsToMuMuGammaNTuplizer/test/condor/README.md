Submiting jobs locally in condor for smaller studies

### Run Lumi Json
Get Lumi-Jason Files  for a datset  from DAS
Set the dataset names in `misc/getJsonFromDasOutputs.py`
```bash
python misc/getJsonFromDasOutputs.py
```

### Running the Ntuplizer over MC / Data files
**proxy needs to set for accesing the files from grid**

The file list have to be given as input to the `misc/condorJobMaker.py`.

Usage
```
./misc/condorJobMaker.py \
        <cmsConfigName> \
        <fileList to run over> \
        <result destination> \
        <max number of jobs> \
        <number of iles per job> \
        <events to be analysed per job> \
        <tag> \
        <max mterialize option for Condor submission (see [1])  >
```

see `misc/jobMaker.sh` for some examples [ shuold be able to run it  out of the box ]


[1] [Submitting Lots Of Jobs in Condor](https://htcondor.readthedocs.io/en/latest/users-manual/submitting-a-job.html#submitting-lots-of-jobs)
