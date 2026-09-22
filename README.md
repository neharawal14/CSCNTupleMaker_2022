# CSCNTupleMaker

### To install
#### This branch has the updated code for processing RAW-RECo dataset and produce ntuples for CSC longevity study
```
Setup your CMSSW environment ; for Run3 I would suggest use CMSSW_16_1_3
```
cmsrel CMSSW_12_4_6
cd CMSSW_12_4_6/src/
git init
git remote add origin git@github.com:neharawal14/CSCNTupleMaker_2022.git
git fetch origin
git checkout -b new_branch origin/new_branch # optional, if you want to make further development it is better to use a create your own local branch with a new name
git submodule init
git submodule update
scram b -j 8
```

#### The submodule points to the branch dev_submodule_UFCSCSoftware. You can change the branch to which submodule points for your own submodule development using
```
git remote set-url origin new-remote-url
```

Important part of code is in : UFCSCSoftware/UFCSCRootMaker/
The file src/UFCSCRootMaker.cc : macro to generate ntuple for longevity study
The file test/UFCSCRootMaker_template.py : 
(1) template with how to process and execute the macro. To process each dataset separately, the important thing that will change in this file is "GlobalTag" in this file. 
(2) Execute command "cmsRun UFCSCSoftware/UFCSCRootMaker/test/UFCSCRootMaker_template.py" to run the code src/UFCSCRootMaker.cc over a single input root file. 
While to run over whole dataset, one need to submit crab jobs.

## For submitting jobs : an example case of 2016: 
1. Write dataset name in file '2016_dataset.txt'
2. Modify the Global tag in the python file 'UFCSCSoftware/UFCSCRootMaker/test/UFCSCRootMaker_template.py'
3. Modify 'SubmicrabJobs.py' to include dataset name
4. Modify 'crabConfigTemplate.py' to include the lumi mask to produce the dataset
5. Make sure to have correct golden json file to process events : eg: "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
5. Submit using the script - makeNTuple_2016.sh
