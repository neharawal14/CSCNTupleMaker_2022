# CSCNTupleMaker

### To install
#### This branch has the updated code for processing 2022 data
```
cmsrel CMSSW_12_4_6
cd CMSSW_12_4_6/src/
git init
git remote add origin git@github.com:neharawal14/CSCNTupleMaker_2022.git
git fetch origin
git checkout origin/dev_2022
git checkout -b <branch_name_you_choose> # optional, if you want to make further development it is better to use make a local branch with this command
git submodule init
git submodule update
scram b -j 8
```

#### The submodule points to the branch dev_submodule_UFCSCSoftware. You can change the branch to which submodule points for your own submodule development using
```
git remote set-url origin new-remote-url
```

## For submitting jobs
1. Write dataset name in file '2016_dataset.txt'
2. Modify the Globla tag in the python file 'UFCSCSoftware/UFCSCRootMaker/test/UFCSCRootMaker_template.py'
3. Modify 'SubmicrabJobs.py' to include dataset name
4. Modify 'crabConfigTemplate.py' to include the lumi mask to produce the dataset
5. Submit using the script - makeNTuple_2016.sh
