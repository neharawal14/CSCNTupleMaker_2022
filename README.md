# CSCNTupleMaker

### To install
```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src/
git init
git remote add origin git@github.com:UF-CSC/CSCNTupleMaker.git
git fetch origin
git checkout origin/master
git checkout -b <branch_name_you_choose> # optional, if you want to make further development it is better to use make a local branch with this command
git submodule init
git submodule update
scram b -j 8
```
PS : Updated code for 2022 is in dev_2022 branch
