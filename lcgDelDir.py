import os, sys, fnmatch

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

#dirName = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_Feb20/SingleElectron/crab_SingleElectron_Run2015C_50ns-16Dec2015-v1/'
#dirName = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_Feb20/SingleMuon/crab_SingleMuon_Run2015C_50ns-16Dec2015-v1/'
#dirName = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_Feb20/DoubleMuon/crab_DoubleMuon_Run2015C_50ns-16Dec2015-v1/'
#dirName = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_Feb20/DoubleEG/crab_DoubleEG_Run2015C_50ns-16Dec2015-v1/'
#dirName = '/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_Feb20/MuonEG/crab_MuonEG_Run2015C_50ns-16Dec2015-v1/'
dirName = '/cms/data/store/user/hmei/rootfiles_2017/2016SingleMu_BCDEFGH_promptReco/'

for filename in find_files(dirName, '*'):
  print 'Found C source:', filename
  os.system('lcg-del -b -l -D srmv2 srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN='+str(filename))

# to copy
#lcg-cp -v -b -D srmv2 file:/test   srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN=/cms/data/store/user/dsperka/

subdirs = [x[0] for x in os.walk(dirName)]

for i in range(len(subdirs)):
  subdir = subdirs[len(subdirs)-1-i]
  print subdir
  os.system('lcg-del -b -l -d -D srmv2 srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN='+subdir)
