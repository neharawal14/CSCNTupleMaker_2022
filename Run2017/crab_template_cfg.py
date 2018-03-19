from WMCore.Configuration import Configuration

version = 'v1'
runEra = 'Run2017E'

#____________________________________________________________||
config = Configuration()
config.section_("General")

#____________________________________________________________||
config.General.requestName = 'CSCNTuple'+version+'-SingleMuon'+runEra
config.General.transferLogs=True

#____________________________________________________________||
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'UFCSCRootMaker.py'
#config.JobType.scriptExe = 'crab_script.sh'
#config.JobType.inputFiles = ['crab_script.py','haddnano.py','keep_and_drop.txt',] #hadd nano will not be needed once nano tools are in cmssw
#config.JobType.sendPythonFolder  = True

#____________________________________________________________||
config.section_("Data")
config.Data.inputDataset = '/SingleMuon/'+runEra+'-PromptReco-v1/RECO'
config.Data.inputDBS = 'global'
config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.splitting = 'EventBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 2000
config.Data.outLFNDirBase = '/store/user/klo/CSC/CSCNTuple/SingleMuon'+runEra+'-'+version+'/'
config.Data.publication = False
config.Data.outputDatasetTag = runEra+'-PromptReco-'+version

#____________________________________________________________||
config.section_("Site")
config.Site.storageSite = "T2_US_Florida"
