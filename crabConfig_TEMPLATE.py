from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
#from FWCore.PythonUtilities.LumiList import LumiList

#import config, getUsernameFromSiteDB
config = Configuration()
config.section_('General')
config.General.workArea = 'resultsAna_JOBTAG/'
config.General.requestName = 'OUTFILENAME'
#config.General.requestName = 'OUTFILENAME_missingLumis3'
config.General.transferOutputs = True
config.General.transferLogs=True
config.General.failureLimit=1

config.section_('JobType')
config.JobType.scriptExe = 'submitFileCrab.sh'
config.JobType.psetName = 'CFGFILE'
config.JobType.pluginName = 'Analysis'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['OUTFILENAME.root']

config.section_('Data')
config.Data.inputDataset = 'DATASETNAME'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
if('Run2016' in 'DATASETNAME'):
  #config.Data.lumiMask = 'resultsAna_JOBTAG/crab_OUTFILENAME_missingLumis/results/notFinishedLumis.json'
  config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 4
elif('Run2017' in 'DATASETNAME'):
  #config.Data.lumiMask = 'resultsAna_JOBTAG/crab_OUTFILENAME/results/notFinishedLumis.json'
  config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 4
elif('Run2018' in 'DATASETNAME'):
  #config.Data.lumiMask = 'resultsAna_JOBTAG/crab_OUTFILENAME_missingLumis2/results/notFinishedLumis.json'
  config.Data.lumiMask = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
  config.Data.splitting = 'FileBased'
  #config.Data.splitting = 'Automatic'
  config.Data.unitsPerJob = 4
else:
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 1
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/nrawal/CSC_NTuples_SingleIsoMuonTrigger/JOBTAG/'
#config.Data.outputDatasetTag = 'CSCLocalReco_OUTFILENAME'
config.Data.ignoreLocality = False
config.Data.allowNonValidInputDataset = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_Florida'
