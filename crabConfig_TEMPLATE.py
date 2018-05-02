from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = Configuration()

config.section_('General')
config.General.workArea = 'resultsAna_JOBTAG/'
config.General.requestName = 'OUTFILENAME'
config.General.transferOutputs = True
config.General.transferLogs=True
config.General.failureLimit=1

config.section_('JobType')
config.JobType.scriptExe = 'submitFileCrab.sh'
#config.JobType.inputFiles = ['/scratch/osghpc/dsperka/Run2/HZZ4l/CMSSW_7_6_3_patch2/src/ZZMatrixElement/MEKD','/scratch/osghpc/dsperka/Run2/HZZ4l/CMSSW_7_6_3_patch2/src/KinZfitter/KinZfitter/ParamZ1','/scratch/osghpc/dsperka/Run2/HZZ4l/CMSSW_7_6_3_patch2/src/KinZfitter/HelperFunction/hists']
config.JobType.psetName = 'CFGFILE'
config.JobType.pluginName = 'Analysis'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['OUTFILENAME.root']

config.section_('Data')
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
#config.Data.inputDBS = 'phys03'
config.Data.inputDataset = 'DATASETNAME'
if ('Run2017' in 'DATASETNAME'):
  config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 5
else:
  config.Data.splitting = 'FileBased'
  config.Data.unitsPerJob = 5
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/%s/rootfiles_2017/JOBTAG/' % (getUsernameFromSiteDB())
#config.Data.outputDatasetTag = 'CSCLocalReco_OUTFILENAME'
config.Data.ignoreLocality = False
config.Data.allowNonValidInputDataset = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_Florida'

