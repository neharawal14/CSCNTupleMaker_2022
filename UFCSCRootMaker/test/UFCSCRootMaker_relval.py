######################################################################
#                                                                    #
# This version works with CMSSW_6_2_X                                #
#                                                                    #
######################################################################
import FWCore.ParameterSet.Config as cms

########## Options ############
isDATA = bool(False)
isRAW = bool(True)
isDIGI = bool(True)
isSIM = bool(True)
isGEN = bool(True)
isLocalRECO = bool(True)
isFullRECO = bool(True)

addMuonInfo = bool(True)
addTrackInfo = bool(True)
addRecHitInfo = bool(True)
addSegmentInfo = bool(True)
addTriggerInfo = bool(True)
addDigiInfo = bool(True)
addTimeMonitoringInfo = bool(True)
addCalibrationInfo = bool(True)

maxEvents = -1

MCGlobalTag = "POSTLS161_V15::All"
#MCGlobalTag = "MC_61_V11::All"
#MCGlobalTag = "STAR17_61_V1A::All"
#MCGlobalTag = "DES17_61_V5::All"
DataGlobalTag = "FT_53_V21_AN3::All"

doDebug = bool(False)
###############################

### Debug Printing ###
if not isDATA :
    print "Sample Type: MC"
else :
    print "Sample Type: Data"
                
#####################
process = cms.Process("UFCSCRootMaker")

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.suppressWarning.append('classByHitsGlb') # kill stupid RPC hit associator warning
#process.MessageLogger.cerr.FwkJob.limit=1
#process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(1))
                                                       
######################### Frontier Conditions #########################
# Conditions data for calibration and alignment                       #
# are defined in the Offline Conditions Database (ORCOF),             #
# which is read in CMSSW applications via Frontier caching servers.   #
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions  #
#######################################################################
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')


if not isDATA:
    process.GlobalTag.globaltag=MCGlobalTag 
    
else:
    process.GlobalTag.globaltag=DataGlobalTag


process.out = cms.OutputModule("PoolOutputModule",
                               #fileName = cms.untracked.string('test.root'),
                               # save only events passing the full path
                               outputCommands = cms.untracked.vstring('drop *')
                               )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("/cms/data/store/user/patTuples/CSC/mc/CMSSW_6_1_2/GEN-SIM-RAW-DIGI-RECO/VBF_HToZZTo4L_14TeV_pythia6_PU35_bx50.root")
                                   )

# Primary Vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
#process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
#                                                  src = cms.InputTag('offlinePrimaryVertices'),
#                                                  cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
#                                                  filter = cms.bool(True)
#                                                  )


process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             # processingMode = cms.untracked.string('RunsAndLumis'),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = cms.untracked.vstring(),      
                             )

if isDATA:
    process.source.fileNames = cms.untracked.vstring(
#        'file:/cms/data/store/data/Run2012B/SingleMu/RECO/22Jan2013-v1/20002/00351421-C271-E211-A5C9-90E6BA19A20A.root',
#        'file:/cms/data/store/data/Run2012B/SingleMu/RECO/22Jan2013-v1/20040/F04E9749-2D74-E211-9EE4-E0CB4E19F9BC.root',
#        'file:/cms/data/store/data/Run2012B/SingleMu/RECO/22Jan2013-v1/20040/F2AF82F1-2C74-E211-9CD1-001EC9D80781.root'

    
        )
else:
    process.source.fileNames = cms.untracked.vstring()
    process.source.fileNames.extend( [
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/0250F3F1-3E4B-E211-A9DB-0026189437FE.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/0289C438-874D-E211-9FB2-002618943975.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/1A442014-434B-E211-9319-003048678FB8.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/1C205C60-BB4E-E211-89BA-002618943948.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/2C1245D4-B955-E211-B9DC-003048FFD744.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/2E473814-234B-E211-95D0-003048678FA6.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/34561B73-3F4B-E211-AF53-00261894380B.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/3C5CEFE8-374B-E211-BADB-003048FFD752.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/42209C76-244B-E211-8CFA-002618943877.root',
        '/store/mc/Summer12/VBF_HToZZTo4L_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU35_POSTLS161_V12-v1/10000/46D7C4F6-AB4B-E211-B8A9-003048679180.root'
    ]
    )
    
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False)
                                     #,SkipEvent = cms.untracked.vstring('ProductNotFound')
                                     )

### Physics Declared Filter (for data)
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

### No scraping
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(doDebug),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.2)
                                  )


process.LumiCorrectionSource = cms.ESSource("LumiCorrectionSource",
                                            #authpath=cms.untracked.string('/afs/cern.ch'file:/cms/lumi/DB'),
                                            #connect=cms.string('oracle:/'file:/cms_orcon_adg'file:/cms_lumi_prod')
                                            connect=cms.string('frontier://LumiCalc_LUMI_PROD')
                                            #normtag=cms.untracked.string('HFV2a')
                                            #datatag=cms.untracked.string('v3')
                                            )


# HB + HE noise filtering
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE        = 999999.
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = 999999
process.HBHENoiseFilter.minIsolatedNoiseSumEt       = 999999.

#EventCount                                                                                                                             
process.nEventsTotal = cms.EDProducer("EventCountProducer")

#RM
process.load("UFCSCSoftware.UFCSCRootMaker.cscRootMaker_cfi")
process.cscRootMaker.isFullRECO = cms.untracked.bool(isFullRECO)
process.cscRootMaker.isLocalRECO = cms.untracked.bool(isLocalRECO)
process.cscRootMaker.isGEN = cms.untracked.bool(isGEN)
process.cscRootMaker.isSIM = cms.untracked.bool(isSIM)
process.cscRootMaker.isRAW = cms.untracked.bool(isRAW)
process.cscRootMaker.isDIGI = cms.untracked.bool(isDIGI)
process.cscRootMaker.isDATA = cms.untracked.bool(isDATA)
process.cscRootMaker.addMuons = cms.untracked.bool(addMuonInfo)
process.cscRootMaker.addTracks = cms.untracked.bool(addTrackInfo)
process.cscRootMaker.addRecHits = cms.untracked.bool(addRecHitInfo)
process.cscRootMaker.addSegments = cms.untracked.bool(addSegmentInfo)
process.cscRootMaker.addTrigger = cms.untracked.bool(addTriggerInfo)
process.cscRootMaker.addDigis = cms.untracked.bool(addDigiInfo)
process.cscRootMaker.addTimeMonitoring = cms.untracked.bool(addTimeMonitoringInfo)
process.cscRootMaker.addCalibrations = cms.untracked.bool(addCalibrationInfo)

process.cscRootMaker.hltTagSrc = cms.untracked.InputTag('TriggerResults::RECO')
process.cscRootMaker.clctDigiTagSrc = cms.untracked.InputTag('simCscTriggerPrimitiveDigis::DIGI2RAW')
process.cscRootMaker.alctDigiTagSrc = cms.untracked.InputTag('simCscTriggerPrimitiveDigis::DIGI2RAW')
process.cscRootMaker.corrlctDigiTagSrc = cms.untracked.InputTag('simCscTriggerPrimitiveDigis::DIGI2RAW')

process.p = cms.Path(
    process.nEventsTotal
#    *process.goodOfflinePrimaryVertices
#    *process.noscraping
    *process.cscRootMaker
    )

if not isDATA:
    process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")
    process.totalKinematicsFilter.tolerance=5.0
    process.p.replace(process.nEventsTotal,process.nEventsTotal*process.totalKinematicsFilter)

if isRAW:
    process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.p.replace(process.cscRootMaker,process.muonCSCDigis*process.gtDigis*process.cscRootMaker)

if isDATA:
    process.p.replace(process.nEventsTotal,process.hltPhysicsDeclared*process.nEventsTotal)

