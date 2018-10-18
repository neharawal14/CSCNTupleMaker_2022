######################################################################
#                                                                    #
# This version works with CMSSW_6_2_X                                #
#                                                                    #
######################################################################
import FWCore.ParameterSet.Config as cms
doUnpacking = bool(False)

########## Options ############
isDATA = bool(True)
isRAW = bool(False)
isDIGI = bool(False)

isSIM = bool(False)
isGEN = bool(False)

isLocalRECO = bool(True)
isFullRECO = bool(True)

addMuonInfo = bool(True)
addTrackInfo = bool(False)
addRecHitInfo = bool(True)
addSegmentInfo = bool(True)
addTriggerInfo = bool(False)

addDigiInfo = bool(False)

addTimeMonitoringInfo = bool(True)
addCalibrationInfo = bool(False)

maxEvents = 200

DataGlobalTag='92X_dataRun2_Prompt_v11'
doDebug = bool(False)

if not isDATA :
    print "Sample Type: MC"
else :
    print "Sample Type: Data"
                
process = cms.Process("UFCSCRootMaker")

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.suppressWarning.append('classByHitsGlb') # kill stupid RPC hit associator warning
                                                    

process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

if not isDATA:
    process.GlobalTag.globaltag=MCGlobalTag 
    
else:
    process.GlobalTag.globaltag=DataGlobalTag

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *')
                               )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("CSCNTuple.root")
                                   )

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                                  src = cms.InputTag('offlinePrimaryVertices'),
                                                  cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                                  filter = cms.bool(True)
                                                  )


process.source = cms.Source ("PoolSource",
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = cms.untracked.vstring(),      
                             )

if isDATA:
    process.source.fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2017E/SingleMuon/RECO/PromptReco-v1/000/304/778/00001/86F81A31-9BB0-E711-932A-02163E01433C.root'

    )
else:
    process.source.fileNames = cms.untracked.vstring(DUMMYFILELIST)
    process.source.fileNames.extend([
        ]
    )
    
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False)
                                     )

process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(doDebug),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.2)
                                  )
                                                                                                    
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsAfter = cms.EDProducer("EventCountProducer")

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

process.unpack = cms.Sequence(process.muonCSCDigis*process.gtDigis)

process.p = cms.Path(
    process.nEventsTotal
#    *process.goodOfflinePrimaryVertices
#    *process.noscraping
    *process.cscRootMaker
    *process.nEventsAfter
    )

print process.p

if doUnpacking:
    process.p.replace(process.nEventsTotal,process.unpack*process.nEventsTotal)

