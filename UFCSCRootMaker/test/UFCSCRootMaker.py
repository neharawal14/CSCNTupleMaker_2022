######################################################################
#                                                                    #
# This version works with CMSSW_6_2_X                                #
#                                                                    #
######################################################################
import FWCore.ParameterSet.Config as cms

########## Options ############
isDATA = bool(False)
isRAW = bool(False)
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

#MCGlobalTag = "POSTLS161_V15::All"
MCGlobalTag = "START70_V6::All"
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
process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
#process.load('Geometry.CommonDetUnit.globalTrackingGeometry_cfi')
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
                                   fileName = cms.string("/cms/data/store/user/patTuples/CSC/mc/CMSSW_7_0_0/TTbar/GEN-SIM-DIGI-RECO/TTbar.root")
                                   )

# Primary Vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                                  src = cms.InputTag('offlinePrimaryVertices'),
                                                  cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                                  filter = cms.bool(True)
                                                  )


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
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/00743452-B498-E311-AD84-02163E00EAC9.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/026226E0-9E98-E311-8853-02163E00EAE4.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/06393950-A098-E311-8B37-02163E00E9E5.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/06BE841B-A098-E311-988C-02163E009DD5.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/081FC825-AE98-E311-8A43-02163E00E5BD.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/0A5AE4A3-9F98-E311-AD61-02163E009E7E.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/0CDD688F-A098-E311-8B2F-02163E00E83C.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/0E9553D3-9E98-E311-820E-02163E00EABA.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/0EC738A7-9F98-E311-A370-02163E00E74B.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/12CB8A61-9F98-E311-AE6F-0025904B5B7C.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/1CD8DE68-A098-E311-8372-02163E00E801.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/2803B46B-A098-E311-9212-02163E008CCF.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/28EE6FEC-A098-E311-8044-02163E00E7B3.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/2CAB4C22-A098-E311-8012-02163E00E5CB.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/2CBAC741-A198-E311-A39E-02163E00EA40.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/2ED2B3AF-B498-E311-AE43-02163E00E8B1.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/302C70C1-9E98-E311-A0B8-02163E00EAAE.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/30BED60A-AD98-E311-8614-0030489452D6.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/349443AD-A098-E311-990D-02163E008D93.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/36D93354-9F98-E311-BA86-02163E009ECB.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/3C6C264A-A098-E311-B546-02163E00E97D.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/42107A47-B498-E311-8532-02163E00E945.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/467EAE87-A098-E311-9D1D-02163E008CDC.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/488A31DF-9E98-E311-A65B-02163E00E8A2.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/4A5C92E5-9F98-E311-82B8-02163E008F58.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5220A329-A298-E311-9ACD-02163E00EAEA.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/54872428-A098-E311-9725-02163E00E5ED.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/549D5FDB-A098-E311-9EB9-02163E00E9E4.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5664999C-A098-E311-8147-02163E00A12C.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/56C7A969-B098-E311-9263-02163E007A05.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/58276759-A198-E311-B161-02163E00E9B3.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/58B937D0-9F98-E311-AA92-02163E00E751.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/58CE04A8-9F98-E311-89ED-02163E008F7A.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5A62251E-B498-E311-9B0F-0025904B2C5A.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5AE7F332-A198-E311-9448-02163E00E7E5.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5CA887ED-9E98-E311-9C57-02163E00E902.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/5E99838E-AA98-E311-A8BB-02163E00EAAB.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/601517C5-9E98-E311-8E99-02163E00E797.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/60F568CA-9E98-E311-9E0E-02163E00EA1B.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/6232B1FC-9E98-E311-8387-02163E00EB0D.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/70EA2674-A098-E311-B063-02163E00E757.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/70F19E7D-9F98-E311-AD75-02163E00EA93.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/74C788BD-C398-E311-82D9-003048D373C6.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/7601CF0E-B498-E311-BB40-02163E00EACA.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/76651B78-9F98-E311-82CA-0025904B0F42.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/7831A06E-A198-E311-85E7-002590496FE6.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/7A75D8AA-B498-E311-99E8-00215AEDFE6A.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/807C723B-A098-E311-9507-02163E00E936.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/82E8070B-A098-E311-89A3-02163E009DCB.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/82FF3480-A098-E311-9DAB-02163E00EA42.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/84D88B02-A198-E311-BEAE-02163E00EA77.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/8620EEC1-9E98-E311-8C31-02163E00E8C8.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/88024B88-A198-E311-8F9A-02163E00E99B.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/8A03BB4E-9F98-E311-A6BC-02163E00E9A7.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/8A7661A2-A098-E311-A0C6-02163E00E7EA.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/8CE08D2D-B498-E311-B268-0025904B1FBE.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/90E03631-9F98-E311-8FF0-02163E00E904.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/9253F6BF-9F98-E311-8FC2-02163E008D90.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/92576CD6-A098-E311-9F8F-02163E00E5E1.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/96FD8653-B498-E311-BD8A-02163E00E64F.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/9880027A-A198-E311-84AF-02163E00E9EB.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/9A47414A-A198-E311-A4EA-02163E00E783.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/9ED1D793-A198-E311-A5FD-02163E008F84.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/A49F1D5D-B498-E311-962F-0025904B57A2.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/ACB618BC-9F98-E311-AF90-02163E00EA15.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/B0898749-9F98-E311-AE5C-02163E00E62A.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/B08ABB36-A098-E311-84F6-02163E00A104.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/B0940E0D-9F98-E311-B2B5-02163E00E847.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/B44B8859-A098-E311-8D5D-02163E00EA6E.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/B8F63623-A198-E311-9CBB-02163E00E9A0.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/BA45D396-9F98-E311-8A92-02163E00EACE.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/BC23BF34-9F98-E311-883F-02163E00E74C.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/BE52745B-A198-E311-8694-02163E009DB4.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/C21136C3-A098-E311-AF8C-02163E008DB8.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/C609E842-A098-E311-8ECD-0025904A90CA.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/C60DA802-A498-E311-A2F6-02163E00EB46.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/CEEB1C5E-A298-E311-B29E-02163E008D94.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/D0E66E98-A098-E311-9AF6-02163E009F82.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/D67BB079-A198-E311-8E3D-02163E00E7AD.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/DA48C77F-9F98-E311-839C-02163E00E870.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/DE059063-B498-E311-9E04-02163E00E7F0.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E0814BE8-A598-E311-A0D6-02163E008EBE.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E0ECE688-9F98-E311-86C0-02163E00E802.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E2F28418-A198-E311-AAF0-02163E00E715.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E44FAB93-B498-E311-A185-003048FEC15C.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E6AA16F3-A098-E311-A4ED-02163E008DB5.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E6F7FD6F-9F98-E311-9C40-003048945560.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/E8EF37CC-A098-E311-B7DF-02163E00E945.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/EC754A3B-9F98-E311-A9BB-02163E008D9D.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/F005B45A-A098-E311-A857-0025904949BE.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/F4FF87B5-A098-E311-BF56-02163E00E736.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/F864F812-A198-E311-92B4-02163E00CCC6.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/FA551A46-9F98-E311-9F06-02163E00CFE2.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/FA8B290A-A098-E311-BA45-002590494DA0.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/FE43B2B9-A198-E311-8938-0025904B2B20.root',
        '/store/relval/CMSSW_7_0_0/RelValTTbar/GEN-SIM-DIGI-RECO/START70_V6_FastSim-v2/00000/FE5B07FE-B398-E311-8C66-02163E00EA95.root'
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

if isDATA:
    process.p.replace(process.nEventsTotal,process.hltPhysicsDeclared*process.nEventsTotal)

