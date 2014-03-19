# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("HIJets")

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.minBiasFiles_p1_cfi import *


process.source = cms.Source("PoolSource",
   #fileNames = minBiasFiles_p1
   fileNames = cms.untracked.vstring(
     #"root://eoscms///store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root"
    '/store/mc/UpgFall13d/HToTauTau_125_14TeV_powheg_pythia6/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FAF4B2D5-0539-E311-B5C3-002618FDA279.root'
   )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    #reportEvery = cms.untracked.int32(500),
    reportEvery = cms.untracked.int32(10),
    limit = cms.untracked.int32(10000000)
)      
       
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# Load geometry
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
                            
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


process.load("Configuration.StandardSequences.Services_cff")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff") ###check this for MC!
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")


process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.p = cms.Path(
    process.RawToDigi+
    process.SLHCCaloTrigger
    )

# bug fix for missing HCAL TPs in MC RAW
process.p.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")


# --- Run the calo local reconstruction
process.towerMaker.hbheInput = cms.InputTag("hbheprereco")
process.towerMakerWithHO.hbheInput = cms.InputTag("hbheprereco")
process.reconstruction_step = cms.Path( process.calolocalreco )


# --- Produce the  HLT HeavyIon jets :
process.load("RecoHI.HiJetAlgos.HiRecoJets_TTI_cff")
process.hireco = cms.Path( process.hiRecoJets )

# --- Put them into "L1Jets"
process.L1JetsFromHIHLTJets = cms.EDProducer("L1JetsFromHIHLTJets",
        ETAMIN = cms.double(0),
        ETAMAX = cms.double(3.),
        HIJetsInputTag = cms.InputTag("iterativeConePu5CaloJets")
)
process.pL1Jets = cms.Path( process.L1JetsFromHIHLTJets )


# --- Produce L1TkJets
process.L1TkJetsHI = process.L1TkJets.clone()
process.L1TkJetsHI.L1CentralJetInputTag = cms.InputTag("L1JetsFromHIHLTJets")
process.L1TkJetsHI.JET_HLTETA = cms.bool(True)
process.pJets = cms.Path( process.L1TkJetsHI )



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "HIJets.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')
process.Out.outputCommands.append('keep *_iterativeConePu5CaloJets_*_*')
process.Out.outputCommands.append('keep *_L1JetsFromHIHLTJets_*_*')
process.Out.outputCommands.append('keep *_L1TkJets_*_*')


process.FEVToutput_step = cms.EndPath(process.Out)





