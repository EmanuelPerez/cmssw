# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("L1TkCalo")

#
# This configuration runs over a file that contains the L1Tracks
# and the tracker digis.
# It creates HT and MHT, from the L1Jets and from the HLT HI jets.
#
 

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
   'file:example_w_Tracks_and_vertex.root'	# created by runL1Tracks_and_L1PrimaryVertexProducer_cfg.py
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


# ---------------------------------------------------------------------------

# ---  Run the SLHCCaloSequence  to produce the L1Jets

process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.L1CalibFilterTowerJetProducer.pTCalibrationThreshold = cms.double(40) # applies calibration only to > 40GeV L1 jets

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

#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#
# --- Load the L1TkCaloSequence :

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkCaloSequence_cff")

# -- Produce L1TkJets, HT and MHT from the L1Jets :
process.L1TkCaloL1Jets = cms.Path( process.L1TkCaloSequence )

# -- Produce the HLT JI Jets and L1TkJets, HT and MHT from  these jets :
process.L1TkCaloHIJets = cms.Path( process.L1TkCaloSequenceHI )

#
# ---------------------------------------------------------------------------




process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TkCalo.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append( 'keep *_*_*_L1EG' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')
process.Out.outputCommands.append('keep *_L1TkHTMiss*_*_*')
process.Out.outputCommands.append('keep *_L1TkJets*_*_*')
process.Out.outputCommands.append('keep *_iterativeConePu5CaloJets_*_*')
#process.Out.outputCommands.append('keep *_L1TowerJetPUSubtractedProducer_PUSubCen8x8_*')
process.Out.outputCommands.append('keep *_L1CalibFilterTowerJetProducer_CalibratedTowerJets_*')
process.Out.outputCommands.append('keep *_L1CalibFilterTowerJetProducer_UncalibratedTowerJets_*')

process.Out.outputCommands.append('keep *_*_*_L1TkCalo')


#process.Out.outputCommands.append('keep *')

process.FEVToutput_step = cms.EndPath(process.Out)



