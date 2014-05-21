import FWCore.ParameterSet.Config as cms
process = cms.Process("Ele")

process.load("FWCore.MessageService.MessageLogger_cfi")

################################################################################
# Example configuratiom file : 
# Here we run the L1EG algorithms (old stage-2 and new clustering),
# we unpack the L1EG objects that were created during the L1 step
# of the central production (i.e. the Run-1 algorithms), and we
# create L1TkElectron particles (matching with L1Tracks).
# Two output collections are created :
#    - L1TkElectrons
#    - L1TkElectrons that are isolated with respect to the L1Tracks
################################################################################

# list of files
file_names = cms.untracked.vstring(
# '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00D6C34E-0339-E311-836A-002618943880.root',
 '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FEDB1C0F-FF38-E311-A659-0025905938D4.root')

# input Events 
process.source = cms.Source("PoolSource",
   fileNames = file_names,
   skipEvents = cms.untracked.uint32(0) 
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

# ---------------------------------------------------------------------------
#
# ---- Run the L1Tracking :

# ---- redo the stubs. Stubs were produced during the central production
#      and are present on the DIGI files, but the "z-matching" condition
#      was enforced. Here we redo the stubs without the z-matching.
#      This leads to better tracking efficiencies.

process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.pStubs = cms.Path( process.L1TkStubsFromPixelDigis )

# L1Tracking
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')

process.BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.BeamSpotFromSim*process.L1Tracks)


# ---------------------------------------------------------------------------
#
# --- Run the L1EG algorithm of Jean-Baptiste (new clustering algorithm).
# --- Note thet the efficiency is poor at very high PU...
#
# --- This also runs the "old" stage-2 algorithm

process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.pSLHCCalo = cms.Path(
    process.RawToDigi+
    process.SLHCCaloTrigger
)
# bug fix for missing HCAL TPs in MC RAW
process.pSLHCCalo.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode    = cms.bool(True)
process.valRctDigis.hcalDigis         = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

# run L1Reco to produce the L1EG objects corresponding
# to the current trigger
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.L1Reco = cms.Path( process.l1extraParticles )

# ---------------------------------------------------------------------------
#
# --- test L1TkElectronTrack


# "electrons" :
import SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkElectronTrackProducer_cfi
# no Isolation applied
process.L1TkElectrons = SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkElectronTrackProducer_cfi.L1TkElectrons.clone()
process.pElectrons = cms.Path( process.L1TkElectrons )

# Isolated (w.r.t. L1Tracks) electrons :
process.L1TkIsoElectrons = process.L1TkElectrons.clone()
process.L1TkIsoElectrons.IsoCut = cms.double(0.1)
process.pElectronsIso = cms.Path( process.L1TkIsoElectrons)


process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TrackElectron.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_SLHCL1ExtraParticles_EGamma_*' )
process.Out.outputCommands.append( 'keep *_SLHCL1ExtraParticles_IsoEGamma_*' )
process.Out.outputCommands.append( 'keep *_L1TkElectrons_*_*' )
process.Out.outputCommands.append( 'keep *_L1TkIsoElectrons_*_*' )
process.Out.outputCommands.append( 'keep *_genParticles_*_*')

process.Out.outputCommands.append( 'keep SimTracks_g4SimHits_*_*') 
process.Out.outputCommands.append('keep *_L1Tracks_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.pSLHCCalo,process.TT_step,process.pElectrons,process.pElectronsIso, process.FEVToutput_step)




