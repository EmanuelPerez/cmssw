import FWCore.ParameterSet.Config as cms

process = cms.Process("TRA")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.singleElectronFiles_cfi import *
from SLHCUpgradeSimulations.L1TrackTriggerObjects.ttbarFiles_p1_cfi import *


process.source = cms.Source("PoolSource",
    #fileNames = singleElectronFiles
    fileNames = ttbarFiles_p1
)


# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
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

# --- now we runn the L1Track producer :

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')

process.pL1Tracks = cms.Path( process.BeamSpotFromSim*process.L1Tracks )

#
# ---------------------------------------------------------------------------


process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_w_Tracks.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_TRA' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_*gen*_*_*')
process.Out.outputCommands.append('keep *_*Gen*_*_*')
process.Out.outputCommands.append('keep *_rawDataCollector_*_*')
process.Out.outputCommands.append('keep *_L1TkStubsFromPixelDigis_StubsPass_*')

process.FEVToutput_step = cms.EndPath(process.Out)




