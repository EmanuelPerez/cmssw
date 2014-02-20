import FWCore.ParameterSet.Config as cms

process = cms.Process("VTX")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.ttbarFiles_p1_cfi import *

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = ttbarFiles_p1
    fileNames = cms.untracked.vstring('/store/group/comm_trigger/L1TrackTrigger/BE5D_620_SLHC6/singleMu/PU140/SingleMuMinus_BE5D_PU140.root')
)


# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# ---- Run the L1Tracking :

# ---- redo the stubs. Stubs were produced during the central production
#      and are present on the DIGI files, but the "z-matching" condition
#      was enforced. Here we redo the stubs without the z-matching.
#      This leads to better tracking efficiencies.

process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
#process.pStubs = cms.Path( process.L1TkStubsFromPixelDigis )

# --- now we runn the L1Track producer :

process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

#process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")

#process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
#process.L1Tracks.geometry = cms.untracked.string('BE5D')

#process.pL1Tracks = cms.Path( process.BeamSpotFromSim*process.L1Tracks )


# --- Run the L1PrimaryVertex producer :

# the vtx is calculated from tracks that have | z | < ZMAX and chi2 < CHI2MAX.
# The vtx maximises e.g. Sum (PT^2)  where the sum runs over tracks that
# are within | z - z_track | < DeltaZ  of the tested vertex.

# The primary vertex producer has to be run together with the L1Tracking,
# or on a file on which the tracker digis have been kept.
# This is because one needs to access the number of stubs in the PS
# modules for the L1Tracks, when one reconstructs the vertex.
# To get this information, one must have access to the tracker digis.

# the configuration parameters below are not yet optimized !

process.L1TrackPrimaryVertex = cms.EDProducer('L1TrackPrimaryVertexProducer',
     #L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
     L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
     ZMAX = cms.double ( 25. ) ,	# in cm
     CHI2MAX = cms.double( 100. ),
     DeltaZ = cms.double( 0.1 ),    	# in cm.   
     PTMINTRA = cms.double( 2.), 	# PTMIN of L1Tracks, in GeV
     nStubsmin = cms.int32( 4 ) ,	# minimum number of stubs
     nStubsPSmin = cms.int32( 3 )       # minimum number of stubs in PS modules 
)

process.p = cms.Path( process.L1TrackPrimaryVertex )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_w_Tracks_and_vertex.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_*gen*_*_*')
process.Out.outputCommands.append('keep *_*Gen*_*_*')
process.Out.outputCommands.append('keep *_rawDataCollector_*_*')
process.Out.outputCommands.append('keep *_L1TkStubsFromPixelDigis_StubsPass_*')

process.FEVToutput_step = cms.EndPath(process.Out)




