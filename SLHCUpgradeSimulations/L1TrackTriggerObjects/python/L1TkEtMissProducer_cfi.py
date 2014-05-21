import FWCore.ParameterSet.Config as cms

L1TkEtMiss = cms.EDProducer('L1TkEtMissProducer',
     L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
     L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),      
     PTMINTRA = cms.double( 2. ),       # in GeV
     DeltaZ = cms.double( 1.0 ),       # in cm
     nStubsmin = cms.int32( 4 )
)
