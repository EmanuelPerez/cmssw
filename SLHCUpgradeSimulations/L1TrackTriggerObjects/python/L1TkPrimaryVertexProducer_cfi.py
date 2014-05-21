import FWCore.ParameterSet.Config as cms

L1TrackPrimaryVertex = cms.EDProducer('L1TkFastVertexProducer',
     L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),
     PTMINTRA = cms.double( 2.),        # PTMIN of L1Tracks, in GeV
     nStubsmin = cms.int32( 4 ) ,       # minimum number of stubs
     nStubsPSmin = cms.int32( 3 ),       # minimum number of stubs in PS modules 
     PTSAT = cms.double( 50. )		# in GeV. Tracks with PT above PTSAT are considered as mismeasured
					# and their PT is set to the saturation value PTSAT.
)
