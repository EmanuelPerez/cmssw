import FWCore.ParameterSet.Config as cms

process = cms.Process("Ele")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#
# This runs over a file that already contains the L1Tracks.
#
# It creates L1TkEmParticles, starting from the Run-1 like L1EG objects that
# are obtained from the process.L1Reco below.
# If you want to run the stage-2 algorithms instead, see e.g.
# EGamma_FullExample.py.
#


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:example_w_Tracks_and_vertex.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/m1_TTbar_BE5D.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/TTbar_BE5D_97.root'
    )
)


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# --- creates l1extra objects for L1EGamma (here for the Run1 trigger !)

        # raw2digi to get the gct digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )
	# run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )


# --- Now run the L1TkEmParticleProducer 

# "photons" :

process.L1TkPhotons = cms.EDProducer("L1TkEmParticleProducer",
	label = cms.string("EGIsoTrk"),	# labels the collection of L1TkEmParticleProducer that is produced.
                                                # e.g. EG or IsoEG if all objects are kept, or
                                                # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                                # objects that pass a cut RelIso < RelIsoCut are written
                                                # into the new collection.
        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),      # input L1EG collection
						# When the standard sequences are used :
                                                #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
                                                #   - for the "old stage-2" algo (2x2 clustering), use 
                                                #     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
        ETmin = cms.double( -1 ),               # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
	RelativeIsolation = cms.bool( True ),	# default = True. The isolation variable is relative if True,
						# else absolute.
        IsoCut = cms.double( 0.1 ), 		# Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
						# When RelativeIsolation = False, IsoCut is in GeV.
	   # Determination of the isolation w.r.t. L1Tracks :
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	ZMAX = cms.double( 25. ),	# in cm
	CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 2. ),	# in GeV
	DRmin = cms.double( 0.07),
	DRmax = cms.double( 0.30 ),
	PrimaryVtxConstrain = cms.bool( False ),  # default = False
	DeltaZMax = cms.double( 999. ),	   # in cm. Used only when PrimaryVtxConstrain = True
        L1VertexInputTag = cms.InputTag("NotUsed"),	# Used only when PrimaryVtxConstrain = True
)
process.pPhotons = cms.Path( process.L1TkPhotons )


# --- Now run the L1TkElectronProducers : one for the algoritm that
# --- matches with L1Tracks, the other for the algorithm that matches
# --- with stubs.


# "electrons" from L1Tracks :

#process.L1TkElectronsTrack = cms.EDProducer("L1TkElectronTrackProducer",
#        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
#        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
#        label = cms.string("NonIsolated")
#)
#process.pElectronsTrack = cms.Path( process.L1TkElectronsTrack )
#
#
## "electrons" from stubs :
#
#process.L1TkElectronsStubs = cms.EDProducer("L1TkElectronStubsProducer",
#        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
#        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
#        label = cms.string("NonIsolated")
#)
#process.pElectronsStubs = cms.Path( process.L1TkElectronsStubs )



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append('keep *')

process.Out.outputCommands.append('keep *_L1TkPhotons_*_*')
process.Out.outputCommands.append('keep *_L1TkElectronsStubs_*_*')
process.Out.outputCommands.append('keep *_L1TkElectronsTrack_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




