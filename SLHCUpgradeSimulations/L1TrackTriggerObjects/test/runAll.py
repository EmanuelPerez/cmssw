import FWCore.ParameterSet.Config as cms

process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.singleElectronFiles_cfi import *

#
# This runs over a file that already contains the L1Tracks and the L1TkPrimaryVertex.
#
# This runs the L1EG algorithms (stage-2 and new clustering), and 
# creates L1TkEmParticles, L1TkElectrons, L1TrkMET.
#
# It also unpacks the L1Jets (Run-1 algorithm) that have been created
# during the centrakl production. They are used to create L1TkJets and
# L1TkHTMiss - these are just technical templates so far.
#


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:example_w_Tracks_and_vertex.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/zmatchingOff/m1_TTbar_BE5D.root'
    )
)



# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# ---------------------------------------------------------------------------
#
# --- Run the L1TrackEtMiss producer :

# needs the L1TkPrimaryVertex
# TrkMET calculated from all tracks that have |z - z_vertex | < DeltaZ,
# pt > PTMINTRA and at least ( >= ) nStubsmin stubs.
# To check the latter condition, the file must contain the stubs.

process.L1TkEtMiss = cms.EDProducer('L1TkEtMissProducer',
     L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
     L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),
     PTMINTRA = cms.double( 2. ),	# in GeV
     DeltaZ = cms.double( 0.2 ),       # in cm
     nStubsmin = cms.int32( 4 )
)

process.pEtMiss = cms.Path( process.L1TkEtMiss )

#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#
# --- Run the L1EG algorithm of Jean-Baptiste (new clustering algorithm).
# --- Note thet the efficiency is poor at very high PU...
#
# --- This also runs the "old" stage-2 algorithm

process.load('Configuration/StandardSequences/L1HwVal_cff')
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
        
        # run L1Reco to produce the L1EG objects corresponding
        # to the current (Run-1) trigger 
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )


# ---------------------------------------------------------------------------
#
# --- test L1TrackEmParticle


# "photons" :

process.L1TkPhotons = cms.EDProducer("L1TkEmParticleProducer",
        label = cms.string("EGIsoTrk"), # labels the collection of L1TkEmParticleProducer that is produced.
                                                # e.g. EG or IsoEG if all objects are kept, or
                                                # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                                # objects that pass a cut RelIso < RelIsoCut are written
                                                # into the new collection.
        #L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","EGamma"),      # input L1EG collection
	L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
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
        RelativeIsolation = cms.bool( True ),   # default = True. The isolation variable is relative if True,
                                                # else absolute.
        IsoCut = cms.double( 0.1 ),             # Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
                                                # When RelativeIsolation = False, IsoCut is in GeV.
           # Determination of the isolation w.r.t. L1Tracks :
     	L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
        ZMAX = cms.double( 25. ),       # in cm
        CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 2. ),    # in GeV
        DRmin = cms.double( 0.07),
        DRmax = cms.double( 0.30 ),
        PrimaryVtxConstrain = cms.bool( False ),  # default = False
        DeltaZMax = cms.double( 999. ),    # in cm. Used only when PrimaryVtxConstrain = True
        L1VertexInputTag = cms.InputTag("NotUsed"),     # Used only when PrimaryVtxConstrain = True
)
process.pPhotons = cms.Path( process.L1TkPhotons )


# --- Now run the L1TkElectronProducers : one for the algoritm that
# --- matches with L1Tracks, the other for the algorithm that matches
# --- with stubs.


# "electrons" from L1Tracks :

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkElectronTrackProducer_cfi")
process.L1TkElectrons.L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated")

process.pElectronsTrack = cms.Path( process.L1TkElectrons )

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
#

#
# ---------------------------------------------------------------------------

# test the L1TkJetProducer
# The collection of (Run 1) L1Jets have been created above, by unpacking the
# gctDigis and running process.L1Reco.

# --- Now run the L1TkJetProducer 

process.L1TkJets = cms.EDProducer("L1TkJetProducer",
        L1CentralJetInputTag = cms.InputTag("l1extraParticles","Central"),      # for Run-1 algos
     	L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
        # cuts on the tracks used to determined the zvertex of the jet (examples) :
        #ZMAX = cms.double( 25. ),      # in cm
        #CHI2MAX = cms.double( 100. ),
        #PTMINTRA = cms.double( 2. ),   # in GeV
)
process.pJets = cms.Path( process.L1TkJets )

#
# ---------------------------------------------------------------------------

# HT and MHT from L1TkHTMissProducer
# The collection of (Run 1) L1Jets have been created above, by unpacking the
# gctDigis and running process.L1Reco.

#process.L1TkHTMiss = cms.EDProducer("L1TkHTMissProducer",
#	L1TkJetInputTag = cms.InputTag("L1TkJets","Central"),
#	DeltaZ = cms.double( 999 ),   #  in mm. Here dummy cut, since I dont have the zvtx of the jets 
#	PrimaryVtxConstrain = cms.bool( True ),
#	L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex")
#)
#process.pHTM = cms.Path( process.L1TkHTMiss )


# ---------------------------------------------------------------------------
#
# --- Run the analyzer

process.ana = cms.EDAnalyzer( 'L1TrackTriggerObjectsAnalyzer' ,
    L1VtxInputTag = cms.InputTag("L1TrackPrimaryVertex"),
    L1TkEtMissInputTag = cms.InputTag("L1TkEtMiss","MET"),
    L1TkElectronsInputTag = cms.InputTag("L1TkElectronsTrack","EG"),
    L1TkPhotonsInputTag = cms.InputTag("L1TkPhotons","EGIsoTrk"),
    L1TkJetsInputTag = cms.InputTag("L1TkJets","Central"),
    L1TkHTMInputTag = cms.InputTag("L1TkHTMiss")
)


# root file with histograms produced by the analyzer
# (mostly empty currently)

filename = "ana.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True)
)

process.pAna = cms.Path( process.ana )

#
# ---------------------------------------------------------------------------






process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_all.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_ALL' )
process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_*gen*_*_*')
process.Out.outputCommands.append('keep *_*Gen*_*_*')
process.Out.outputCommands.append('keep *_rawDataCollector_*_*')
process.Out.outputCommands.append('keep *_L1TkStubsFromPixelDigis_StubsPass_*')

process.FEVToutput_step = cms.EndPath(process.Out)




