import FWCore.ParameterSet.Config as cms

process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#
# This runs over a file that already contains the L1Tracks and the L1TrackPrimaryVertex.
# see runL1Tracks_and_L1PrimaryVertexProducer_cfg.py
#
# This configuration file:
#    - runs the SLHCCalo sequence, i.e. it runs the L1EG algorithms
#      (both the "old stage 2" and the "new clustering) and
#      the L1Jet algorithm 
#    - creates L1TkEmParticles and  L1TkElectrons.
#    - creates L1TkJets and L1TkHTMiss (i.e. HT and MHT) using for
#      calo jets the L1Jets created by the SLHCCalo sequence.
#    - creates the TrkMET 
#    - creates L1TkMuons using as input the L1Muons that were
#      produced dueing the central production 
#
#    - and runs a trivial analyzer that prints the L1TkObjects
#

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:example_w_Tracks_and_vertex.root'	 	# created by runL1Tracks_and_L1PrimaryVertexProducer_cfg.py
    )
)



# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# ---------------------------------------------------------------------------
#
# --- Run the L1TrackEtMiss producer :

# needs the L1TrackPrimaryVertex.
# TrkMET calculated from all tracks that have |z - z_vertex | < DeltaZ,
# pt > PTMINTRA and at least ( >= ) nStubsmin stubs.
# To check the latter condition, the file must contain the stubs.

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkEtMissProducer_cfi")
process.pEtMiss = cms.Path( process.L1TkEtMiss )

#
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#
# --- Run the L1EG algorithm of Jean-Baptiste (new clustering algorithm).
# --- Note thet the efficiency is poor at very high PU...
# --- This also runs the "old" stage-2 algorithm
#
# --- The sequence also creates stage-2 L1Jets.
#

process.load('Configuration/StandardSequences/L1HwVal_cff')
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
        

# ---------------------------------------------------------------------------
#
# --- produce the L1TrackEmParticles


# "photons" :
process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkEmParticleProducer_cfi")
process.pPhotons = cms.Path( process.L1TkPhotons )

# ---------------------------------------------------------------------------


# --- Now run the L1TkElectronProducers : one for the algoritm that
# --- matches with L1Tracks, the other for the algorithm that matches
# --- with stubs. Here we run only the matching with L1Tracks.

# "electrons" from L1Tracks :

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkElectronTrackProducer_cfi")
process.pElectronsTrack = cms.Path( process.L1TkElectrons )

# ---------------------------------------------------------------------------

# the L1TkJetProducer

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkJetProducer_cfi")
process.pJets = cms.Path( process.L1TkJets )

process.pJets = cms.Path( process.L1TkJets )


#
# ---------------------------------------------------------------------------

# HT and MHT from L1TkHTMissProducer, using the L1TkJets just created above :

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkHTMissProducer_cfi")
process.pHTMCalo = cms.Path( process.L1TkHTMissCalo )		# calo only, no vtx constraint
process.pHTMVtx  = cms.Path( process.L1TkHTMissVtx )	    # HT and MHT from jets that come from the same vertex


# ---------------------------------------------------------------------------

# L1TkMuons

        # first, run L1Reco to produce the L1Muon objects corresponding
        # to the current (Run-1) trigger. The simulation of the Run-1
	# L1 trigger has been run in the central production. L1Reco
	# uses the GT digis, that have been obtained from the unpakcing
	# step above.
	 
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )


# Now produce the L1TkMuons, by matching the L1Muons just created with L1Tracks.

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkMuonProducer_cfi")
process.pMuons = cms.Path( process.L1TkMuons )
 

# ---------------------------------------------------------------------------

#
# --- Run the analyzer

process.ana = cms.EDAnalyzer( 'L1TrackTriggerObjectsAnalyzer' ,
    L1VtxInputTag = cms.InputTag("L1TrackPrimaryVertex"),
    L1TkEtMissInputTag = cms.InputTag("L1TkEtMiss","MET"),
    L1TkElectronsInputTag = cms.InputTag("L1TkElectronsTrack","EG"),
    L1TkPhotonsInputTag = cms.InputTag("L1TkPhotons","IsoTrk"),
    L1TkJetsInputTag = cms.InputTag("L1TkJets","Central"),
    L1TkHTMInputTag = cms.InputTag("L1TkHTMissVtx")
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




