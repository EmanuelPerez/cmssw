import FWCore.ParameterSet.Config as cms
process = cms.Process("Ele")

process.load("FWCore.MessageService.MessageLogger_cfi")

################################################################################
# Example configuratiom file : 
#
# Here we run the L1EG algorithms (old stage-2 and new clustering),
# and we create L1TkElectron objects starting from the "old stage-2" L1EGs.
#
# The L1Tracking is also run here.
#
################################################################################

# list of files
file_names = cms.untracked.vstring(
 #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/MinBias/BE5D/m1_MinBias_BE5D.root')
 '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00D6C34E-0339-E311-836A-002618943880.root',
 '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FEDB1C0F-FF38-E311-A659-0025905938D4.root')
#)
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
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
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

process.L1TkElectrons = cms.EDProducer("L1TkElectronTrackProducer",
        #label = cms.string("ElecTrk"),
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	label = cms.string("EG"),	# labels the collection of L1TkEmParticleProducer that is produced.
                                        # e.g. EG or IsoEG if all objects are kept, or
                                        # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                        # objects that pass a cut RelIso < RelIsoCut are written
                                        # into the new collection.
        L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","EGamma"),      # input EGamma collection
					# When the standard sequences are used :
                                                #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
                                                #   - for the "old stage-2" algo (2x2 clustering), use 
                                                #     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
        ETmin = cms.double( -1. ),	        # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
	# Track selection :
        CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 12. ),   # in GeV
        TrackEGammaDeltaPhi = cms.double(0.08),  # Delta Phi cutoff to match Track with L1EG objects
        TrackEGammaDeltaR = cms.double(0.08),   # Delta R cutoff to match Track with L1EG objects
        TrackEGammaDeltaEta = cms.double(0.08), # Delta Eta cutoff to match Track with L1EG objects
                                                # are considered. 
	RelativeIsolation = cms.bool( True ),	# default = True. The isolation variable is relative if True,
						# else absolute.
        IsoCut = cms.double( -0.15 ), 		# Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
						# When RelativeIsolation = False, IsoCut is in GeV.
						# NOTE: TRKBASED-ISOLATION IS NOT READY TO USE YET.
        # Determination of the isolation w.r.t. L1Tracks :
	DRmin = cms.double( 0.06),
	DRmax = cms.double( 0.5 ),
	DeltaZ = cms.double( 1.0 )    # in cm. Used for tracks to be used isolation calculation
)
process.pElectrons = cms.Path( process.L1TkElectrons )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TrackElectron.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_SLHCL1ExtraParticles_EGamma_*' )
process.Out.outputCommands.append( 'keep *_L1TkElectrons_*_*' )
process.Out.outputCommands.append( 'keep *_genParticles_*_*')
process.Out.outputCommands.append('keep *_generator_*_*')

#process.Out.outputCommands.append( 'keep *_L1TkElectrons_ElecTrk_*' )
#process.Out.outputCommands.append( 'keep SimTracks_g4SimHits_*_*') 
#process.Out.outputCommands.append('keep *')

process.FEVToutput_step = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.pSLHCCalo,process.TT_step,process.pElectrons,process.FEVToutput_step)




