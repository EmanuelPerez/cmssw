import FWCore.ParameterSet.Config as cms

process = cms.Process("ETMISS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:example_w_Tracks_and_vertex.root'   # created by runL1Tracks_and_L1PrimaryVertexProducer_cfg.py
    )
)


# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

        # run L1Reco to produce the L1ETM object corresponding
        # to the current trigger
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

	# now produce the TrkMET :
process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkEtMissProducer_cfi")
process.p = cms.Path( process.L1TkEtMiss )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_withEtMiss.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append( 'keep *_*_*_ETMISS' )
process.Out.outputCommands.append('keep *_generator_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




