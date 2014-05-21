import FWCore.ParameterSet.Config as cms

process = cms.Process("Muo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500),
    limit = cms.untracked.int32(10000000)
)

#
# This runs over a file that already contains the L1Tracks.
#
# Creates L1TkMuons starting from the Run-1 L1Muons.
#

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
     '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/Muon/BE5D/zmatchingOff/Muon_BE5D.root'
    #'file:example_w_Tracks_and_vertex.root'
    )
)


process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# --- creates l1extra objects for L1Muons 

        # raw2digi to get the GT digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )
	# run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )


# --- Now run the L1TkMuonParticleProducer 

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkMuonProducer_cfi")
process.pMuons = cms.Path( process.L1TkMuons )
process.pMuons = cms.Path( process.L1TkMuons )



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TkMuons.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append('keep *')

process.Out.outputCommands.append('keep *_L1TkMuons*_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')
process.Out.outputCommands.append('keep *_gen*_*_*')


process.FEVToutput_step = cms.EndPath(process.Out)




