import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST_forElec/SingleElectron/BE5D/zmatchingOff/L1TkElectrons_SingleElectron_BE5D_v0.root'
    #'file:L1TrackElectron.root'
    #'file:example_withEtMiss.root'
    #'file:example_w_Tracks.root'
    )
)




process.ana = cms.EDAnalyzer( 'ElectronInefficiencies' ,
    #L1TkElectronsInputTag = cms.InputTag("L1TkElectrons","EG"),
    L1TkElectronsInputTag = cms.InputTag("L1TkElectronsLoose","EG"),
            L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","EGamma"),      # input EGamma collection
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	BREM_CUT = cms.double( 0.10 )
)


# root file with histograms produced by the analyzer
filename = "ana.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path( process.ana )




