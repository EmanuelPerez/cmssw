import FWCore.ParameterSet.Config as cms

process = cms.Process("Tau")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500),
    limit = cms.untracked.int32(10000000)
)

#
# This runs over a file that already contains the L1Tracks.
#
# Creates L1TkTaus starting from the stage-2 L1CaloTaus.
#

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
     '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/HTauTau/BE5D/zmatchingOff/m1_HTauTau_BE5D.root'
    #'file:example_w_Tracks_and_vertex.root'
    )
)


process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# -------------------------------------------------------------------------

# -- Run the sequence SLHCCaloTrigger which creates "stage-2" L1Taus.
# -- Two collections are created, ("SLHCL1ExtraParticles","Taus")
# -- and ("SLHCL1ExtraParticles","IsoTaus").
# -- Some time ago we had checked that ("SLHCL1ExtraParticles","Taus")
# -- was reasonable; I do not think we ever looked into the "IsoTaus".

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

# -------------------------------------------------------------------------



# -------------------------------------------------------------------------

# --- Now run the L1TkTauParticle producer(s)

process.L1TkTauFromCalo = cms.EDProducer("L1TkTauFromCaloProducer",
	L1TausInputTag = cms.InputTag("SLHCL1ExtraParticles","Taus"),
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	ZMAX = cms.double( 25. ),	# in cm
	CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 2. ),	# in GeV
	DRmax = cms.double( 0.5 ),
     	nStubsmin = cms.int32( 5 )        # minimum number of stubs
)
process.pTaus = cms.Path( process.L1TkTauFromCalo )

# -------------------------------------------------------------------------



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TkTaus.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)


process.Out.outputCommands.append('keep *_L1TkTau*_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_*_*')
process.Out.outputCommands.append('keep *_gen*_*_*')


process.FEVToutput_step = cms.EndPath(process.Out)




