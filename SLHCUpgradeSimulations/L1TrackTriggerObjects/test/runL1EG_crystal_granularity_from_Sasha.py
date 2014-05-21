#
#
#
#     	- Runs the stage-2 L1Calo algorithms (both old 2x2 clustering and	
#     	  the new clustering from LLR).
#     	- Runs the crystal-based L1EG algorithm of Sasha Savin.
#	- Runs the producer that matches the L1EH objects with the crystal-based clusters.
#

# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("L1EG")

from SLHCUpgradeSimulations.L1TrackTriggerObjects.minBiasFiles_p1_cfi import *


process.source = cms.Source("PoolSource",
   fileNames = minBiasFiles_p1
   #fileNames = cms.untracked.vstring("root://eoscms///store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    #reportEvery = cms.untracked.int32(500),
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(10000000)
)      
       
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# Load geometry
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
                            
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


process.load("Configuration.StandardSequences.Services_cff")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff") ###check this for MC!
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")


# --------------------------------------------------------------------------------------------
#
# ----    Runs the stage-2 L1EG algorithms
#

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
	# to the current trigger
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.L1Reco = cms.Path( process.l1extraParticles )

# --------------------------------------------------------------------------------------------
#
# ----    Produce the L1EGCrystal clusters (code of Sasha Savin)

	# first you need the ECAL RecHIts :
process.reconstruction_step = cms.Path( process.calolocalreco )

process.L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
   DEBUG = cms.bool(False)
)
process.pSasha = cms.Path( process.L1EGammaCrystalsProducer )

# --------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------
#
# ----  Match the L1EG stage-2 objects created by the SLHCCaloTrigger sequence above
#	with the crystal-level clusters.
#	This produces a new collection of L1EG objects, starting from the original
#	L1EG collection. The eta and phi of the L1EG objects is corrected using the
#	information of the xtal level clusters.

process.l1ExtraCrystalProducer = cms.EDProducer("L1ExtraCrystalPosition",
   eGammaSrc = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
   eClusterSrc = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster")
)


process.egcrystal_producer = cms.Path(process.l1ExtraCrystalProducer)



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1EG.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

#process.Out.outputCommands.append( 'keep *_*_*_L1EG' )
process.Out.outputCommands.append('keep *_generator_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_*_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_*_*')
process.Out.outputCommands.append('keep *_L1EGammaCrystalsProducer_*_*')

process.Out.outputCommands.append('keep *_l1extraParticles_*_*')
process.Out.outputCommands.append('keep *_l1ExtraCrystalProducer_*_*')


#process.Out.outputCommands.append('keep *')

process.FEVToutput_step = cms.EndPath(process.Out)



