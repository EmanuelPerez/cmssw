import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.L1HwVal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_forTTI_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
  '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/Electrons/NoPU_CSM/SingleElectronFlatPt0p2To50_RAW-DIGI1.root'
                    #'file:./SingleElectronFlatPt0p2To50_RAW-DIGI1.root'
                  #  'file:./test.root'
                              )
        )

process.L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
   EtminForStore = cms.double( -999.),
   debug = cms.untracked.bool(False),
   useECalEndcap = cms.bool(True)
)

#################################################################################################
# Analyzer
#################################################################################################
process.NtupleMaker = cms.EDAnalyzer('L1PixelTrigger',
            tau = cms.InputTag("SLHCL1ExtraParticles","Taus"),
            egamma = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
            egcrystal = cms.InputTag("L1EGammaCrystalsProducer","EGammaCrystal")
            )

process.p = cms.Path(
        process.RawToDigi+
            process.valHcalTriggerPrimitiveDigis+
                process.SLHCCaloTrigger+
                    process.siPixelRecHits+
                     process.calolocalreco+
                      process.L1EGammaCrystalsProducer+
                        process.NtupleMaker
                        )



#################################################################################################
# Output file
#################################################################################################
#process.TFileService = cms.Service("TFileService", fileName = cms.string('Single_ElectronFlatPt0_50_with_eGammaposition_620_SLHC10_1.root') )
process.TFileService = cms.Service("TFileService", fileName = cms.string('620_patch2.root') )

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')


from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
process.siPixelRecHits = siPixelRecHits

process.raw2digi_step = cms.Path(process.RawToDigi)
process.mix.digitizers = cms.PSet(process.theDigitizersValid)

# bug fix for missing HCAL TPs in MC RAW
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("valHcalTriggerPrimitiveDigis")

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI

#call to customisation function cust_2023TTI imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023TTI(process)

