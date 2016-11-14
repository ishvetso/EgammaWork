import FWCore.ParameterSet.Config as cms

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring16DR80/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/02934727-2210-E611-AAC4-2C600CAFEF7C.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("CommonTools.PileupAlgos.Puppi_cff")

from RecoEgamma.EgammaIsolationAlgos.egmPhotonIsolationPUPPI_cff import egmPhotonIsolationAODPUPPI

process.egmPhotonIsolationAODPUPPI = egmPhotonIsolationAODPUPPI

process.ntupler = cms.EDAnalyzer('SimplePhotonNtupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #                                    
                                 rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 photons = cms.InputTag("gedPhotons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 #Value maps from CITK
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationAODPUPPI:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationAODPUPPI:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationAODPUPPI:gamma-DR030-"),
                                 genInfo = cms.InputTag("generator")
                                )			   

process.analysis = cms.Path(process.puppi + process.egmPhotonIsolationAODPUPPI +  process.ntupler)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_AOD.root")
                                  )
