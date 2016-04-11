import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)


from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = process.pfAllChargedHadrons.clone()
process.pfNoPileUpCandidates.pdgId.extend(process.pfAllNeutralHadronsAndPhotons.pdgId)

from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import egmGedGsfElectronPFNoPileUpIsolation
process.ElectronIsolation = egmGedGsfElectronPFNoPileUpIsolation.clone()

process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("genParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlinePrimaryVertices"),
				 electrons = cms.InputTag("gedGsfElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("ElectronIsolation", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("ElectronIsolation", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("ElectronIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)



process.electrons = cms.Path(process.particleFlowTmpPtrs + process.pfParticleSelectionSequence + process.pfNoPileUpCandidates + process.ElectronIsolation + process.ntupler)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15DR76/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/08DC4220-16A7-E511-AF59-1CC1DE19286E.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
