import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False)

from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import egmGedGsfElectronPFNoPileUpIsolation
process.ElectronIsolation = egmGedGsfElectronPFNoPileUpIsolation.clone()
process.ElectronIsolation.srcToIsolate = cms.InputTag("slimmedElectrons")
process.ElectronIsolation.srcForIsolationCone = cms.InputTag("packedPFCandidates")



process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("prunedGenParticles"),
				 pileup = cms.InputTag("slimmedAddPileupInfo"),
				 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
				 electrons = cms.InputTag("slimmedElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("ElectronIsolation", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("ElectronIsolation", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("ElectronIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)

process.electrons = cms.Path(process.ElectronIsolation  + process.ntupler)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/02EA1AC3-40B8-E511-8AE8-3417EBE7009F.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_miniAOD.root")
                                  )
