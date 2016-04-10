import FWCore.ParameterSet.Config as cms

process = cms.Process( "PhotonIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv1/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/167B87FB-84AC-E511-A144-0CC47A745298.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from RecoEgamma.EgammaIsolationAlgos.egmPhotonIsolationMiniAOD_cff import egmPhotonIsolationMiniAOD

process.egmPhotonIsolationMiniAOD = egmPhotonIsolationMiniAOD.clone()


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
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:gamma-DR030-"),
                                 # 
                                 # Locations of files with the effective area constants.
                                 # The constants in these files below are derived for PHYS14 MC.
                                 #
                                 effAreaChHadFile = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                 effAreaNeuHadFile= cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                 effAreaPhoFile   = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),
                                 genInfo = cms.InputTag("generator")
                                )			   

process.analysis = cms.Path(process.egmPhotonIsolationMiniAOD +  process.ntupler)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_miniAOD.root")
                                  )
