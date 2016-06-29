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
    fileNames = cms.untracked.vstring('file:///afs/cern.ch/user/i/ishvetso/eos/cms/store/mc/RunIISpring16DR80/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/6AC91184-18FA-E511-9906-A4BADB1E6F7A.root')
)

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = process.pfAllChargedHadrons.clone()
process.pfNoPileUpCandidates.pdgId.extend(process.pfAllNeutralHadronsAndPhotons.pdgId)

process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)
from RecoEgamma.EgammaIsolationAlgos.egmPhotonIsolationAOD_cff import egmPhotonIsolationAOD

process.egmPhotonIsolationAOD = egmPhotonIsolationAOD.clone()

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
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:gamma-DR030-"),
                                 # 
                                 # Locations of files with the effective area constants.
                                 # The constants in these files below are derived for PHYS14 MC.
                                 #
                                 effAreaChHadFile = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfChargedHadrons_25ns_NULLcorrection.txt"),
                                 effAreaNeuHadFile= cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfNeutralHadrons_25ns_90percentBased.txt"),
                                 effAreaPhoFile   = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfPhotons_25ns_90percentBased.txt"),
                                 genInfo = cms.InputTag("generator"),
                                 vertices = cms.InputTag("offlinePrimaryVertices")
                                )			   

process.analysis = cms.Path(process.particleFlowTmpPtrs +  process.pfParticleSelectionSequence + process.pfNoPileUpCandidates + process.egmPhotonIsolationAOD +  process.ntupler)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_AOD.root")
                                  )
