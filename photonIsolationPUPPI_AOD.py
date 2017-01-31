import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = process.pfAllChargedHadrons.clone()
process.pfNoPileUpCandidates.pdgId.extend(process.pfAllNeutralHadronsAndPhotons.pdgId)



process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///afs/cern.ch/user/i/ishvetso/eos/cms/store/mc/RunIISpring16DR80/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/70000/38413A28-080C-E611-A5B4-02163E012D0D.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#Loading PUPPI sequences
process.load("CommonTools.PileupAlgos.Puppi_cff")

process.puppi.candName = cms.InputTag('particleFlow')
process.puppi.vertexName = cms.InputTag('offlinePrimaryVertices')
process.puppi.puppiForLeptons = False


process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)


process.egmPhotonIsolation = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("gedPhotons"),
			  srcForIsolationCone = cms.InputTag('pfNoPileUpCandidates'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    ),
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h0'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    ),
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('gamma'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    )
    )
  )	

process.egmPhotonIsolationPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
			  srcToIsolate = cms.InputTag("gedPhotons"),
			  srcForIsolationCone = cms.InputTag('particleFlow'),
              puppiValueMap = cms.InputTag('puppi'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    ),
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h0'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    ),
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('gamma'),
				      miniAODVertexCodes = cms.vuint32(2,3),
				      vertexIndex = cms.int32(0),
				    )
    )
  )

	   

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

                                 candidates = cms.InputTag("particleFlow"),
                                 vertices = cms.InputTag("offlinePrimaryVertices"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),                                 
                                 #CITK
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolation:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolation:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolation:gamma-DR030-"),
                                 #PUPPI
                                 phoChargedIsolation_PUPPI = cms.InputTag("egmPhotonIsolationPUPPI:h+-DR030-"),
                                 phoNeutralHadronIsolation_PUPPI = cms.InputTag("egmPhotonIsolationPUPPI:h0-DR030-"),
                                 phoPhotonIsolation_PUPPI = cms.InputTag("egmPhotonIsolationPUPPI:gamma-DR030-"),
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
                                  genInfo = cms.InputTag("generator"),
                                )			  

process.analysis = cms.Path(process.particleFlowTmpPtrs +  process.pfParticleSelectionSequence + process.pfNoPileUpCandidates +  process.puppi + process.egmPhotonIsolation + process.egmPhotonIsolationPUPPI + process.ntupler)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
