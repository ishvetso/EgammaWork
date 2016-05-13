import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv1/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/06F2069A-40BC-E511-BD4F-34E6D7BEAF28.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#Loading PUPPI sequences
process.load("CommonTools.PileupAlgos.Puppi_cff")

process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.puppiForLeptons = False



process.egmPhotonIsolationAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("slimmedPhotons"),
			  srcForIsolationCone = cms.InputTag('packedPFCandidates'),
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

process.egmPhotonIsolationAODPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
			  srcToIsolate = cms.InputTag("slimmedPhotons"),
			  srcForIsolationCone = cms.InputTag('packedPFCandidates'),
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

                                 candidates = cms.InputTag("packedPFCandidates"),
                                 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),                                 
                                 #CITK
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationAOD:gamma-DR030-"),
                                 #PUPPI
                                 phoChargedIsolation_PUPPI = cms.InputTag("egmPhotonIsolationAODPUPPI:h+-DR030-"),
                                 phoNeutralHadronIsolation_PUPPI = cms.InputTag("egmPhotonIsolationAODPUPPI:h0-DR030-"),
                                 phoPhotonIsolation_PUPPI = cms.InputTag("egmPhotonIsolationAODPUPPI:gamma-DR030-"),
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

process.analysis = cms.Path(  process.puppi + process.egmPhotonIsolationAOD + process.egmPhotonIsolationAODPUPPI   + process.ntupler)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
