import FWCore.ParameterSet.Config as cms

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/RunII_preparation/samples/WW_74X.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

#Loading PUPPI sequences
process.load("CommonTools.PileupAlgos.Puppi_cff")
process.load("EgammaWork.ElectronNtupler.pfNoLeptons_cfi")

process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.puppiForLeptons = True

process.egmPhotonIsolationMiniAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("slimmedPhotons"),
			  srcForIsolationCone = cms.InputTag('packedPFCandidates'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(3),
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

process.egmPhotonIsolationMiniAODPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
			  srcToIsolate = cms.InputTag("slimmedPhotons"),
			  srcForIsolationCone = cms.InputTag('packedPFCandidates'),
        puppiValueMap = cms.InputTag('puppi'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(3),
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
			   
process.photonIDValueMapProducer = cms.EDProducer('PhotonIDValueMapProducer',
                                          # The module automatically detects AOD vs miniAOD, so we configure both
                                          #
                                          # AOD case
                                          #
                                          ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                          eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                          esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                                          particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
                                          vertices = cms.InputTag("offlinePrimaryVertices"),
                                          pfCandidates = cms.InputTag("particleFlow"),
                                          src = cms.InputTag('gedPhotons'),
                                          #
                                          # miniAOD case
                                          #
                                          ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEBRecHits"),
                                          eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEERecHits"),
                                          esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedESRecHits"),
                                          verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                          pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
                                          # there is no need for the isolation map here, for miniAOD it is inside packedPFCandidates
                                          srcMiniAOD = cms.InputTag('slimmedPhotons'),
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
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 #
                                 # ValueMap names from the producer upstream
                                 #
                                 full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                 phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                 phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                 phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                 
                                 #CITK
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:gamma-DR030-"),
                                 #PUPPI
                                 phoChargedIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:h+-DR030-"),
                                 phoNeutralHadronIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:h0-DR030-"),
                                 phoPhotonIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:gamma-DR030-"),
                                 # 
                                 # Locations of files with the effective area constants.
                                 # The constants in these files below are derived for PHYS14 MC.
                                 #
                                 effAreaChHadFile = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                 effAreaNeuHadFile= cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                 effAreaPhoFile   = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt")
                                )			  

process.analysis = cms.Path(process.puppi + process.egmPhotonIsolationMiniAOD + process.egmPhotonIsolationMiniAODPUPPI + process.photonIDValueMapProducer )


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple_miniAOD.root'),
                               outputCommands = cms.untracked.vstring('keep *')
                               )                            

                           
process.outpath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
