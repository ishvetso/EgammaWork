import FWCore.ParameterSet.Config as cms

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/RunII_preparation/samples/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8_PHYS14.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Geometry_cff")

process.GlobalTag.globaltag = 'PHYS14_25_V2::All'

process.egmPhotonIsolationMiniAOD = cms.EDProducer(
    "CITKPFIsolationSumProducer",
    srcToIsolate = cms.InputTag("slimmedPhotons"),
    srcForIsolationCone = cms.InputTag('packedPFCandidates'),
    isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(1,2,3),
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

process.treeMaker = cms.EDAnalyzer("PhotonValidator",
				   photons = cms.InputTag("slimmedPhotons"),
				   phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
				   phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
				   phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
				   #CITK
				   phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h+-DR030-"),
				   phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h0-DR030-"),
				   phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:gamma-DR030-")
      
    )			   

process.analysis = cms.Path(process.egmPhotonIsolationMiniAOD + process.photonIDValueMapProducer + process.treeMaker)


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
