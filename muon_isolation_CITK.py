import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.load("CommonTools.PileupAlgos.Puppi_cff")
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.MuonIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("slimmedMuons"),
					    srcForIsolationCone = cms.InputTag('packedPFCandidates'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(2,3) )
								      )
					)

process.MuonIsolationPUPPI = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("slimmedMuons"),
					    srcForIsolationCone = cms.InputTag('puppi'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('MuonIsolation'), 
									coneSize = cms.double(0.4),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(2,3) )
								      )
					)
process.ntupler = cms.EDAnalyzer('MuonNtupler',
				 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
				 muons = cms.InputTag("slimmedMuons"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("MuonIsolation", "h+-DR040-"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("MuonIsolation", "h0-DR040-"),
				 ValueMaps_Photons_src = cms.InputTag("MuonIsolation", "gamma-DR040-"),
				 ValueMaps_ChargedHadrons_PUPPI_src = cms.InputTag("MuonIsolationPUPPI", "h+-DR040-"),
				 ValueMaps_NeutralHadrons_PUPPI_src = cms.InputTag("MuonIsolationPUPPI", "h0-DR040-"),
				 ValueMaps_Photons_PUPPI_src = cms.InputTag("MuonIsolationPUPPI", "gamma-DR040-"),
				 genInfo = cms.InputTag("generator")
								
				)

process.electrons = cms.Path(process.puppi + process.MuonIsolation + process.MuonIsolationPUPPI + process.ntupler  )

#process.maxEvents.input = 1000
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/EgammaWork/test_samples/DY_2.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

'''process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)'''

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
