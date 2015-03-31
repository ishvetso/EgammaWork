import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.ElectronIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("slimmedElectrons"),
					    srcForIsolationCone = cms.InputTag('packedPFCandidates'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.015),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.0),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.08),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(2,3) )
								      )
					)
process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 packed = cms.InputTag("packedGenParticles"),
				 pruned = cms.InputTag("prunedGenParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
				 electrons = cms.InputTag("slimmedElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("ElectronIsolation", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("ElectronIsolation", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("ElectronIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)

process.electrons = cms.Path(process.ElectronIsolation  + process.ntupler)

#process.maxEvents.input = 1000
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:////afs/cern.ch/work/i/ishvetso/RunII_preparation/Synchronization_March2015/miniAOD/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8_synch_exercise.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

'''process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)
'''
process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
