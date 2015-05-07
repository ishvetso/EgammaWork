import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("CommonTools.PileupAlgos.Puppi_cff")
process.load("EgammaWork.ElectronNtupler.pfNoLeptons_cfi")

process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

									
									
									
									
process.ntupler = cms.EDAnalyzer('Test_pt',
				 puppi = cms.InputTag("puppi"),
				 pf = cms.InputTag("packedPFCandidates"),
				)

process.electrons = cms.Path(process.puppi + process.ntupler)


process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_updated28April2014/test_samples/ttbar.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_no_hack.root")
                                  )
