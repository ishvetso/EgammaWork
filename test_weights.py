import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.load("CommonTools.PileupAlgos.Puppi_cff")


process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

									
process.SelectedElectrons = cms.EDFilter('ElectronSelector',
				 electron_src = cms.InputTag("slimmedElectrons"),
				 cand_src = cms.InputTag("packedPFCandidates"),
				)									
									
									
process.ElectronTree = cms.EDAnalyzer('PFCandidateAnalyzer',
				 leptons = cms.InputTag("SelectedElectrons"),
				 cand_src = cms.InputTag("packedPFCandidates"),
				 puppiValueMap = cms.InputTag("puppi"),
				)

process.pfNoLeptons = cms.EDProducer("PdgIDSelector",
    cand_src = cms.InputTag("packedPFCandidates"),
    pdgIDs = cms.vuint32(211,321,999211,2212,111,130,310,2112,22)
)



process.electrons = cms.Path(process.SelectedElectrons + process.puppi  + process.ElectronTree)


process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/EgammaWork/test_samples/DY.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_ele_puppi.root")
                                  )
