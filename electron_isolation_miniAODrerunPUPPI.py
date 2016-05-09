import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("CommonTools.PileupAlgos.Puppi_cff")
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False)

IsoConeDefinitions = cms.VPSet(
        cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
                  coneSize = cms.double(0.3),
                  VetoConeSizeBarrel = cms.double(0.0),
                  VetoConeSizeEndcaps = cms.double(0.015),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
                  coneSize = cms.double(0.3),
                  VetoConeSizeBarrel = cms.double(0.0),
                  VetoConeSizeEndcaps = cms.double(0.0),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3) ),
        cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
                  coneSize = cms.double(0.3),
                  VetoConeSizeBarrel = cms.double(0.0),
                  VetoConeSizeEndcaps = cms.double(0.08),
                  isolateAgainst = cms.string('gamma'),
                  miniAODVertexCodes = cms.vuint32(2,3) )
        )

process.egmElectronIsolaionPUPPIMiniAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("slimmedElectrons"),
			  srcForIsolationCone = cms.InputTag('puppi'),
			  isolationConeDefinitions = IsoConeDefinitions
  )



process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("prunedGenParticles"),
				 pileup = cms.InputTag("slimmedAddPileupInfo"),
				 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
				 electrons = cms.InputTag("slimmedElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("egmElectronIsolaionPUPPIMiniAOD", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("egmElectronIsolaionPUPPIMiniAOD", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("egmElectronIsolaionPUPPIMiniAOD", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)

process.electrons = cms.Path(process.puppi + process.egmElectronIsolaionPUPPIMiniAOD + process.ntupler)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/20000/0A21CE5A-8BFC-E511-B1F9-782BCB4086A8.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree_miniAOD.root")
                                  )
