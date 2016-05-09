import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.load("CommonTools.PileupAlgos.Puppi_cff")


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

process.egmElectronIsolaionPUPPIAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("gedGsfElectrons"),
			  srcForIsolationCone = cms.InputTag('puppi'),
			  puppiValueMap = cms.InputTag(''),
			  isolationConeDefinitions = IsoConeDefinitions
  )


process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("genParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlinePrimaryVertices"),
				 electrons = cms.InputTag("gedGsfElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("egmElectronIsolaionPUPPIAOD", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("egmElectronIsolaionPUPPIAOD", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("egmElectronIsolaionPUPPIAOD", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)



process.electrons = cms.Path( process.puppi + process.egmElectronIsolaionPUPPIAOD + process.ntupler)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/20000/021B31F4-33FC-E511-A597-02163E014626.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
