import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("CommonTools.PileupAlgos.Puppi_cff")
process.load("EgammaWork.ElectronNtupler.pfNoLeptons_cfi")

process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.puppiForLeptons = False

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.puppiNoLeptons = process.puppi.clone()
process.puppiNoLeptons.puppiForLeptons = False
process.puppiNoLeptons.candName = cms.InputTag('pfNoLeptons')

process.ElectronIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("SelectedElectrons"),
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

process.ElectronIsolationMapBasedVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("SelectedElectrons"),
					    srcForIsolationCone = cms.InputTag('packedPFCandidates'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithMapBasedVeto'), 
									coneSize = cms.double(0.3),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithMapBasedVeto'), 
									coneSize = cms.double(0.3),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithMapBasedVeto'), 
									coneSize = cms.double(0.3),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(2,3) )
								      )
					)
									
process.ElectronIsolationOnPUPPI = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("SelectedElectrons"),
					    srcForIsolationCone = cms.InputTag('puppi'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.015),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.0),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.08),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(1,2,3) )
								      )
					)

process.ElectronIsolationOnPUPPINoLeptons = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("SelectedElectrons"),
					    srcForIsolationCone = cms.InputTag('puppiNoLeptons'),
					    isolationConeDefinitions = cms.VPSet(
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.015),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h+'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.0),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('h0'),
									miniAODVertexCodes = cms.vuint32(1,2,3) ),
									cms.PSet( isolationAlgo = cms.string('ElectronPFIsolationWithConeVeto'), 
									coneSize = cms.double(0.3),
									VetoConeSizeEndcaps = cms.double(0.08),
									VetoConeSizeBarrel = cms.double(0.0),
									isolateAgainst = cms.string('gamma'),
									miniAODVertexCodes = cms.vuint32(1,2,3) )
								      )
					)									
									
									
									
process.ntupler = cms.EDAnalyzer('ElectronNtuplerTest',
				 pruned = cms.InputTag("prunedGenParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
				 electrons = cms.InputTag("SelectedElectrons"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("ElectronIsolation", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("ElectronIsolation", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("ElectronIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
				 #CITK map based
				 ValueMaps_ChargedHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "h+-DR030-"),
				 ValueMaps_NeutralHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "h0-DR030-"),
				 ValueMaps_Photons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "gamma-DR030-"),
				 #PUPPI
				 ValueMaps_PUPPI_ChargedHadrons_src = cms.InputTag("ElectronIsolationOnPUPPI", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NeutralHadrons_src = cms.InputTag("ElectronIsolationOnPUPPI", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_Photons_src = cms.InputTag("ElectronIsolationOnPUPPI", "gamma-DR030-BarVeto000-EndVeto008"),
				  #PUPPINoLeptons
				 ValueMaps_PUPPI_NoLeptons_ChargedHadrons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NoLeptons_NeutralHadrons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_NoLeptons_Photons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "gamma-DR030-BarVeto000-EndVeto008"),
								
				)

process.SelectedElectrons = cms.EDFilter('ElectronSelector',
				 electron_src = cms.InputTag("slimmedElectrons"),
				 cand_src = cms.InputTag("packedPFCandidates"),
				)


process.electrons = cms.Path( process.SelectedElectrons + process.pfNoLeptons +  process.puppi + process.puppiNoLeptons + process.ElectronIsolation + process.ElectronIsolationMapBasedVeto + process.ElectronIsolationOnPUPPI + process.ElectronIsolationOnPUPPINoLeptons + process.ntupler)


process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/RunII_preparation/samples/WW_74X.root')
    
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
'''

process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)'''

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
