import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 



process.load("EgammaWork.ElectronNtupler.pfNoLeptons_cfi")
process.pfNoLeptons.src = cms.InputTag('particleFlow')
process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = process.pfAllChargedHadrons.clone()
process.pfNoPileUpCandidates.pdgId.extend(process.pfAllNeutralHadronsAndPhotons.pdgId)



process.load("CommonTools.PileupAlgos.Puppi_cff")
process.puppi.candName = cms.InputTag('particleFlow')
process.puppi.vertexName = cms.InputTag('offlinePrimaryVertices')
process.puppi.puppiForLeptons = False



process.puppiNoLeptons = process.puppi.clone()
process.puppiNoLeptons.puppiForLeptons = False
process.puppiNoLeptons.candName = cms.InputTag('pfNoLeptons')


useAOD = True
#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#standard isolation with cone veto
process.ElectronIsolationConeVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('pfNoPileUpCandidates'),
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

#standard isolation with map based veto
process.ElectronIsolationMapBasedVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('pfNoPileUpCandidates'),
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



#puppi isolation with cone veto									
process.ElectronIsolationOnPUPPIConeVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('puppi'),
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

#puppi no leptons with cone veto
process.ElectronIsolationOnPUPPINoLeptonsConeVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('puppiNoLeptons'),
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
									

#puppi isolation with map based veto									
process.ElectronIsolationOnPUPPIMapBasedVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('puppi'),
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

#puppi no leptons with map based veto
process.ElectronIsolationOnPUPPINoLeptonsMapBasedVeto = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('puppiNoLeptons'),
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
									
									
process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("genParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlinePrimaryVertices"),
				 electrons = cms.InputTag("gedGsfElectrons"),
				 cand_src = cms.InputTag("particleFlow"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationConeVeto", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationConeVeto", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_ConeVeto_src = cms.InputTag("ElectronIsolationConeVeto", "gamma-DR030-BarVeto000-EndVeto008"),

				 ValueMaps_ChargedHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "h+-DR030-"),
				 ValueMaps_NeutralHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "h0-DR030-"),
				 ValueMaps_Photons_MapBasedVeto_src = cms.InputTag("ElectronIsolationMapBasedVeto", "gamma-DR030-"),
				 #PUPPI
				 ValueMaps_PUPPI_ChargedHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPIConeVeto", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NeutralHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPIConeVeto", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_Photons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPIConeVeto", "gamma-DR030-BarVeto000-EndVeto008"),

				  #PUPPI
				 ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPIMapBasedVeto", "h+-DR030-"),
				 ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPIMapBasedVeto", "h0-DR030-"),
				 ValueMaps_PUPPI_Photons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPIMapBasedVeto", "gamma-DR030-"),
				  #PUPPINoLeptons
				 ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsConeVeto", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsConeVeto", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsConeVeto", "gamma-DR030-BarVeto000-EndVeto008"),

				 ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsMapBasedVeto", "h+-DR030-"),
				 ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsMapBasedVeto", "h0-DR030-"),
				 ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptonsMapBasedVeto", "gamma-DR030-"),
				 # for electron ids
				 mva_idw80_src = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
				 mva_idw90_src = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  				 #generator info 
				 genInfo = cms.InputTag("generator"),
				 effAreaFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
				)

process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)


process.electrons = cms.Path(process.particleFlowTmpPtrs + process.pfParticleSelectionSequence + process.pfNoPileUpCandidates + process.egmGsfElectronIDSequence + process.pfNoLeptons + process.puppi + process.puppiNoLeptons + process.ElectronIsolationConeVeto + process.ElectronIsolationMapBasedVeto + process.ElectronIsolationOnPUPPIConeVeto + process.ElectronIsolationOnPUPPINoLeptonsConeVeto + process.ElectronIsolationOnPUPPIMapBasedVeto + process.ElectronIsolationOnPUPPINoLeptonsMapBasedVeto + process.ntupler)


process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:///afs/cern.ch/user/i/ishvetso/eos/cms/store/mc/RunIIFall15DR76/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/5AAC31E2-1FC0-E511-994B-02163E00C3AD.root'),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
