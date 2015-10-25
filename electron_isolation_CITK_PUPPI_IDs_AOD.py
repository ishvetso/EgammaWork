import FWCore.ParameterSet.Config as cms

process = cms.Process( "ElectronIsolation" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    src = cms.InputTag("pfNoPileUpIso"),
    pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212, 22, 111, 130, 310, 2112, 11, -11, 13, -13),
    makeClones = cms.bool(True)
)

process.load("CommonTools.PileupAlgos.Puppi_cff")
process.load("EgammaWork.ElectronNtupler.pfNoLeptons_cfi")
process.pfNoLeptons.src = cms.InputTag('pfNoPileUpCandidates')

process.puppi.candName = cms.InputTag('pfNoPileUpCandidates')
process.puppi.vertexName = cms.InputTag('offlinePrimaryVertices')
process.puppi.puppiForLeptons = False

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

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
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ElectronIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
					    srcToIsolate = cms.InputTag("gedGsfElectrons"),
					    srcForIsolationCone = cms.InputTag('pfNoPileUpCandidates'),
					    #puppiValueMap = cms.InputTag('puppiNoLeptons'),
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


									
process.ElectronIsolationOnPUPPI = cms.EDProducer("CITKPFIsolationSumProducer",
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

process.ElectronIsolationOnPUPPINoLeptons = cms.EDProducer("CITKPFIsolationSumProducer",
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
									
									
									
process.ntupler = cms.EDAnalyzer('ElectronNtupler',
				 pruned = cms.InputTag("genParticles"),
				 pileup = cms.InputTag("addPileupInfo"),
				 vertices = cms.InputTag("offlinePrimaryVertices"),
				 electrons = cms.InputTag("gedGsfElectrons"),
				 cand_src = cms.InputTag("pfNoPileUpCandidates"),
				 rho = cms.InputTag("fixedGridRhoFastjetAll"),
				 #CITK
				 ValueMaps_ChargedHadrons_src = cms.InputTag("ElectronIsolation", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_NeutralHadrons_src = cms.InputTag("ElectronIsolation", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_Photons_src = cms.InputTag("ElectronIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
				 #PUPPI
				 ValueMaps_PUPPI_ChargedHadrons_src = cms.InputTag("ElectronIsolationOnPUPPI", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NeutralHadrons_src = cms.InputTag("ElectronIsolationOnPUPPI", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_Photons_src = cms.InputTag("ElectronIsolationOnPUPPI", "gamma-DR030-BarVeto000-EndVeto008"),
				  #PUPPINoLeptons
				 ValueMaps_PUPPI_NoLeptons_ChargedHadrons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "h+-DR030-BarVeto000-EndVeto001"),
				 ValueMaps_PUPPI_NoLeptons_NeutralHadrons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "h0-DR030-BarVeto000-EndVeto000"),
				 ValueMaps_PUPPI_NoLeptons_Photons_src = cms.InputTag("ElectronIsolationOnPUPPINoLeptons", "gamma-DR030-BarVeto000-EndVeto008"),
				 # for electron ids
				 mva_idw80_src = cms.InputTag("egmGsfElectronIDs", "mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
				 mva_idw90_src = cms.InputTag("egmGsfElectronIDs", "mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),
  				 #generator info 
				 genInfo = cms.InputTag("generator")
			)

process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)


process.electrons = cms.Path(process.particleFlowTmpPtrs + process.pfParticleSelectionSequence + process.pfNoPileUpCandidates + process.egmGsfElectronIDSequence + process.pfNoLeptons + process.puppi + process.puppiNoLeptons + process.ElectronIsolation + process.ElectronIsolationOnPUPPI + process.ElectronIsolationOnPUPPINoLeptons  + process.ntupler)


process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/i/ishvetso/EgammaWork/test_samples/ttbar_AOD.root'),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


'''process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string('patTuple.root'),
  outputCommands = cms.untracked.vstring('keep *')
)

process.outpath = cms.EndPath(process.out)
'''
process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )
