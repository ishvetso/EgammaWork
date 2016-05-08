import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process = cms.Process( "PhotonIsoTest" )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
from RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi import *

process.load("CommonTools.ParticleFlow.pfNoPileUpIso_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfNoPileUpCandidates = process.pfAllChargedHadrons.clone()
process.pfNoPileUpCandidates.pdgId.extend(process.pfAllNeutralHadronsAndPhotons.pdgId)



process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(False) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15DR76/GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/08D8F01F-538E-E511-8492-0CC47A78A3F4.root')
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'

useAOD = True

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


#Loading PUPPI sequences
process.load("CommonTools.PileupAlgos.Puppi_cff")

process.puppi.candName = cms.InputTag('particleFlow')
process.puppi.vertexName = cms.InputTag('offlinePrimaryVertices')
process.puppi.puppiForLeptons = False


process.particleFlowTmpPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
src = cms.InputTag('particleFlow')
)


process.egmPhotonIsolationMiniAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
			  srcToIsolate = cms.InputTag("gedPhotons"),
			  srcForIsolationCone = cms.InputTag('pfNoPileUpCandidates'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(2,3),
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

process.egmPhotonIsolationMiniAODPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
			  srcToIsolate = cms.InputTag("gedPhotons"),
			  srcForIsolationCone = cms.InputTag('particleFlow'),
              puppiValueMap = cms.InputTag('puppi'),
			  isolationConeDefinitions = cms.VPSet(
			   cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
				      coneSize = cms.double(0.3),
				      isolateAgainst = cms.string('h+'),
				      miniAODVertexCodes = cms.vuint32(2,3),
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

	   

process.ntupler = cms.EDAnalyzer('SimplePhotonNtupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #                                    
                                 rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 photons = cms.InputTag("gedPhotons"),
                                 genParticles = cms.InputTag("genParticles"),

                                 candidates = cms.InputTag("particleFlow"),
                                 vertices = cms.InputTag("offlinePrimaryVertices"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 photonsMiniAOD = cms.InputTag("slimmedPhotons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),                                 
                                 #CITK
                                 phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h+-DR030-"),
                                 phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:h0-DR030-"),
                                 phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD:gamma-DR030-"),
                                 #PUPPI
                                 phoChargedIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:h+-DR030-"),
                                 phoNeutralHadronIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:h0-DR030-"),
                                 phoPhotonIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI:gamma-DR030-"),
                                 # 
                                 # Locations of files with the effective area constants.
                                 # The constants in these files below are derived for PHYS14 MC.
                                 #
                                 effAreaChHadFile = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
                                 effAreaNeuHadFile= cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
                                 effAreaPhoFile   = cms.FileInPath
                                 ("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),
                                  genInfo = cms.InputTag("generator"),
                                  mva_idw90_src = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90")
                                )			  

process.analysis = cms.Path(process.egmPhotonIDSequence + process.particleFlowTmpPtrs +  process.pfParticleSelectionSequence + process.pfNoPileUpCandidates +  process.puppi + process.egmPhotonIsolationMiniAOD + process.egmPhotonIsolationMiniAODPUPPI   + process.ntupler)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1


'''process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple_miniAOD.root'),
                               outputCommands = cms.untracked.vstring('keep *')
                               )                            

                           
process.outpath = cms.EndPath(process.out)'''

process.TFileService = cms.Service("TFileService",
                                 fileName = cms.string("tree.root")
                                  )