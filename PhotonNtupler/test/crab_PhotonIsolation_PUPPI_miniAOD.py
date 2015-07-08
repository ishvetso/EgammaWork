from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'PhotonIsolation_PUPPI_miniAOD_8July2015'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_PUPPI_746_8July/CMSSW_7_4_6_patch2/src/EgammaWork/photonIsolationPUPPI_miniAOD.py'
config.section_("Data")
config.Data.inputDataset = '/GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.publication = True
#config.Data.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'PhotonIsolation_PUPPI_miniAOD_8July2015'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
