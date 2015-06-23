from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update22June2015'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_741_update19June2015/CMSSW_7_4_2_patch1/src/EgammaWork/electron_isolation_CITK_PUPPI.py'
config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.publication = True
#config.Data.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update22June2015'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
