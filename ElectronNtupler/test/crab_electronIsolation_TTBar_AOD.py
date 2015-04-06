from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'Electron-Isolation_CITK_validation_ttbar_miniAOD_bugfixed'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/Validation_ElectronIsolation/CMSSW_7_3_0/src/EgammaWork/electron_isolation_CITK.py'
config.section_("Data")
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2 
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.publication = True
#config.Data.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'Electron-Isolation_CITK_validation_ttbar_miniAOD_bugfixed'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
