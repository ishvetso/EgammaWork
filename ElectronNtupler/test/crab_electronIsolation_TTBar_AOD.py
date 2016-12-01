from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'Electron-AllIsolations-CITK-ttbar-1December2016'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_741_update19June2015/CMSSW_7_4_2_patch1/src/EgammaWork/electron_isolation_CITK_AOD.py'
config.section_("Data")
config.Data.inputDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2 
config.Data.publication = False
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
