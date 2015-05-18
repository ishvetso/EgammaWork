from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'Photon_Validation_18May'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/Validation_PhotonProducer_74X/CMSSW_7_4_0/src/EgammaWork/photonIsoTest_validation_miniAOD.py'
config.section_("Data")
config.Data.inputDataset = '/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.publication = True
#config.Data.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'Photon_Validation_18May'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
