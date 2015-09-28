from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'PhotonIsolation_checks_17September2015_AOD_fixing_keys'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_final_checks/CMSSW_7_4_7/src/EgammaWork/photonIsolation_AOD.py'
config.section_("Data")
config.Data.inputDataset = '/GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.publication = False
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'