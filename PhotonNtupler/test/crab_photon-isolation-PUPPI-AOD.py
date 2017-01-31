from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'PhotonIsolations_31January2017'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/PUPPI-isolation-moriond-samples-photons/CMSSW_8_0_21/src/EgammaWork/photonIsolationPUPPI_AOD.py'
config.section_("Data")
config.Data.inputDataset = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.publication = False
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
