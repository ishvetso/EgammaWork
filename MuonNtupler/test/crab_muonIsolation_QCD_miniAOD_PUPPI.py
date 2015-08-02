from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'MuonIsolation_First_QCD'
config.General.workArea = 'crab_projects'
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/i/ishvetso/EgammaWork/MuonIsolation/CMSSW_7_4_6_patch2/src/EgammaWork/muon_isolatin_CITK.py'
config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
config.Data.publication = False
#config.Data.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'MuonIsolation_First_QCD'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
