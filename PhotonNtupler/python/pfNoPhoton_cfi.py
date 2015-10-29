import FWCore.ParameterSet.Config as cms

pfNoPhoton = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    src = cms.InputTag("particleFlow"),
    pdgId = cms.vint32(211,321,999211,2212,111,130,310,2112,11,13),
    makeClones = cms.bool(True)
)
