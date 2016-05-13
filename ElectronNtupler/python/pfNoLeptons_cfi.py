import FWCore.ParameterSet.Config as cms

pfNoLeptons = cms.EDFilter("PdgIdCandViewSelector",
    src = cms.InputTag("packedPFCandidates"),
    pdgId = cms.vint32(211,321,999211,2212,111,130,310,2112,22)
)
