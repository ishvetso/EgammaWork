import FWCore.ParameterSet.Config as cms

pfNoLeptons = cms.EDProducer("MyPdgIDPackedCandidateSelector",
    src = cms.InputTag("packedPFCandidates"),
    pdgId = cms.vint32(11,-11, 13,-13)
)

