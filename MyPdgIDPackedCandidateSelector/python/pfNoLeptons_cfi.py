import FWCore.ParameterSet.Config as cms

pfNoLeptons = cms.EDProducer("MyPdgIDPackedCandidateSelector",
    src = cms.InputTag("packedPFCandidates"),
    pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212, 111,130,310,2112, 22)
)

