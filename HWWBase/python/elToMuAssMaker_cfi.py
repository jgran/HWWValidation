import FWCore.ParameterSet.Config as cms

elToMuAssMaker = cms.EDProducer(
    "ElToMuAssMaker",
    aliasPrefix = cms.untracked.string("els"),
    minDR       = cms.double(0.1),
    elsInputTag = cms.InputTag("electronMaker"),
    musInputTag = cms.InputTag("muonMaker")
)


