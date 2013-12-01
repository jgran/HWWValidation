import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDProducer(
    "EventMaker",
    aliasPrefix = cms.untracked.string("evt"),
    dcsTag      = cms.InputTag("scalersRawToDigi"),
    isData      = cms.bool(False)
)

