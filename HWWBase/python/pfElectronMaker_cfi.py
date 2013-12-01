import FWCore.ParameterSet.Config as cms

pfElectronMaker = cms.EDProducer(
    "PFElectronMaker",
    pfCandidatesTag     = cms.InputTag("particleFlow","electrons"),
)
