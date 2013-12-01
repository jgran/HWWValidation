import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer(
    "PFJetMaker",
    pfJetsInputTag                   = cms.InputTag("prunedUncorrectedJets", "pfjet"),
    pfCandidatesTag                  = cms.InputTag("particleFlow"),
    pfJetPtCut                       = cms.double(5.),
    PFJetCorrectorL1FastL2L3         = cms.string("ak5PFL1FastL2L3"),
    PFJetCorrectorL1FastL2L3residual = cms.string("ak5PFL1FastL2L3Residual")
)


