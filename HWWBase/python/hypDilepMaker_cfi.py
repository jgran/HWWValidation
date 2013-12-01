
import FWCore.ParameterSet.Config as cms

hypMaker = cms.EDProducer(
    "HypDilepMaker",
	  aliasPrefix       = cms.untracked.string("hyp"),
    electronsInputTag = cms.InputTag("electronMaker"),
    muonsInputTag     = cms.InputTag("muonMaker"),
    jetsInputTag      = cms.InputTag("pfJetMaker","pfjetsp4"),
    LooseLepton_PtCut = cms.double(10.0),
    TightLepton_PtCut = cms.double(20.0),
    hypJetMaxEtaCut   = cms.double(5.0),
    hypJetMinPtCut    = cms.double(30.0)
)


