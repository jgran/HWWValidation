import FWCore.ParameterSet.Config as cms

bTagPFJetMaker = cms.EDProducer(
        "BTagMaker",
        AliasPrefix                  = cms.string("pfjets"),
	      CaloJetsTag                  = cms.InputTag("prunedUncorrectedJets","pfjet"),        
        referenceCaloJetsTag         = cms.InputTag("ak5PFJets"),
	      trackCountingHighEffBJetTags = cms.InputTag("PFTrackCountingHighEffBJetTags"),
)

