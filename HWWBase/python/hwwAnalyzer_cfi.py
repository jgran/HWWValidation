import FWCore.ParameterSet.Config as cms

hwwAnalyzer = cms.EDAnalyzer(
    "HWWAnalyzer",
    eventInputTag        = cms.InputTag("eventMaker"),
    vertexInputTag       = cms.InputTag("vertexMaker"),
    trackInputTag        = cms.InputTag("trackMaker"),
    electronInputTag     = cms.InputTag("electronMaker"),
    muonInputTag         = cms.InputTag("muonMaker"),
    hypInputTag          = cms.InputTag("hypMaker"),
    conversionInputTag   = cms.InputTag("recoConversionMaker"),
    scInputTag           = cms.InputTag("scMaker"),
    pfCandidateInputTag  = cms.InputTag("pfCandidateMaker"),
    gsfTrackInputTag     = cms.InputTag("gsfTrackMaker"),
    kt6PFInputTag        = cms.InputTag("kt6PFJetsForEGIsolationRhoMaker"),
    elToMuAssInputTag    = cms.InputTag("elToMuAssMaker"),
    pfElToElAssInputTag  = cms.InputTag("pfElToElAssMaker"),
    fastJetInputTag      = cms.InputTag("fastJetMaker"),
    wwRhoDefaultInputTag = cms.InputTag("wwRhoDefaultMaker"),
    wwRhoVoronoiInputTag = cms.InputTag("wwRhoVoronoiMaker"),
    pfmetInputTag        = cms.InputTag("pfmetMaker"),
    trkMetInputTag       = cms.InputTag("trkMetMaker"),
    pfJetInputTag        = cms.InputTag("pfJetMaker"),
    bTagPFJetInputTag    = cms.InputTag("bTagPFJetMaker"),
    mvaJetIdInputTag     = cms.InputTag("mvaJetIdMaker")
)
