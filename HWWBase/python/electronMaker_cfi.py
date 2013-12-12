import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDProducer(
    "ElectronMaker",
    aliasPrefix             = cms.untracked.string("els"),
    electronsInputTag       = cms.InputTag("gedGsfElectrons"),
    trksInputTag            = cms.InputTag("generalTracks"),
    gsftracksInputTag       = cms.InputTag("electronGsfTracks"),
    pfCandsInputTag         = cms.InputTag("particleFlow"),
    vtxInputTag             = cms.InputTag("offlinePrimaryVertices"),
    recoConversionInputTag  = cms.InputTag("allConversions"),
    beamSpotTag             = cms.InputTag("offlineBeamSpot"),
    minAbsDist              = cms.double(0.02),        
    minAbsDcot              = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45)
)

