import FWCore.ParameterSet.Config as cms

gsfTrackMaker = cms.EDProducer(
    "GSFTrackMaker",
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    beamSpotTag       = cms.InputTag("offlineBeamSpot")
)


