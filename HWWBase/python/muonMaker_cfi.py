import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer(
    "MuonMaker",
    aliasPrefix          = cms.untracked.string("mus"),
    muonsInputTag        = cms.InputTag("muons"        ),                         
    pfCandsInputTag      = cms.InputTag("particleFlow"),
    vtxInputTag          = cms.InputTag("offlinePrimaryVertices"),
    cosmicCompat         = cms.InputTag("muons", "cosmicsVeto"),
    muonShower           = cms.InputTag("muons", "muonShowerInformation"),
    beamSpotTag          = cms.InputTag("offlineBeamSpot"),
    pfNoPileUpInputTag_  = cms.InputTag("pfNoPileUp"),
    tevMuonsName         = cms.string("tevMuons")
)


