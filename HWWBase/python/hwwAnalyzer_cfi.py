import FWCore.ParameterSet.Config as cms
from HWWValidation.HWWBase.puJetIDAlgo_cff import *

hwwAnalyzer = cms.EDAnalyzer(
    "HWWAnalyzer",
    primaryVertexInputTag        = cms.InputTag("offlinePrimaryVertices"),
    trackInputTag                = cms.InputTag("generalTracks"),
    electronsInputTag            = cms.InputTag("gedGsfElectrons"),
    gsftrksInputTag              = cms.InputTag("electronGsfTracks"),
    pfCandsInputTag              = cms.InputTag("particleFlow"),
    recoConversionInputTag       = cms.InputTag("allConversions"),
    beamSpotTag                  = cms.InputTag("offlineBeamSpot"),
    muonsInputTag                = cms.InputTag("muons"),                         
    muonShower                   = cms.InputTag("muons", "muonShowerInformation"),
    pfJetsInputTag               = cms.InputTag("ak5PFJets"),
    pfElectronsTag               = cms.InputTag("particleFlow","electrons"),
    gsftracksInputTag            = cms.InputTag("electronGsfTracks"),
    rhoInputTag                  = cms.InputTag("kt6PFJetsDeterministicIso","rho"),
    wwrhoInputTag                = cms.InputTag("kt6PFJets","rho"),
    wwrhovorInputTag             = cms.InputTag("kt6PFJetsForRhoComputationVoronoi","rho"),
    forEGIsoInputTag             = cms.InputTag("kt6PFJetsForEGIsolation","rho"),
    pfmetInputTag                = cms.InputTag("pfMet"),
    #corrPFJetsInputTag           = cms.InputTag("ak5PFJetsL1FastL2L3"),
    jetCorrector                 = cms.string("ak5PFL1FastL2L3"),
    trackCountingHighEffBJetTags = cms.InputTag("PFTrackCountingHighEffBJetTags"),
    cutBased        = cms.bool(False),
    tmvaWeights   	= cms.string("HWWValidation/HWWBase/data/mva_JetID_v1.weights.xml"),
    tmvaMethod    	= cms.string("JetID"),
    version       	= cms.int32(-1),
    tmvaVariables 	= cms.vstring(
      "nvtx",
      "jetPt",
      "jetEta",
      "jetPhi",
      "dZ",
      "d0",
      "beta",
      "betaStar",
      "nCharged",
      "nNeutrals",
      "dRMean",
      "frac01",
      "frac02",
      "frac03",
      "frac04",
      "frac05"
      ),
    tmvaSpectators = cms.vstring(),
    JetIdParams = JetIdParams,
    InputEGammaWeights = cms.vstring(
      "data/Electrons_BDTG_TrigV0_Cat1.weights.xml", 
      "data/Electrons_BDTG_TrigV0_Cat2.weights.xml", 
      "data/Electrons_BDTG_TrigV0_Cat3.weights.xml", 
      "data/Electrons_BDTG_TrigV0_Cat4.weights.xml", 
      "data/Electrons_BDTG_TrigV0_Cat5.weights.xml", 
      "data/Electrons_BDTG_TrigV0_Cat6.weights.xml"
    ),
    InputMuonIsoWeights = cms.vstring(
      "data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml",
      "data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml",
      "data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml",
      "data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml",
      "data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml",
      "data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml" 
    )
)

