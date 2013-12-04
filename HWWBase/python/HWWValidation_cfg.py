import FWCore.ParameterSet.Config as cms

process = cms.Process("HWWValidation")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/nfs-7/userdata/jaehyeok/A6DE4085-8191-E111-BF4E-001E67396D5B.root'),
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START62_V1::All"

process.load("HWWValidation.HWWBase.eventMaker_cfi");
process.load("HWWValidation.HWWBase.vertexMaker_cfi");
process.load("HWWValidation.HWWBase.trackMaker_cfi");
process.load("HWWValidation.HWWBase.electronMaker_cfi");
process.load("HWWValidation.HWWBase.muonMaker_cfi");
process.load("HWWValidation.HWWBase.fastJetSequence_cff")
process.load("HWWValidation.HWWBase.jetCollectionPruner_cfi");
process.load("HWWValidation.HWWBase.pfJetMaker_cfi");
process.load("HWWValidation.HWWBase.hypDilepMaker_cfi");
process.load("HWWValidation.HWWBase.recoConversionMaker_cfi");
process.load("HWWValidation.HWWBase.scMaker_cfi");
process.load("HWWValidation.HWWBase.pfCandidateMaker_cfi");
process.load("HWWValidation.HWWBase.gsfTrackMaker_cfi");
process.load("HWWValidation.HWWBase.pfElectronMaker_cfi");
process.load("HWWValidation.HWWBase.elToMuAssMaker_cfi");
process.load("HWWValidation.HWWBase.pfElToElAssMaker_cfi");
process.load("HWWValidation.HWWBase.pfmetMaker_cfi");
process.load("HWWValidation.HWWBase.trkMetMaker_cfi");
process.load("HWWValidation.HWWBase.bTagPFJetMaker_cfi");
process.load("HWWValidation.HWWBase.mvaJetIdMaker_cfi");
process.load("HWWValidation.HWWBase.hwwAnalyzer_cfi");

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.ak5PFJetsL1FastL2L3 = process.ak5PFJetsL2L3.clone(
    src = 'ak5PFJets',
    correctors = ['ak5PFL1FastL2L3']
)

process.ak5PFJetsL1FastL2L3Residual = process.ak5PFJetsL2L3.clone(
    src = 'ak5PFJets',
    correctors = ['ak5PFL1FastL2L3Residual']
)

# b-tagging general configuration
process.load('RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi')
process.load('RecoBTag.Configuration.RecoBTag_cff')

process.PFJetTracksAssociatorAtVertex = process.ic5PFJetTracksAssociatorAtVertex.clone()
process.PFJetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.PFJetTracksAssociatorAtVertex.tracks = "generalTracks"
process.PFImpactParameterTagInfos = process.impactParameterTagInfos.clone()
process.PFImpactParameterTagInfos.jetTracks = "PFJetTracksAssociatorAtVertex"
process.PFTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone()
process.PFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("PFImpactParameterTagInfos") )

process.p = cms.Path(process.eventMaker*process.vertexMaker*process.trackMaker*
                     process.electronMaker*process.muonMaker*process.prunedUncorrectedJets*
                     process.pfJetMaker*process.hypMaker*
                     process.recoConversionMaker*process.scMaker*process.pfCandidateMaker*
                     process.gsfTrackMaker*process.fastJetSequence*process.elToMuAssMaker*
                     process.pfElectronMaker*process.pfElToElAssMaker*process.pfmetMaker*
                     process.trkMetMaker*process.ak5PFJetsL1FastL2L3*process.ak5PFJetsL1FastL2L3Residual*
                     process.PFJetTracksAssociatorAtVertex*process.PFImpactParameterTagInfos*process.PFTrackCountingHighEffBJetTags*process.bTagPFJetMaker*process.mvaJetIdMaker*process.hwwAnalyzer
)
