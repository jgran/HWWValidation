import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *

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
    #eventsToProcess = cms.untracked.VEventRange('1:1764050')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V15::All"

process.load("hww.HWWValidation.eventMaker_cfi");
process.load("hww.HWWValidation.vertexMaker_cfi");
process.load("hww.HWWValidation.trackMaker_cfi");
process.load("hww.HWWValidation.electronMaker_cfi");
process.load("hww.HWWValidation.muonMaker_cfi");
process.load("hww.HWWValidation.fastJetSequence_cff")
process.load("hww.HWWValidation.jetCollectionPruner_cfi");
process.load("hww.HWWValidation.pfJetMaker_cfi");
process.load("hww.HWWValidation.hypDilepMaker_cfi");
process.load("hww.HWWValidation.recoConversionMaker_cfi");
process.load("hww.HWWValidation.scMaker_cfi");
process.load("hww.HWWValidation.pfCandidateMaker_cfi");
process.load("hww.HWWValidation.gsfTrackMaker_cfi");
process.load("hww.HWWValidation.pfElectronMaker_cfi");
process.load("hww.HWWValidation.elToMuAssMaker_cfi");
process.load("hww.HWWValidation.pfElToElAssMaker_cfi");
process.load("hww.HWWValidation.pfmetMaker_cfi");
process.load("hww.HWWValidation.trkMetMaker_cfi");
process.load("hww.HWWValidation.bTagPFJetMaker_cfi");
process.load("hww.HWWValidation.mvaJetIdMaker_cfi");
process.load("hww.HWWValidation.hwwAnalyzer_cfi");

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
process.load('RecoBTag.SoftLepton.softElectronCandProducer_cfi')

process.PFJetTracksAssociatorAtVertex = process.ic5PFJetTracksAssociatorAtVertex.clone()
process.PFJetTracksAssociatorAtVertex.jets = "ak5PFJets"
process.PFJetTracksAssociatorAtVertex.tracks = "generalTracks"
process.PFImpactParameterTagInfos = process.impactParameterTagInfos.clone()
process.PFImpactParameterTagInfos.jetTracks = "PFJetTracksAssociatorAtVertex"
process.PFTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone()
process.PFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("PFImpactParameterTagInfos") )

#from CMGTools.External.puJetIDAlgo_cff import *

process.p = cms.Path(process.eventMaker*process.vertexMaker*process.trackMaker*
                     process.electronMaker*process.muonMaker*process.prunedUncorrectedJets*
                     process.pfJetMaker*process.hypMaker*
                     process.recoConversionMaker*process.scMaker*process.pfCandidateMaker*
                     process.gsfTrackMaker*process.fastJetSequence*process.elToMuAssMaker*
                     process.pfElectronMaker*process.pfElToElAssMaker*process.pfmetMaker*
                     process.trkMetMaker*process.ak5PFJetsL1FastL2L3*process.ak5PFJetsL1FastL2L3Residual*
                     process.PFJetTracksAssociatorAtVertex*process.PFImpactParameterTagInfos*process.PFTrackCountingHighEffBJetTags*process.bTagPFJetMaker*process.mvaJetIdMaker*process.hwwAnalyzer
)
