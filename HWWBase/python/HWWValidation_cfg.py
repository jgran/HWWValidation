import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HWWValidation")
options = VarParsing.VarParsing ('analysis')

#set default arguments
options.inputFiles= 'root://eoscms//eos/cms/store/relval/CMSSW_6_2_0_pre7_g496p02/RelValProdTTbar/AODSIM/PRE_ST62_V7-v1/00000/CCBB63D5-BCCF-E211-9D8F-002590593878.root'
options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
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
