import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HWWValidation")
options = VarParsing.VarParsing ('analysis')

#set default arguments
#options.inputFiles= 'root://eoscms//eos/cms/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2_gedEG-v2/00000/F6EB426C-385B-E311-AF76-02163E008D87.root'
#options.inputFiles= 'root://eoscms//eos/cms/store/relval/CMSSW_7_0_0_pre11/RelValProdTTbar/AODSIM/START70_V4-v1/00000/D0516C65-766A-E311-B744-00259059642E.root'
#options.inputFiles= 'file:/home/users/jgran/CMSSW_7_0_0_pre9/src/HWWValidation/RelVal_CMSSW_7_0_0_pre9.root'
options.inputFiles= 'file:/home/users/jgran/CMSSW_7_0_0_pre11/src/HWWValidation/RelVal_CMSSW_7_0_0_pre11.root'
#options.maxEvents = -1 # -1 means all events
options.maxEvents = 1000 # 

# get and parse the command line arguments
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    #eventsToProcess = cms.untracked.VEventRange('1:2335')
)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.GlobalTag.globaltag = "START62_V1::All"
process.GlobalTag.globaltag = "START70_V4::All"

process.load("HWWValidation.HWWBase.hwwAnalyzer_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('HWWValidation.HWWBase.puJetIDAlgo_cff')

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

process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.dqmsave_step = cms.Path(process.DQMSaver)

process.myPath = cms.Path(process.ak5PFJetsL1FastL2L3*process.ak5PFJetsL1FastL2L3Residual*
                     process.PFJetTracksAssociatorAtVertex*process.PFImpactParameterTagInfos*
                     process.PFTrackCountingHighEffBJetTags*
                     process.hwwAnalyzer
)

process.schedule = cms.Schedule(
                                process.myPath,
                                process.dqmsave_step
)

process.DQMStore.verbose =  cms.untracked.int32(1)
process.DQMStore.verboseQT =  cms.untracked.int32(1)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

