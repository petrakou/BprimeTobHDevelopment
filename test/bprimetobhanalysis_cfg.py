import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

#from BprimebHAnalysis.BprimeTobHAnalysis.bprimetobhanalysis_cfi import *
#from BprimebHAnalysis.BprimeTobHAnalysis.TTJets_Hadronic.TTJets_Hadronic_00_cfi import *
#from BprimebHAnalysis.BprimeTobHAnalysis.BpBpToBHTWinc.BprimeBprimeTobHtWinc_M_1000_cfi import *
from Bprime_kit.BprimeTobHDevelopment.BpBpToBHTWinc.BprimeBprimeTobHtWinc_M_1500_cfi import *
#from BprimebHAnalysis.BprimeTobHAnalysis.BpBpToBHTWinc.BprimeBprimeTobHtWinc_M_800_cfi import *

options = VarParsing('python')

options.register('outFilename', 'file.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('jetPtMin', 50.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum jet Pt"
    )
options.register('jetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet Pt"
    )
options.register('bJetPtMin', 100.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('fatJetPtMin', 150.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet Pt"
    )
options.register('fatJetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet Pt"
    )
options.register('fatJetMassMin', 100.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet mass"
    )
options.register('fatJetMassMax', 150.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet mass"
    )
options.register('fatJetPrunedMassMin', 75.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet pruned mass"
    )
options.register('fatJetPrunedMassMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet pruned mass"
    )
options.register('fatJetTau2ByTau1Min', 0.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet tau2/tau1"
    )
options.register('subjet1CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet1 b discriminator"
    )
options.register('subjet1CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet1 b discriminator"
    )
options.register('subjet2CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet2 b discriminator"
    )
options.register('subjet2CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet2 b discriminator"
    )
options.register('hTMin', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('hTMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum HT"
    )

options.setDefault('maxEvents', -50000) 

options.parseArguments()

process = cms.Process("BprimebHDevelopment") #eleni

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'), 
    ) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

process.BprimebHDevelopment = cms.EDAnalyzer('BprimeTobHDevelopment',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('ntuple/tree'),
#eleni    InputFiles          = cms.vstring(FileNames), 
    InputFiles          = cms.vstring("/tmp/petrakou/BprimeTobH_v1_1_1_ylE.root",
"/tmp/petrakou/BprimeTobH_v1_2_1_UbX.root",
"/tmp/petrakou/BprimeTobH_v1_3_1_3HI.root",
"/tmp/petrakou/BprimeTobH_v1_4_1_fmJ.root",
"/tmp/petrakou/BprimeTobH_v1_5_1_CSW.root",
"/tmp/petrakou/BprimeTobH_v1_6_1_sy4.root",
"/tmp/petrakou/BprimeTobH_v1_7_1_mlE.root",
"/tmp/petrakou/BprimeTobH_v1_8_1_DhQ.root",
"/tmp/petrakou/BprimeTobH_v1_9_1_n73.root",
"/tmp/petrakou/BprimeTobH_v1_10_1_1SL.root"),
#("/tmp/petrakou/BprimeTobH_v1_7_1_mlE.root"),
    JetPtMin            = cms.double(options.jetPtMin),
    JetPtMax            = cms.double(options.jetPtMax),
    JetAbsEtaMax        = cms.double(2.4),
    BJetPtMin           = cms.double(options.bJetPtMin),
    FatJetPtMin         = cms.double(options.fatJetPtMin),
    FatJetPtMax         = cms.double(options.fatJetPtMax),
    FatJetAbsEtaMax     = cms.double(2.4),
    FatJetMassMin       = cms.double(options.fatJetMassMin),
    FatJetMassMax       = cms.double(options.fatJetMassMax),
    FatJetPrunedMassMin = cms.double(options.fatJetPrunedMassMin),
    FatJetPrunedMassMax = cms.double(options.fatJetPrunedMassMax),
    FatJetTau2ByTau1Min = cms.double(options.fatJetTau2ByTau1Min),
    Subjet1CSVDiscMin   = cms.double(options.subjet1CSVDiscMin),
    Subjet1CSVDiscMax   = cms.double(options.subjet1CSVDiscMax),
    Subjet2CSVDiscMin   = cms.double(options.subjet2CSVDiscMin),
    Subjet2CSVDiscMax   = cms.double(options.subjet2CSVDiscMax),
    HTMin               = cms.double(options.hTMin),
    HTMax               = cms.double(options.hTMax), 
    ) 

process.p = cms.Path(process.BprimebHDevelopment)

