import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

####from BpbH.BprimeTobHDevelopment.BpBpToBHBHinc.BprimeBprimeTobHbHinc_M_1000_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'Dec17_bprimeTobH_test.root',
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
options.register('bJetPtMin', 80., #was 100
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('fatJetPtMin', 300., #was 150
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
options.register('fatJetTau2ByTau1Max', 0.5,
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
options.register('doPUReweighting', True,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Do pileup reweighting"
     )
options.register('genDecay', True,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Match H->bb with gen info"
     )
options.register('doNminus1', False,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Do N-1 selection"
     )
options.register('doAnalysis', True,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Do analysis"
     )
options.register('TriggerStudyOn', False,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Do trigger studies"
     )
options.register('MuTriggerStudyOn', False,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Select events passing the muon triggers"
     )
options.register('doTree', False,
     VarParsing.multiplicity.singleton,
     VarParsing.varType.bool,
     "Fill tree"
     )
 
options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("BprimebHDevelopment")

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

MuHLTPaths  = [1830, 1831, 2290, 2650, 2912, 3313, 3590, 4875, 5201, 5671, ] #HLT_Mu40
HLTPaths  = [3225, 4136, 4137, 5089, 5537, 5538] #HLT_HT750
#HLTPaths = [4457, 4458, 4459, 4893, 5703]       #HLT_PFHT750
#HLTPaths = [4469, 4470, 4471, 5222, 5710, 5711] #HLT_PFJet320
#HLTPaths = [3379, 4637, 4638, 5292, 5763]       #HLT_QuadJet90
#HLTPaths = [4213, 4214, 4833, 5115, 5579]       #HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose

process.BprimebHDevelopment = cms.EDAnalyzer('BprimeTobHDevelopment',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('ntuple/tree'),
####    InputFiles          = cms.vstring(FileNames),
    InputFiles          = cms.vstring("root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-800_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_10_1_xsR.root"),
    MuHLTPaths            = cms.vint32(MuHLTPaths), 
    HLTPaths            = cms.vint32(HLTPaths), 
    DoPUReweighting	    = cms.bool(options.doPUReweighting),
    File_PUDistMC	      = cms.string('pileup_Data_Summer12_53X_S10.root'),
    File_PUDistData	    = cms.string('pileup_Data_Summer12_53X_S10.root'),
    Hist_PUDistMC	      = cms.string('pileup_mc'),
    Hist_PUDistData	    = cms.string('pileup_data'),
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
    FatJetTau2ByTau1Max = cms.double(options.fatJetTau2ByTau1Max),
    Subjet1CSVDiscMin   = cms.double(options.subjet1CSVDiscMin),
    Subjet1CSVDiscMax   = cms.double(options.subjet1CSVDiscMax),
    Subjet2CSVDiscMin   = cms.double(options.subjet2CSVDiscMin),
    Subjet2CSVDiscMax   = cms.double(options.subjet2CSVDiscMax),
    HTMin               = cms.double(options.hTMin),
    HTMax               = cms.double(options.hTMax), 
    GenDecay		= cms.bool(options.genDecay),
    DoNminus1           = cms.bool(options.doNminus1),
    DoAnalysis          = cms.bool(options.doAnalysis),
    TriggerStudyOn      = cms.bool(options.TriggerStudyOn),
    MuTriggerStudyOn    = cms.bool(options.MuTriggerStudyOn),
    DoTree              = cms.bool(options.doTree), 
    ) 

process.p = cms.Path(process.BprimebHDevelopment)

