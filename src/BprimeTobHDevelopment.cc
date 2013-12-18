// -*- C++ -*-
//
// Package:    BprimeTobHDevelopment
// Class:      BprimeTobHDevelopment
// 
/**\class BprimeTobHDevelopment BprimeTobHDevelopment.cc Bprime_kit/BprimeTobHDevelopment/src/BprimeTobHDevelopment.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Eleni Petrakou,27 2-020,+41227674870,
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// Second Author:    Devdatta Majumder 
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TEfficiency.h>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "../../BprimeTobH/interface/format.h"
#include "../../BprimeTobH/interface/TriggerBooking.h"
#include "../../BprimeTobH/interface/Njettiness.hh"
#include "../../BprimeTobH/interface/Nsubjettiness.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#include "../../BprimeTobHDevelopment/interface/JetID.h"
#include "../../BprimeTobHDevelopment/interface/histCummulative.h"

//
// class declaration
//

class BprimeTobHDevelopment : public edm::EDAnalyzer {
  public:
    explicit BprimeTobHDevelopment(const edm::ParameterSet&);
    ~BprimeTobHDevelopment();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void CreateHistos(const TString&) ; 
    void AddHisto(const TString&, const TString&,  const TString&, const int&, const double&, const double&) ; 
    template <class Type>
    void FillHisto(const TString& name, const Type value, const double weight);

    edm::LumiReWeighting LumiWeights_;

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    const std::vector<int>          hltPaths_; 
    const std::vector<int>          muHLTPaths_;
    const int                       doPUReweighting_ ;
    const std::string               file_PUDistMC_ ;
    const std::string               file_PUDistData_ ;
    const std::string               hist_PUDistMC_ ;
    const std::string               hist_PUDistData_ ;
    const double jetPtPreselection_ ;
    const double jetPtMin_ ; 
    const double jetPtMax_ ; 
    const double jetAbsEtaMax_ ;
    const double bjetPtMin_ ; 
    const double fatJetPtMin_ ; 
    const double fatJetPtMax_ ; 
    const double fatJetAbsEtaMax_ ;
    const double fatJetMassMin_ ;
    const double fatJetMassMax_ ; 
    const double fatJetPrunedMassMin_ ;
    const double fatJetPrunedMassMax_ ; 
    const double fatJetTau2ByTau1Max_ ; //eleni It was Min till now (?) 
    const double subjet1CSVDiscMin_ ; 
    const double subjet1CSVDiscMax_ ; 
    const double subjet2CSVDiscMin_ ; 
    const double subjet2CSVDiscMax_ ; 
    const double HTMin_ ; 
    const double HTMax_ ; 

    TChain*            chain_;

    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    GenInfoBranches    GenInfo;
    JetInfoBranches    GenJetInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;
    JetInfoBranches    SubJetInfo;
    LepInfoBranches    LepInfo;

    edm::Service<TFileService> fs; 

    bool isData_ ; 
    float evtwt_ ; 
    double puweight_ ;

    bool genDecay_ ;
    bool doNminus1_ ;
    bool doAnalysis_ ;
    bool TriggerStudyOn_ ;
    bool MuTriggerStudyOn_ ;
    bool doTree_ ;

    //// Histograms

    TH1D* h_cutflow ; 

    std::map<TString, TH1D*> hmap_1d ;  

    TGraphAsymmErrors* nHiggsJets_Eff ;
    TGraphAsymmErrors* nBJets_Eff ;
    TGraphAsymmErrors* HT_Eff ;
    TGraphAsymmErrors* pt_HiggsJets_Eff ;
    TGraphAsymmErrors* pt_BJets_Eff ;

    TGraphAsymmErrors* TriggerEff_HT ;
    TGraphAsymmErrors* TriggerEff_AK5HT ;
    TGraphAsymmErrors* TriggerEff_HiggsJetpT ;
    TGraphAsymmErrors* TriggerEff_JetpT ;
    TGraphAsymmErrors* TriggerEff_JetEta ;

    TH1D* h_GenNumHiggsTobb ;
    TH1D* h_DR_GenHjet ;

    TH1D* h_nHiggsJets_NoCuts ;
    TH1D* h_nBJets_NoCuts ;
    TH1D* h_HT_NoCuts ;
    TH1D* h_pt_HiggsJets_NoCuts ;
    TH1D* h_pt_BJets_NoCuts ;

    TH1D* h_nHiggsJets_NoHiggsJetpTCut ;
    TH1D* h_nBJets_NoHiggsJetpTCut ;
    TH1D* h_HT_NoHiggsJetpTCut ;
    TH1D* h_pt_HiggsJets_NoHiggsJetpTCut ;
    TH1D* h_pt_BJets_NoHiggsJetpTCut ;

    TH1D* h_nHiggsJets_NoBJetpTCut ;
    TH1D* h_nBJets_NoBJetpTCut ;
    TH1D* h_HT_NoBJetpTCut ;
    TH1D* h_pt_HiggsJets_NoBJetpTCut ;
    TH1D* h_pt_BJets_NoBJetpTCut ;

    TH1D* h_nHiggsJets_NoNumberHiggsJetsCut ;
    TH1D* h_nBJets_NoNumberHiggsJetsCut ;
    TH1D* h_HT_NoNumberHiggsJetsCut ;
    TH1D* h_pt_HiggsJets_NoNumberHiggsJetsCut ;
    TH1D* h_pt_BJets_NoNumberHiggsJetsCut ;

    TH1D* h_nHiggsJets_NoNumberBJetsCut ;
    TH1D* h_nBJets_NoNumberBJetsCut ;
    TH1D* h_HT_NoNumberBJetsCut ;
    TH1D* h_pt_HiggsJets_NoNumberBJetsCut ;
    TH1D* h_pt_BJets_NoNumberBJetsCut ;

    TH1D* h_nHiggsJets_NoHTCut ;
    TH1D* h_nBJets_NoHTCut ;
    TH1D* h_HT_NoHTCut ;
    TH1D* h_pt_HiggsJets_NoHTCut ;
    TH1D* h_pt_BJets_NoHTCut ;

    //// For trigger efficiencies 
    TH1D* h_TrigStudySel_before_HT ;
    TH1D* h_TrigStudySel_before_AK5HT ;
    TH1D* h_TrigStudySel_before_HiggsJet_pT ;
    TH1D* h_TrigStudySel_before_bJet_pT ;
    TH1D* h_TrigStudySel_before_bJet_Eta ;
    TH1D* h_TrigStudySel_HT ;
    TH1D* h_TrigStudySel_AK5HT ;
    TH1D* h_TrigStudySel_HiggsJet_pT ;
    TH1D* h_TrigStudySel_bJet_pT ;
    TH1D* h_TrigStudySel_bJet_Eta ;

    TH1D* h_varBins_TrigStudySel_before_HT ;
    TH1D* h_varBins_TrigStudySel_before_AK5HT ;
    TH1D* h_varBins_TrigStudySel_HT ;
    TH1D* h_varBins_TrigStudySel_AK5HT ;
    TH1D* h_varBins_TrigStudySel_failHLT_HT ;
    TH1D* h_varBins_TrigStudySel_failHLT_AK5HT ;

    TH2D* h2_HT_AK5HT;
    TH2D* h2_HT_AK5HT_passHLT;
    TH2D* h2_HT_AK5HT_failHLT;

    //// The tree
    TTree *tree;
    float final_nHiggsJets, final_nBJets, final_HiggsJets_pt, final_BJets_pt, final_HT, final_evtwt, final_HiggsJets_eta, final_BJets_eta;
    float final_AK5HT, final_HiggsJets_phi, final_BJets_phi, final_Bprime_mass, final_DR_bHiggs;
    float final_HJet_Mass, final_HJet_PrunedMass, final_HJet_Nsubjetiness, final_Subjet_MaxPt, final_Subjet_MinPt, final_Subjet_CSV_MaxPt, final_Subjet_CSV_MinPt, final_Subjet_MaxCSV, final_Subjet_MinCSV, final_Subjet_DR, final_bJet_CSV;

};

//
// constructors and destructor
//
BprimeTobHDevelopment::BprimeTobHDevelopment(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
  muHLTPaths_(iConfig.getParameter<std::vector<int> >("MuHLTPaths")),
  hltPaths_(iConfig.getParameter<std::vector<int> >("HLTPaths")),
  doPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")),
  file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
  file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
  hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
  hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),
  jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
  jetPtMax_(iConfig.getParameter<double>("JetPtMax")),
  jetAbsEtaMax_(iConfig.getParameter<double>("JetAbsEtaMax")),
  bjetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
  fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
  fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")),
  fatJetAbsEtaMax_(iConfig.getParameter<double>("FatJetAbsEtaMax")),
  fatJetMassMin_(iConfig.getParameter<double>("FatJetMassMin")),
  fatJetMassMax_(iConfig.getParameter<double>("FatJetMassMax")), 
  fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
  fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
  fatJetTau2ByTau1Max_(iConfig.getParameter<double>("FatJetTau2ByTau1Max")),
  subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
  subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
  subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
  subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
  HTMin_(iConfig.getParameter<double>("HTMin")), 
  HTMax_(iConfig.getParameter<double>("HTMax")),
  isData_(0),
  evtwt_(1),
  puweight_(1),
  jetPtPreselection_(30.),
  genDecay_(iConfig.getParameter<bool>("GenDecay")), // Note: This applies to any other set flag. 
  doNminus1_(iConfig.getParameter<bool>("DoNminus1")),
  doAnalysis_(iConfig.getParameter<bool>("DoAnalysis")),
  TriggerStudyOn_(iConfig.getParameter<bool>("TriggerStudyOn")), // Note: For this to be meaningful, the doAnalysis_ has to be set as well.
								 // NOTE: When this is set, then the trigger bin in cutflow won't get filled.
  MuTriggerStudyOn_(iConfig.getParameter<bool>("MuTriggerStudyOn")), // Note: This has an effect only when TriggerStudyOn_ is set. 
  doTree_(iConfig.getParameter<bool>("DoTree"))
{ 

   if (doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_) ;
}


BprimeTobHDevelopment::~BprimeTobHDevelopment() { 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BprimeTobHDevelopment::beginJob() { 

  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());

    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    f->Close();
  }

  EvtInfo.Register(chain_);
  VtxInfo.Register(chain_);
  GenInfo.Register(chain_);
  GenJetInfo.Register(chain_,"GenJetInfo");
  JetInfo.Register(chain_,"JetInfo");
  FatJetInfo.Register(chain_,"FatJetInfo");
  SubJetInfo.Register(chain_,"SubJetInfo");
  LepInfo.Register(chain_);

  if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  TH1::SetDefaultSumw2(kTRUE);

// NOTE: Devdatta's binning below: 4000 and 2000 instead of 200 and 100 bins. 

  h_GenNumHiggsTobb	       = fs->make<TH1D>("h_nGenHiggsTobb"             ,"N_{Gen HiggsTobb}"             ,5  ,-0.5,4.5 );
  h_DR_GenHjet		       = fs->make<TH1D>("h_DR_GenHjet"                ,"DR_{GenH,Hjet}"                ,100 ,0.  ,5.);

  h_nHiggsJets_NoCuts          = fs->make<TH1D>("h_nHJets_NoCuts"                  ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoCuts              = fs->make<TH1D>("h_nBJets_NoCuts"                  ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoCuts                  = fs->make<TH1D>("h_HT_NoCuts"                      ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoCuts        = fs->make<TH1D>("h_HiggsJet_Pt_NoCuts"             ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoCuts            = fs->make<TH1D>("h_BJet_Pt_NoCuts"                 ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

  h_nHiggsJets_NoHiggsJetpTCut  = fs->make<TH1D>("h_nHJets_NoHiggsJetpTCut"         ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoHiggsJetpTCut      = fs->make<TH1D>("h_nBJets_NoHiggsJetpTCut"         ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoHiggsJetpTCut          = fs->make<TH1D>("h_HT_NoHiggsJetpTCut"             ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoHiggsJetpTCut = fs->make<TH1D>("h_HiggsJet_Pt_NoHiggsJetpTCut"   ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoHiggsJetpTCut    = fs->make<TH1D>("h_BJet_Pt_NoHiggsJetpTCut"        ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

  h_nHiggsJets_NoBJetpTCut      = fs->make<TH1D>("h_nHJets_NoBJetpTCut"             ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoBJetpTCut          = fs->make<TH1D>("h_nBJets_NoBJetpTCut"             ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoBJetpTCut              = fs->make<TH1D>("h_HT_NoBJetpTCut"                 ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoBJetpTCut    = fs->make<TH1D>("h_HiggsJet_Pt_NoBJetpTCut"        ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoBJetpTCut        = fs->make<TH1D>("h_BJet_Pt_NoBJetpTCut"            ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

  h_nHiggsJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_nHJets_NoNumberHiggsJetsCut"        ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_nBJets_NoNumberHiggsJetsCut"            ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoNumberHiggsJetsCut     = fs->make<TH1D>("h_HT_NoNumberHiggsJetsCut"                ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_HiggsJet_Pt_NoNumberHiggsJetsCut" ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_BJet_Pt_NoNumberHiggsJetsCut"         ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

  h_nHiggsJets_NoNumberBJetsCut = fs->make<TH1D>("h_nHJets_NoNumberBJetsCut"                ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoNumberBJetsCut     = fs->make<TH1D>("h_nBJets_NoNumberBJetsCut"                ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoNumberBJetsCut         = fs->make<TH1D>("h_HT_NoNumberBJetsCut"                    ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoNumberBJetsCut = fs->make<TH1D>("h_HiggsJet_Pt_NoNumberBJetsCut"         ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoNumberBJetsCut   = fs->make<TH1D>("h_BJet_Pt_NoNumberBJetsCut"               ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

  h_nHiggsJets_NoHTCut          = fs->make<TH1D>("h_nHJets_NoHTCut"                  ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoHTCut              = fs->make<TH1D>("h_nBJets_NoHTCut"                  ,"N_{b jets}"                 ,21  ,-0.5,20.5 );
  h_HT_NoHTCut                  = fs->make<TH1D>("h_HT_NoHTCut"                      ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoHTCut        = fs->make<TH1D>("h_HiggsJet_Pt_NoHTCut"             ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoHTCut            = fs->make<TH1D>("h_BJet_Pt_NoHTCut"                 ,"b jet p_{T} [GeV]"          ,100 ,0.  ,2000.);

// NOTE: Devdatta's binning below: 4000 and 2000 instead of 100 and 50 bins (except for eta). 

  h_TrigStudySel_before_HT	= fs->make<TH1D>("h_TrigStudySel_before_HT"                     ,"H_{T}[GeV]"                 ,100 ,0.  ,4000.);
  h_TrigStudySel_before_AK5HT      = fs->make<TH1D>("h_TrigStudySel_before_AK5HT"                     ,"H_{T}[GeV]"                 ,100 ,0.  ,4000.);
  h_TrigStudySel_before_HiggsJet_pT = fs->make<TH1D>("h_TrigStudySel_before_HiggsJet_pT"        ,"Higgs jet p_{T} [GeV]"      ,50 ,0.  ,2000.);
  h_TrigStudySel_before_bJet_pT = fs->make<TH1D>("h_TrigStudySel_before_bJet_pT"     ,"b jet p_{T} [GeV]"          ,50 ,0.  ,2000.);
  h_TrigStudySel_before_bJet_Eta = fs->make<TH1D>("h_TrigStudySel_before_bJet_Eta"     ,"b jet #eta"          ,50 ,-2.5  ,2.5);
  h_TrigStudySel_HT		= fs->make<TH1D>("h_TrigStudySel_HT"                 ,"H_{T}[GeV]"                 ,100 ,0.  ,4000.);
  h_TrigStudySel_AK5HT             = fs->make<TH1D>("h_TrigStudySel_AK5HT"                 ,"H_{T}[GeV]"                 ,100 ,0.  ,4000.);
  h_TrigStudySel_HiggsJet_pT	= fs->make<TH1D>("h_TrigStudySel_HiggsJet_pT"        ,"Higgs jet p_{T} [GeV]"      ,50 ,0.  ,2000.);
  h_TrigStudySel_bJet_pT	= fs->make<TH1D>("h_TrigStudySel_bJet_pT"            ,"b jet p_{T} [GeV]"          ,50 ,0.  ,2000.);
  h_TrigStudySel_bJet_Eta = fs->make<TH1D>("h_TrigStudySel_bJet_Eta"     ,"b jet #eta"          ,50 ,-2.5  ,2.5);

  const int nbins(12) ;
  const double binEdges[nbins+1] = {0., 100., 200, 300., 400., 500., 600., 800., 1000., 1500., 2000., 3000., 4000.} ;
  h_varBins_TrigStudySel_before_HT    = fs->make<TH1D>("h_varBins_TrigStudySel_before_HT"    ,"HT without HLT; H_{T}(Higgs + b jets) [GeV]"       ,nbins ,binEdges) ;
  h_varBins_TrigStudySel_before_AK5HT = fs->make<TH1D>("h_varBins_TrigStudySel_before_AK5HT" ,"HT(AK5 jets) without HLT; H_{T}(AK5 jets) [GeV]" ,nbins ,binEdges) ;
  h_varBins_TrigStudySel_HT           = fs->make<TH1D>("h_varBins_TrigStudySel_HT"           ,"HT without HLT; H_{T}(Higgs + b jets) [GeV]"       ,nbins ,binEdges) ;
  h_varBins_TrigStudySel_AK5HT        = fs->make<TH1D>("h_varBins_TrigStudySel_AK5HT"        ,"HT(AK5 jets) without HLT; H_{T}(AK5 jets) [GeV]" ,nbins ,binEdges) ;
  h_varBins_TrigStudySel_failHLT_HT    = fs->make<TH1D>("h_varBins_TrigStudySel_failHLT_HT"    ,"HT without HLT; H_{T}(Higgs + b jets) [GeV]"       ,nbins ,binEdges) ;
  h_varBins_TrigStudySel_failHLT_AK5HT = fs->make<TH1D>("h_varBins_TrigStudySel_failHLT_AK5HT" ,"HT(AK5 jets) without HLT; H_{T}(AK5 jets) [GeV]" ,nbins ,binEdges) ;

  h2_HT_AK5HT = fs->make<TH2D>("h2_HT_AK5HT", "H_{T} (Higgs + b jets) vs. H_{T} (AK5 jets); H_{T} (Higgs + b jets) [GeV]; H_{T} (AK5 jets) [GeV]; Events/50 GeV" , 80, 0., 4000., 80, 0., 4000) ;
  h2_HT_AK5HT_passHLT = fs->make<TH2D>("h2_HT_AK5HT_passHLT", "H_{T} (Higgs + b jets) vs. H_{T} (AK5 jets) (Events passing HLT); H_{T} (Higgs + b jets) [GeV]; H_{T} (AK5 jets) [GeV]; Events/50 GeV" , 80, 0., 4000., 80, 0., 4000) ;
  h2_HT_AK5HT_failHLT = fs->make<TH2D>("h2_HT_AK5HT_failHLT", "H_{T} (Higgs + b jets) vs. H_{T} (AK5 jets) (Events failing HLT); H_{T} (Higgs + b jets) [GeV]; H_{T} (AK5 jets) [GeV]; Events/50 GeV" , 80, 0., 4000., 80, 0., 4000) ;

//  h_HiggsPt                    = fs->make<TH1D>("h_HiggsPt"                   ,"Higgs p_{T} [GeV]"          ,100 ,0.  ,2000.);
//  h_HiggsPtMatchedJet          = fs->make<TH1D>("h_HiggsPtMatchedJet"         ,"Matched Higgs p_{T} [GeV]"  ,100 ,0.  ,2000.);
  h_cutflow                    = fs->make<TH1D>("h_cutflow"                   ,"Cut flow"                   ,20  ,0.  ,20.  ); 

//  teff_HiggsJetMatch           = fs->make<TEfficiency>("teff_HiggsJetMatch" ,"Higgs-Jet matching efficiency" ,100 ,0.  ,2000.) ; 

//  h_HiggsPt                    -> Sumw2() ; 
//  h_HiggsPtMatchedJet          -> Sumw2() ; 

//eleni?  h_cutflow -> Sumw2() ; 

  h_cutflow -> GetXaxis() -> SetBinLabel(1,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2,"TriggerSel") ; 
/*eleni  h_cutflow -> GetXaxis() -> SetBinLabel(3,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5,"HTSel") ; */
  h_cutflow -> GetXaxis() -> SetBinLabel(3,"FatJetSel") ;
  h_cutflow -> GetXaxis() -> SetBinLabel(4,"HiggsJetSel") ;
  h_cutflow -> GetXaxis() -> SetBinLabel(5,"JetSel") ;
  h_cutflow -> GetXaxis() -> SetBinLabel(6,"BJetsSel") ;
  h_cutflow -> GetXaxis() -> SetBinLabel(7,"HTSel") ;
  h_cutflow -> GetXaxis() -> SetBinLabel(8,"TriggerStudy") ;

  for (int ii = 1; ii <= 8; ++ii) 
    CreateHistos(h_cutflow->GetXaxis()->GetBinLabel(ii)) ;

  CreateHistos("TrigStudySel"); //eleni
 
/* eleni?
  h_GenNumHiggsTobb -> Sumw2() ;
  h_DR_GenHjet      -> Sumw2() ;

  h_nHiggsJets_NoCuts -> Sumw2() ;
  h_nBJets_NoCuts -> Sumw2() ;
  h_HT_NoCuts -> Sumw2() ;
  h_pt_HiggsJets_NoCuts -> Sumw2() ;
  h_pt_BJets_NoCuts -> Sumw2() ;

  h_nHiggsJets_NoHiggsJetpTCut -> Sumw2() ;
  h_nBJets_NoHiggsJetpTCut -> Sumw2() ;
  h_HT_NoHiggsJetpTCut -> Sumw2() ;
  h_pt_HiggsJets_NoHiggsJetpTCut -> Sumw2() ;
  h_pt_BJets_NoHiggsJetpTCut -> Sumw2() ;

  h_nHiggsJets_NoBJetpTCut -> Sumw2() ;
  h_nBJets_NoBJetpTCut -> Sumw2() ;
  h_HT_NoBJetpTCut -> Sumw2() ;
  h_pt_HiggsJets_NoBJetpTCut -> Sumw2() ;
  h_pt_BJets_NoBJetpTCut -> Sumw2() ;

  h_nHiggsJets_NoNumberHiggsJetsCut -> Sumw2() ;
  h_nBJets_NoNumberHiggsJetsCut -> Sumw2() ;
  h_HT_NoNumberHiggsJetsCut -> Sumw2() ;
  h_pt_HiggsJets_NoNumberHiggsJetsCut -> Sumw2() ;
  h_pt_BJets_NoNumberHiggsJetsCut -> Sumw2() ;

  h_nHiggsJets_NoNumberBJetsCut -> Sumw2() ;
  h_nBJets_NoNumberBJetsCut -> Sumw2() ;
  h_HT_NoNumberBJetsCut -> Sumw2() ;
  h_pt_HiggsJets_NoNumberBJetsCut -> Sumw2() ;
  h_pt_BJets_NoNumberBJetsCut -> Sumw2() ;

  h_nHiggsJets_NoHTCut -> Sumw2() ;
  h_nBJets_NoHTCut -> Sumw2() ;
  h_HT_NoHTCut -> Sumw2() ;
  h_pt_HiggsJets_NoHTCut -> Sumw2() ;
  h_pt_BJets_NoHTCut -> Sumw2() ;

  h_TrigStudySel_before_HT -> Sumw2() ;
  h_TrigStudySel_before_AK5HT -> Sumw2() ;
  h_TrigStudySel_before_HiggsJet_pT -> Sumw2() ;
  h_TrigStudySel_before_bJet_pT -> Sumw2() ;
  h_TrigStudySel_before_bJet_Eta -> Sumw2() ;
  h_TrigStudySel_HT -> Sumw2() ;
  h_TrigStudySel_AK5HT -> Sumw2() ;
  h_TrigStudySel_HiggsJet_pT -> Sumw2() ;
  h_TrigStudySel_bJet_pT -> Sumw2() ;
  h_TrigStudySel_bJet_Eta -> Sumw2() ;
*/

   tree = new TTree("T","Branches for TMVA use");

   // Event kinematics
   tree->Branch("final_evtwt", &final_evtwt, "final_evtwt/F");
   tree->Branch("final_nHiggsJets", &final_nHiggsJets, "final_nHiggsJets/F");
   tree->Branch("final_nBJets", &final_nBJets, "final_nBJets/F");
   tree->Branch("final_HT", &final_HT, "final_HT/F");
   tree->Branch("final_AK5HT", &final_AK5HT, "final_AK5HT/F");
   tree->Branch("final_HiggsJets_pt", &final_HiggsJets_pt, "final_HiggsJets_pt/F");
   tree->Branch("final_HiggsJets_eta", &final_HiggsJets_eta, "final_HiggsJets_eta/F");
   tree->Branch("final_HiggsJets_phi", &final_HiggsJets_phi, "final_HiggsJets_phi/F");
   tree->Branch("final_BJets_pt", &final_BJets_pt, "final_BJets_pt/F");
   tree->Branch("final_BJets_eta", &final_BJets_eta, "final_BJets_eta/F");
   tree->Branch("final_BJets_phi", &final_BJets_phi, "final_BJets_phi/F");
   tree->Branch("final_DR_bHiggs", &final_DR_bHiggs, "final_DR_bHiggs/F");
   tree->Branch("final_Bprime_mass", &final_Bprime_mass, "final_Bprime_mass/F");

   // Subjets, b-taggings and the like
   tree->Branch("final_HJet_Mass", &final_HJet_Mass, "final_HJet_Mass/F");
   tree->Branch("final_HJet_PrunedMass", &final_HJet_PrunedMass, "final_HJet_PrunedMass/F");
   tree->Branch("final_HJet_Nsubjetiness", &final_HJet_Nsubjetiness, "final_HJet_Nsubjetiness/F");
   tree->Branch("final_Subjet_MaxPt", &final_Subjet_MaxPt, "final_Subjet_MaxPt/F");
   tree->Branch("final_Subjet_MinPt", &final_Subjet_MinPt, "final_Subjet_MinPt/F");
   tree->Branch("final_Subjet_CSV_MaxPt", &final_Subjet_CSV_MaxPt, "final_Subjet_CSV_MaxPt/F");
   tree->Branch("final_Subjet_CSV_MinPt", &final_Subjet_CSV_MinPt, "final_Subjet_CSV_MinPt/F");
   tree->Branch("final_Subjet_MaxCSV", &final_Subjet_MaxCSV, "final_Subjet_MaxCSV/F");
   tree->Branch("final_Subjet_MinCSV", &final_Subjet_MinCSV, "final_Subjet_MinCSV/F");
   tree->Branch("final_Subjet_DR", &final_Subjet_DR, "final_Subjet_DR/F");
   tree->Branch("final_bJet_CSV", &final_bJet_CSV, "final_bJet_CSV/F");

   // Event shape and other advanced coolness
   //tree->Branch("", &, "/F");

   //tree->Branch("", &, "/F");

  return ;  

}

void BprimeTobHDevelopment::CreateHistos(const TString& cutname) {

   AddHisto(cutname ,"_nPVtx_NoPUWt"               ,"N(PV), No PU weight"       ,50     ,-0.5     ,49.5    ) ;
   AddHisto(cutname ,"_nPVtx_PUWt"                 ,"N(PV)"                     ,50     ,-0.5     ,49.5    ) ;
   AddHisto(cutname ,"_nJets"                      ,"N(AK5 jets)"               ,20     ,-0.5     ,19.5    ) ;
   AddHisto(cutname ,"_nBJets"                     ,"N(b jets)"                 ,20     ,-0.5     ,19.5    ) ;
   AddHisto(cutname ,"_nFatJets"                   ,"N(fat jets)"               ,20     ,-0.5     ,19.5    ) ;
   AddHisto(cutname ,"_nHJets"                     ,"N(Higgs jets)"             ,20     ,-0.5     ,19.5    ) ;
   AddHisto(cutname ,"_FatJets_Pt"                 ,"p_{T}(fat jets)"           ,1000   ,0.       ,1000.   ) ;
   AddHisto(cutname ,"_FatJets_Eta"                ,"#eta(fat jets)"            ,50     ,-4.      ,4.      ) ;
   AddHisto(cutname ,"_FatJets_Mass"               ,"Fat jet mass [GeV]"        ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_FatJets_MassPruned"         ,"Fat jet pruned mass [GeV]" ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_FatJets_tau2ByTau1"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ;
   AddHisto(cutname ,"_FatJets_tau3ByTau2"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ;
   AddHisto(cutname ,"_FatJets_tau3ByTau1"         ,"Fat jet #tau2/#tau1"       ,20     ,0.       ,1.      ) ;
   AddHisto(cutname ,"_FatJets_CombinedSVBJetTags" ,"Fat jet CSV discriminator" ,20     ,0.       ,1.      ) ;

   AddHisto(cutname ,"_SubJet1_Pt"                 ,"SubJet1 p_{T} [GeV]"       ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_SubJet1_Eta"                ,"SubJet1 #eta"              ,50     ,-4.      ,4.      ) ;
   AddHisto(cutname ,"_SubJet1_Mass"               ,"SubJet1 mass [GeV]"        ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_SubJet1_CombinedSVBJetTags" ,"SubJet1 CSV discriminator" ,20     ,0.       ,1.      ) ;

   AddHisto(cutname ,"_SubJet2_Pt"                 ,"SubJet2 p_{T} [GeV]"       ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_SubJet2_Eta"                ,"SubJet2 #eta"              ,50     ,-4.      ,4.      ) ;
   AddHisto(cutname ,"_SubJet2_Mass"               ,"SubJet2 mass [GeV]"        ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_SubJet2_CombinedSVBJetTags" ,"SubJet2 CSV discriminator" ,20     ,0.       ,1.      ) ;

   AddHisto(cutname ,"_HiggsJet_Pt"                ,"p_{T} (Higgs jet)[GeV]"    ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_HiggsJet_Eta"               ,"#eta (Higgs jet)"          ,50     ,-4.      ,4.      ) ;
   AddHisto(cutname ,"_HiggsJet_Mass"              ,"Mass (Higgs jet) [GeV]"    ,100    ,0.       ,200.    ) ;

   AddHisto(cutname ,"_BJet_Pt"                    ,"p_{T} (b jet)[GeV]"        ,100    ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_BJet_Eta"                   ,"#eta (b jet)"              ,50     ,-4.      ,4.      ) ;

   AddHisto(cutname ,"_HT"                         ,"H_{T}[GeV]"                 ,200   ,0.       ,4000.   ) ;
   AddHisto(cutname ,"_AK5HT"                      ,"H_{T}R (AK5 jets) [GeV]"    ,200   ,0.       ,4000.   ) ;
   AddHisto(cutname ,"_bprimePt"                   ,"b' p_{T} [GeV]"             ,100   ,0.       ,2000.   ) ;
   AddHisto(cutname ,"_bprimeMass"                 ,"b' mass [GeV]"              ,40    ,0.       ,2000.   ) ;

  return ; 
}

void BprimeTobHDevelopment::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbins, const double& min, const double& max) { 

  TH1D* h1d ;
  h1d = fs->make<TH1D>(cutname+histname, cutname+histtitle, nbins, min, max);
  h1d -> Sumw2() ;
  hmap_1d[cutname+histname] = h1d ;

  return ; 

}

template <class Type>
void BprimeTobHDevelopment::FillHisto(const TString& name, const Type value, const double weight){

  hmap_1d[name]->Fill(double(value),weight);

  return ; 

}


// ------------ method called for each event  ------------
void BprimeTobHDevelopment::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  JetID jetIDTight(JetID::FIRSTDATA,JetID::TIGHT, JetInfo) ; 
  JetID fatjetIDLoose(JetID::FIRSTDATA,JetID::LOOSE, FatJetInfo) ;
  pat::strbitset retak5 = jetIDTight.getBitTemplate() ;
  pat::strbitset retca8 = fatjetIDLoose.getBitTemplate() ;

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  if ( doNminus1_ ) {
  for(int entry=0; entry<maxEvents_; entry++) {
//  for(int entry=0; entry<1000; entry++) {

    //// Event variables 
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector> bprimes ; 

    int njets(0) ; 
    int  nGoodVtxs(0) ;

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << "A- " << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1;
    if ( !isData_ ) evtwt_    = EvtInfo.Weight ;
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ;
    evtwt_ *= puweight_ ;

    //// Select good vertices 
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { cout << endl << "Vtx!" << endl; continue ; }

    //// Select events with at least one H->bb decay 
    if( genDecay_ ) {
        int iHbb(0);
        for( int igen=0; igen<GenInfo.Size; igen++){
            if( GenInfo.Status[igen] != 3 || GenInfo.PdgID[igen] != 25 || GenInfo.nDa[igen] < 2) continue;
            if( abs(GenInfo.Da0PdgID[igen])!=5 ) continue;
	    if( GenInfo.Da0PdgID[igen]/GenInfo.Da1PdgID[igen]!=-1 ) continue;
            iHbb++;
        }
        if( iHbb==0 ) continue;
    }

    bool passHLT(false);
    for ( std::vector<int>::const_iterator ihlt = hltPaths_.begin(); ihlt != hltPaths_.end(); ++ihlt ) {
      if (EvtInfo.TrgBook[*ihlt] == 1) { 
         passHLT = true ;
         break ; 
      }
      else passHLT = false ;
    }
    
    if ( !passHLT ) continue;

    bool CA8pre(false);
    for (int ifj=0; ifj < FatJetInfo.Size; ++ifj) {
      if ( FatJetInfo.Pt[ifj] > 100. ) {
	CA8pre = true;
	break;
      }
    }
    if( !CA8pre ) continue;

    int AK5preCount(0);
    for (int ij=0; ij < JetInfo.Size; ++ij) {
      if ( JetInfo.Pt[ij] > 30. ) AK5preCount++;
    }
    if( AK5preCount < 2 ) continue;


    // 0. NO CUTS ------------------------- 

    double HT0(0) ;
    std::vector<TLorentzVector>higgsJets0 ;
    std::vector<TLorentzVector>bJets0 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
//      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) { }
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ )
        continue; //// apply jet eta cut
//      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 )
//        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ )
        continue; //// apply pruned jet mass cut
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet],
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. )
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1],
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2],
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) )
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_
          && fatjet_p4.Mag() < fatJetMassMax_
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) {

        higgsJets0.push_back(fatjet_p4) ;

      }

    } //// Loop over fat jets


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) {

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) { }
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ;
//      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ;
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ;

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet],
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ;
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets0.begin(); ihig != higgsJets0.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ;
          break ;
        }
        else {
          isJetNotHiggs = true ;
        }
      }
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet],
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets0.push_back(bjet_p4) ;

      } //// Select b-tagged AK5 jets

    } //// Loop over AK5 jets

    for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets0.begin(); ihig != higgsJets0.end(); ++ihig) {
      HT0 += ihig->Pt() ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib = bJets0.begin(); ib != bJets0.end(); ++ib) {
      HT0 += ib->Pt() ;
    }

    h_nHiggsJets_NoCuts->Fill(higgsJets0.size(),evtwt_) ;
    h_nBJets_NoCuts->Fill(bJets0.size(),evtwt_) ;
    h_HT_NoCuts->Fill(HT0,evtwt_) ;

    if ( higgsJets0.size()!=0 && bJets0.size()!=0 ) {
      for (std::vector<TLorentzVector>::const_iterator ihig_ = higgsJets0.begin(); ihig_ != higgsJets0.end(); ++ihig_) {
        h_pt_HiggsJets_NoCuts->Fill(ihig_->Pt(),evtwt_) ;
      }
      for (std::vector<TLorentzVector>::const_iterator ib_ = bJets0.begin(); ib_ != bJets0.end(); ++ib_) {
        h_pt_BJets_NoCuts->Fill(ib_->Pt(),evtwt_) ;
      }
    }
    njets = 0;


    // 1. NO HIGGS JETS PT CUT ---------------------------------------- 

    double HT1(0) ;
    std::vector<TLorentzVector>higgsJets1 ;
    std::vector<TLorentzVector>bJets1 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
//      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) { }
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets1.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets1.begin(); ihig != higgsJets1.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets1.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


    if (higgsJets1.size() >= 1) {

      if (bJets1.size() >= 2 ) { 

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets1.begin(); ihig != higgsJets1.end(); ++ihig) { 
          HT1 += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets1.begin(); ib != bJets1.end(); ++ib) { 
          HT1 += ib->Pt() ; 
        }

	if ( HT1 > HTMin_ && HT1 < HTMax_ ) {
	  h_nHiggsJets_NoHiggsJetpTCut->Fill(higgsJets1.size(),evtwt_) ;
	  h_nBJets_NoHiggsJetpTCut->Fill(bJets1.size(),evtwt_) ;
	  h_HT_NoHiggsJetpTCut->Fill(HT1,evtwt_) ;

	  for (std::vector<TLorentzVector>::const_iterator ihig_ = higgsJets1.begin(); ihig_ != higgsJets1.end(); ++ihig_) {
	    h_pt_HiggsJets_NoHiggsJetpTCut->Fill(ihig_->Pt(),evtwt_) ;
	  }
	  for (std::vector<TLorentzVector>::const_iterator ib_ = bJets1.begin(); ib_ != bJets1.end(); ++ib_) {
	    h_pt_BJets_NoHiggsJetpTCut->Fill(ib_->Pt(),evtwt_) ;
	  }
	}

      } //// If at least two b-jets
    } //// If at least one Higgs jet
    njets = 0;


    // 2. NO B-JETS PT CUT ---------------------------------------- 

    double HT2(0) ;
    std::vector<TLorentzVector>higgsJets2 ;
    std::vector<TLorentzVector>bJets2 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ )
	continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets2.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
//      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) {}
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets2.begin(); ihig != higgsJets2.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets2.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


    if (higgsJets2.size() >= 1) {

      if (bJets2.size() >= 2 ) { 

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets2.begin(); ihig != higgsJets2.end(); ++ihig) { 
          HT2 += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets2.begin(); ib != bJets2.end(); ++ib) { 
          HT2 += ib->Pt() ; 
        }

  if ( HT1 > HTMin_ && HT1 < HTMax_ ) {
    h_nHiggsJets_NoBJetpTCut->Fill(higgsJets2.size(),evtwt_) ;
    h_nBJets_NoBJetpTCut->Fill(bJets2.size(),evtwt_) ;
    h_HT_NoBJetpTCut->Fill(HT2,evtwt_) ;

    for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets2.begin(); ihig != higgsJets2.end(); ++ihig) {
        h_pt_HiggsJets_NoBJetpTCut->Fill(ihig->Pt(),evtwt_) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib = bJets2.begin(); ib != bJets2.end(); ++ib) {
        h_pt_BJets_NoBJetpTCut->Fill(ib->Pt(),evtwt_) ;
    }
  }
      } //// If at least two b-jets
    } //// If at least one Higgs jet
    njets = 0;


    // 3. NO CUT ON NUMBER OF HIGGS JETS ---------------------------------------- 

    double HT3(0) ;
    std::vector<TLorentzVector>higgsJets3 ;
    std::vector<TLorentzVector>bJets3 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
	continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets3.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets3.begin(); ihig != higgsJets3.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets3.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


      if (bJets3.size() >= 2 ) { 

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets3.begin(); ihig != higgsJets3.end(); ++ihig) { 
          HT3 += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets3.begin(); ib != bJets3.end(); ++ib) { 
          HT3 += ib->Pt() ; 
        }

  if ( HT1 > HTMin_ && HT1 < HTMax_ ) {
    h_nHiggsJets_NoNumberHiggsJetsCut->Fill(higgsJets3.size(),evtwt_) ;
    h_nBJets_NoNumberHiggsJetsCut->Fill(bJets3.size(),evtwt_) ;
    h_HT_NoNumberHiggsJetsCut->Fill(HT3,evtwt_) ;

    for (std::vector<TLorentzVector>::const_iterator ihig__ = higgsJets3.begin(); ihig__ != higgsJets3.end(); ++ihig__) {
        h_pt_HiggsJets_NoNumberHiggsJetsCut->Fill(ihig__->Pt(),evtwt_) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib__ = bJets3.begin(); ib__ != bJets3.end(); ++ib__) {
        h_pt_BJets_NoNumberHiggsJetsCut->Fill(ib__->Pt(),evtwt_) ;
    }
  }
      } //// If at least two b-jets
    njets = 0;


    // 4. NO CUT ON NUMBER OF B-JETS ---------------------------------------- 

    double HT4(0) ;
    std::vector<TLorentzVector>higgsJets4 ;
    std::vector<TLorentzVector>bJets4 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets4.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets4.begin(); ihig != higgsJets4.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets4.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


    if (higgsJets4.size() >= 1) {

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets4.begin(); ihig != higgsJets4.end(); ++ihig) { 
          HT4 += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets4.begin(); ib != bJets4.end(); ++ib) { 
          HT4 += ib->Pt() ; 
        }

  if ( HT1 > HTMin_ && HT1 < HTMax_ ) {
    h_nHiggsJets_NoNumberBJetsCut->Fill(higgsJets4.size(),evtwt_) ;
    h_nBJets_NoNumberBJetsCut->Fill(bJets4.size(),evtwt_) ;
    h_HT_NoNumberBJetsCut->Fill(HT4,evtwt_) ;

    for (std::vector<TLorentzVector>::const_iterator ihig_ = higgsJets4.begin(); ihig_ != higgsJets4.end(); ++ihig_) {
        h_pt_HiggsJets_NoNumberBJetsCut->Fill(ihig_->Pt(),evtwt_) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib_ = bJets4.begin(); ib_ != bJets4.end(); ++ib_) {
        h_pt_BJets_NoNumberBJetsCut->Fill(ib_->Pt(),evtwt_) ;
    }
  }
    } //// If at least one Higgs jet
    njets = 0;


    // 5. NO HT CUT ---------------------------------------- 

    double HT5(0) ;
    std::vector<TLorentzVector>higgsJets5 ;
    std::vector<TLorentzVector>bJets5 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets5.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets5.begin(); ihig != higgsJets5.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

        bJets5.push_back(bjet_p4) ; 

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


    if (higgsJets5.size() >= 1) {

      if (bJets5.size() >= 2 ) { 

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets5.begin(); ihig != higgsJets5.end(); ++ihig) { 
          HT5 += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets5.begin(); ib != bJets5.end(); ++ib) { 
          HT5 += ib->Pt() ; 
        }

    h_nHiggsJets_NoHTCut->Fill(higgsJets5.size(),evtwt_) ;
    h_nBJets_NoHTCut->Fill(bJets5.size(),evtwt_) ;
    h_HT_NoHTCut->Fill(HT5,evtwt_) ;

    for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets5.begin(); ihig != higgsJets5.end(); ++ihig) {
        h_pt_HiggsJets_NoHTCut->Fill(ihig->Pt(),evtwt_) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib = bJets5.begin(); ib != bJets5.end(); ++ib) {
        h_pt_BJets_NoHTCut->Fill(ib->Pt(),evtwt_) ;
    }

      } //// If at least two b-jets
    } //// If at least one Higgs jet


  } //// entry loop 


  pt_HiggsJets_Eff = getCummulative(h_pt_HiggsJets_NoHiggsJetpTCut, "Down");
  pt_BJets_Eff = getCummulative(h_pt_BJets_NoBJetpTCut, "Down");
  nHiggsJets_Eff = getCummulative(h_nHiggsJets_NoNumberHiggsJetsCut, "Down");
  nBJets_Eff = getCummulative(h_nBJets_NoNumberBJetsCut, "Down");
  HT_Eff = getCummulative(h_HT_NoHTCut, "Down");

  pt_HiggsJets_Eff->SetName("Eff_pt_HiggsJets");
  pt_BJets_Eff->SetName("Eff_pt_BJets");
  nHiggsJets_Eff->SetName("Eff_nHiggsJets");
  nBJets_Eff->SetName("Eff_nBJets");
  HT_Eff->SetName("Eff_HT");

  pt_HiggsJets_Eff->Write();
  pt_BJets_Eff->Write();
  nHiggsJets_Eff->Write();
  nBJets_Eff->Write();
  HT_Eff->Write();

 } // doNminus1_


//==============================================================================================  


  // Analysis loop 
  // (should be in sync with the Analyzer module)
  // Also used for trigger studies 

 if ( doAnalysis_ ) {
  for(int entry=0; entry<maxEvents_; entry++) {
//  for(int entry=0; entry<1000; entry++) {

    //// Event variables 
    std::vector<TLorentzVector>fatJets ;
    std::vector<TLorentzVector>higgsJets ; 
    std::vector<TLorentzVector>jets ; 
    std::vector<TLorentzVector>bJets ; 
    std::vector<TLorentzVector>GenHiggsbb ;
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    int njets(0) ; 
    bool passHLT(false) ;
    bool   passHLTTest(false) ;
    double HT(0) ; 
    double AK5HT(0.);
    int  nGoodVtxs(0) ;

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << "B- "  << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    isData_ = EvtInfo.McFlag ? 0 : 1;
    if ( !isData_ ) evtwt_  = EvtInfo.Weight ;
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    nGoodVtxs = 0 ;
    //// Select good vertices 
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }

    //// Select events with at least one H->bb decay 
    int iHbb(0);
    for( int igen=0; igen<GenInfo.Size; igen++){
	//if( GenInfo.PdgID[igen] == 25 ) cout << endl << "  " << igen <<*/ ".) " << GenInfo.Status[igen] << ". " 
	//<< GenInfo.PdgID[igen] << ": " << GenInfo.Da0PdgID[igen] << ", " << GenInfo.Da1PdgID[igen];
        if( GenInfo.Status[igen] != 3 || GenInfo.PdgID[igen] != 25 || GenInfo.nDa[igen] < 2) continue;
        if( abs(GenInfo.Da0PdgID[igen])!=5 ) continue;
        if( GenInfo.Da0PdgID[igen]/GenInfo.Da1PdgID[igen]!=-1 ) continue;
        iHbb++;
	TLorentzVector genH_p4;
	genH_p4.SetPtEtaPhiM(GenInfo.Pt[igen], GenInfo.Eta[igen], GenInfo.Phi[igen], GenInfo.Mass[igen]);
        GenHiggsbb.push_back(genH_p4) ;
    }
    h_GenNumHiggsTobb->Fill(iHbb, evtwt_*puweight_); 
    if( iHbb==0 && genDecay_ ) continue;

    h_cutflow -> Fill("AllEvents", 1) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nFatJets"), FatJetInfo.Size, evtwt_*puweight_) ; 

    if ( TriggerStudyOn_ && MuTriggerStudyOn_ ) {
      for ( std::vector<int>::const_iterator ihlt = muHLTPaths_.begin(); ihlt != muHLTPaths_.end(); ++ihlt ) {
        if (EvtInfo.TrgBook[*ihlt] == 1) {
          passHLT = true ;
          break ;
        }
        else passHLT = false ;
      }

      if ( !passHLT ) continue ;
    }

    for ( std::vector<int>::const_iterator ihlt = hltPaths_.begin(); ihlt != hltPaths_.end(); ++ihlt ) { 
      if (EvtInfo.TrgBook[*ihlt] == 1) { 
         passHLTTest = true ; 
         break ; 
      }
      else passHLTTest = false ; 
    }

    if (!TriggerStudyOn_) {
      if ( passHLTTest ) {
        FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ;
        FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ;
        FillHisto(TString("TriggerSel")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ;
        FillHisto(TString("TriggerSel")+TString("_nFatJets"), FatJetInfo.Size, evtwt_*puweight_) ;
        h_cutflow -> Fill("TriggerSel", 1) ; 
      } else continue ; 
    }   

/* deleted by Devdatta 
    bool CA8pre(false);
    for (int ifj=0; ifj < FatJetInfo.Size; ++ifj) {
      if ( FatJetInfo.Pt[ifj] > 100. ) {
        CA8pre = true;
        break;
      }
    }
    if( !CA8pre ) continue;

    int AK5preCount(0);
    for (int ij=0; ij < JetInfo.Size; ++ij) {
      if ( JetInfo.Pt[ij] > 30. ) AK5preCount++;
    }
    if( AK5preCount < 2 ) continue;
*/

    evtwt_ *= puweight_ ; 

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
        continue; //// apply jet pT cut
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

       FillHisto(TString("TriggerSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_MassPruned")         ,FatJetInfo.MassPruned[ifatjet] ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_CombinedSVBJetTags") ,FatJetInfo.CombinedSVBJetTags[ifatjet] ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_tau2ByTau1")         ,FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau2")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau2[ifatjet] ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau1")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ;


      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dR = subjet1_p4.DeltaR(subjet2_p4);
      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations

       FillHisto(TString("TriggerSel")+TString("_SubJet1_Pt") ,subjet1_p4.Pt() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet1_Eta") ,subjet1_p4.Eta() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet1_Mass") ,subjet1_p4.Mag() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ;

       FillHisto(TString("TriggerSel")+TString("_SubJet2_Pt") ,subjet2_p4.Pt() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet2_Eta") ,subjet2_p4.Eta() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet2_Mass") ,subjet2_p4.Mag() ,evtwt_)  ;
       FillHisto(TString("TriggerSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ;


       //// Selecting fat jets
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_) {

	  fatJets.push_back(fatjet_p4) ;

         //// Higgs tagging
         if ( SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 
           FillHisto(TString("TriggerSel")+TString("_HiggsJet_Pt")   ,fatjet_p4.Pt() ,evtwt_)  ;
           FillHisto(TString("TriggerSel")+TString("_HiggsJet_Eta")  ,fatjet_p4.Eta() ,evtwt_)  ;
           FillHisto(TString("TriggerSel")+TString("_HiggsJet_Mass") ,fatjet_p4.Mag() ,evtwt_)  ;

           float DRgh(5.);
           if ( genDecay_ ) {
  	     for (std::vector<TLorentzVector>::const_iterator igenH = GenHiggsbb.begin(); igenH != GenHiggsbb.end(); ++igenH) {
	       if( fatjet_p4.DeltaR(*igenH) < DRgh )
                 DRgh = fatjet_p4.DeltaR(*igenH);
	     }
	     h_DR_GenHjet->Fill(DRgh ,evtwt_);
	   }

	   if ( !genDecay_ || ( genDecay_ && DRgh < 0.5 ) )
	        higgsJets.push_back(fatjet_p4) ;

        } //// Higgs tagging

      } //// Selecting fat jets

    } //// Loop over fat jets 


//eleni?    float AK5HT(0.);
    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;


      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);
      AK5HT += JetInfo.Pt[ijet];

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
           break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      jets.push_back(jet_p4) ; 

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.679) {

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

	bJets.push_back(bjet_p4) ;

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 

    FillHisto(TString("TriggerSel")+TString("_nJets"), njets, evtwt_) ;
    FillHisto(TString("TriggerSel")+TString("_nBJets"),bJets.size(), evtwt_) ;
    FillHisto(TString("TriggerSel")+TString("_nHJets"),higgsJets.size(), evtwt_) ;
    FillHisto(TString("TriggerSel")+TString("_AK5HT") ,AK5HT ,evtwt_)  ;

     if (fatJets.size() >= 1) {
       FillHisto(TString("FatJetSel")+TString("_nJets"), njets, evtwt_) ;
       FillHisto(TString("FatJetSel")+TString("_nBJets"), bJets.size(), evtwt_) ;
       FillHisto(TString("FatJetSel")+TString("_nHJets"), higgsJets.size(), evtwt_) ;
       FillHisto(TString("FatJetSel")+TString("_AK5HT")   ,AK5HT ,evtwt_)  ;
       for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) {
         FillHisto(TString("FatJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ;
         FillHisto(TString("FatJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ;
       }
       h_cutflow -> Fill("FatJetSel", 1) ;
     }

    if (higgsJets.size() >= 1) {

       FillHisto(TString("HiggsJetSel")+TString("_nJets"), njets, evtwt_) ;
       FillHisto(TString("HiggsJetSel")+TString("_nBJets"), bJets.size(), evtwt_) ;
       FillHisto(TString("HiggsJetSel")+TString("_nBJets"), bJets.size(), evtwt_) ;
       FillHisto(TString("HiggsJetSel")+TString("_AK5HT")   ,AK5HT ,evtwt_)  ;
       for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) {
         FillHisto(TString("HiggsJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ;
         FillHisto(TString("HiggsJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ;
       }
      h_cutflow -> Fill("HiggsJetSel", 1) ; 

      if (bJets.size() >= 2 ) { 

        float maxHiggspT(0.), maxbpT(0.), maxbEta(0.);
        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) {
          HT += ihig->Pt() ;
          if ( ihig->Pt() > maxHiggspT ) maxHiggspT = ihig->Pt();
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) {
          HT += ib->Pt() ;
          if ( ib->Pt() > maxbpT ) {
            maxbpT = ib->Pt();
            maxbEta = ib->Eta();
          }
        }

        h_cutflow -> Fill("BJetsSel", 1) ;
        FillHisto(TString("BJetsSel")+TString("_nJets"), JetInfo.Size, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_nHJets"),higgsJets.size(), evtwt_) ;
        FillHisto(TString("BJetsSel")+TString("_nBJets"),bJets.size(), evtwt_) ;
        FillHisto(TString("BJetsSel")+TString("_AK5HT") ,AK5HT ,evtwt_)  ;
        FillHisto(TString("BJetsSel")+TString("_HT")    ,HT ,evtwt_)  ;

// NOTE: The trigger study part moved to before the HT cut! 
        //// Trigger studies
        if ( TriggerStudyOn_ ) {
          h_TrigStudySel_before_HT -> Fill(HT ,evtwt_)  ;
          h_TrigStudySel_before_AK5HT -> Fill(AK5HT ,evtwt_)  ;
          h_TrigStudySel_before_HiggsJet_pT -> Fill(maxHiggspT ,evtwt_)  ;
          h_TrigStudySel_before_bJet_pT -> Fill(maxbpT ,evtwt_)  ;
          h_TrigStudySel_before_bJet_Eta -> Fill(maxbEta ,evtwt_)  ;
          h_varBins_TrigStudySel_before_HT -> Fill(HT ,evtwt_)  ;
          h_varBins_TrigStudySel_before_AK5HT -> Fill(AK5HT ,evtwt_)  ;
          h2_HT_AK5HT->Fill(HT, AK5HT, evtwt_) ;
          if (passHLTTest) {
            h_cutflow -> Fill("TriggerStudy", 1) ;
            h_TrigStudySel_HT -> Fill(HT ,evtwt_)  ;
            h_TrigStudySel_AK5HT -> Fill(AK5HT ,evtwt_)  ;
            h_TrigStudySel_HiggsJet_pT -> Fill(maxHiggspT ,evtwt_)  ;
            h_TrigStudySel_bJet_pT -> Fill(maxbpT ,evtwt_)  ;
            h_TrigStudySel_bJet_Eta -> Fill(maxbEta ,evtwt_)  ;
            h_varBins_TrigStudySel_HT -> Fill(HT ,evtwt_)  ;
            h_varBins_TrigStudySel_AK5HT -> Fill(AK5HT ,evtwt_)  ;
            h2_HT_AK5HT_passHLT->Fill(HT, AK5HT, evtwt_) ;
          }
          else {
            h_varBins_TrigStudySel_failHLT_HT -> Fill(HT ,evtwt_)  ;
            h_varBins_TrigStudySel_failHLT_AK5HT -> Fill(AK5HT ,evtwt_)  ;
            h2_HT_AK5HT_failHLT->Fill(HT, AK5HT, evtwt_) ;
          }
        }

        if (HT < HTMin_ || HT > HTMax_) continue ;

        FillHisto(TString("HTSel")+TString("_nJets"), njets, evtwt_) ;
        FillHisto(TString("HTSel")+TString("_nBJets"), bJets.size(), evtwt_) ;
        FillHisto(TString("HTSel")+TString("_nHJets"), higgsJets.size(), evtwt_) ;
        FillHisto(TString("HTSel")+TString("_HT") ,HT ,evtwt_)  ;
        FillHisto(TString("HTSel")+TString("_AK5HT") ,AK5HT ,evtwt_)  ;
        if ( !isData_ ) h_cutflow -> Fill("HTSel", 1) ;

/* Replaced by Devdatta's block above.
	//// Trigger studies 
	h_TrigStudySel_before_HT -> Fill(HT ,evtwt_)  ;
        h_TrigStudySel_before_AK5HT -> Fill(AK5HT ,evtwt_)  ;
        h_TrigStudySel_before_HiggsJet_pT -> Fill(maxHiggspT ,evtwt_)  ;
	h_TrigStudySel_before_bJet_pT -> Fill(maxbpT ,evtwt_)  ;
        h_TrigStudySel_before_bJet_Eta -> Fill(maxbEta ,evtwt_)  ;

	bool passHLT_(false);
	if ( TriggerStudyOn_ ) {
	  for ( std::vector<int>::const_iterator ihlt = hltPaths_.begin(); ihlt != hltPaths_.end(); ++ihlt ) {
	    if (EvtInfo.TrgBook[*ihlt] == 1) {
              passHLT_ = true ;
      	      break ;
	    }
	  }
	}
	if ( passHLT_ ) {
	  h_cutflow -> Fill("TriggerStudy", 1) ;
	  h_TrigStudySel_HT -> Fill(HT ,evtwt_)  ;
	  h_TrigStudySel_AK5HT -> Fill(AK5HT ,evtwt_)  ;
	  h_TrigStudySel_HiggsJet_pT -> Fill(maxHiggspT ,evtwt_)  ;
	  h_TrigStudySel_bJet_pT -> Fill(maxbpT ,evtwt_)  ;
	  h_TrigStudySel_bJet_Eta -> Fill(maxbEta ,evtwt_)  ;
	}
*/

        //// Reconstruct b' candidates

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          const TLorentzVector* closestB_p4 ;
          double deltaR(TMath::Pi()) ; 
          for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
            if ( ihig->DeltaR(*ib) < deltaR) {
              deltaR = ihig->DeltaR(*ib) ; 
              closestB_p4 = &(*ib) ; 
            }
          }
          if (deltaR < TMath::Pi()) {
            bprimes.push_back(*ihig + *closestB_p4) ; 

            FillHisto(TString("HTSel")+TString("_bprimePt") ,(*ihig + *closestB_p4).Pt() ,evtwt_)  ;
            FillHisto(TString("HTSel")+TString("_bprimeMass") ,(*ihig + *closestB_p4).Mag() ,evtwt_)  ;

          }

        } //// Reconstruct b' candidates

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  } //// entry loop 

  cout << endl << "h_TrigStudySel_before_HT: " << h_TrigStudySel_before_HT->GetEntries() << " entries, h_TrigStudySel_HT: " << h_TrigStudySel_HT->GetEntries() << endl;

  TriggerEff_HT = new TGraphAsymmErrors(h_TrigStudySel_HT, h_TrigStudySel_before_HT, "cl=0.683 b(1,1) mode") ;
  TriggerEff_HT->SetName("TriggerEff_HT") ;
  TriggerEff_HT->Write() ;

  TriggerEff_AK5HT = new TGraphAsymmErrors(h_TrigStudySel_AK5HT, h_TrigStudySel_before_AK5HT, "cl=0.683 b(1,1) mode") ;
  TriggerEff_AK5HT->SetName("TriggerEff_AK5HT") ;
  TriggerEff_AK5HT->Write() ;

  TriggerEff_HiggsJetpT = new TGraphAsymmErrors(h_TrigStudySel_HiggsJet_pT, h_TrigStudySel_before_HiggsJet_pT, "cl=0.683 b(1,1) mode") ;
  TriggerEff_HiggsJetpT->SetName("TriggerEff_HiggHiggssJetpT") ;
  TriggerEff_HiggsJetpT->Write() ; 

  TriggerEff_JetpT = new TGraphAsymmErrors(h_TrigStudySel_bJet_pT, h_TrigStudySel_before_bJet_pT, "cl=0.683 b(1,1) mode") ;
  TriggerEff_JetpT->SetName("TriggerEff_JetpT") ;
  TriggerEff_JetpT->Write() ;

  TriggerEff_JetEta = new TGraphAsymmErrors(h_TrigStudySel_bJet_Eta, h_TrigStudySel_before_bJet_Eta, "cl=0.683 b(1,1) mode") ;
  TriggerEff_JetEta->SetName("TriggerEff_JetEta") ;
  TriggerEff_JetEta->Write() ;

 } // doAnalysis_


//==============================================================================================  


  // Analysis loop for filling the mini tree 
  // Should be in sync with the Analysis module, 
  // but with the changes needed for having looser samples. 

 if ( doTree_ ) {
  for(int entry=0; entry<maxEvents_; entry++) {
//  for(int entry=0; entry<1000; entry++) {

    //// Event variables 
    std::vector<TLorentzVector>fatJets ;
    std::vector<TLorentzVector>higgsJets ; 
    std::vector<int>higgsJetsIndex ;
    std::vector<TLorentzVector>jets ; 
    std::vector<TLorentzVector>bJets ; 
    std::vector<int>bJetsIndex ;
    std::vector<TLorentzVector>GenHiggsbb ;
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    int njets(0) ; 
    bool passHLT(false) ;
    double HT(0) ; 
    int  nGoodVtxs(0) ;

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << "C- "  << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    isData_ = EvtInfo.McFlag ? 0 : 1;
    if ( !isData_ ) evtwt_  = EvtInfo.Weight ;
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    nGoodVtxs = 0 ;
    //// Select good vertices 
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }

    //// Select events with at least one H->bb decays
    if( genDecay_ ) {
        int iHbb(0);
        for( int igen=0; igen<GenInfo.Size; igen++){
            if( GenInfo.Status[igen] != 3 || GenInfo.PdgID[igen] != 25 || GenInfo.nDa[igen] < 2) continue;
            if( abs(GenInfo.Da0PdgID[igen])!=5 ) continue;
            if( GenInfo.Da0PdgID[igen]/GenInfo.Da1PdgID[igen]!=-1 ) continue;
            iHbb++;
            TLorentzVector genH_p4;
            genH_p4.SetPtEtaPhiM(GenInfo.Pt[igen], GenInfo.Eta[igen], GenInfo.Phi[igen], GenInfo.Mass[igen]);
            GenHiggsbb.push_back(genH_p4) ;
        }
        if( iHbb==0 ) continue;
    }


    for ( std::vector<int>::const_iterator ihlt = hltPaths_.begin(); ihlt != hltPaths_.end(); ++ihlt ) { 
      if (EvtInfo.TrgBook[*ihlt] == 1) { 
         passHLT = true ; 
         break ; 
      }
      else passHLT = false ; 
    }
 
    if ( !passHLT ) continue ; 

    bool CA8pre(false);
    for (int ifj=0; ifj < FatJetInfo.Size; ++ifj) {
      if ( FatJetInfo.Pt[ifj] > 100. ) {
        CA8pre = true;
        break;
      }
    }
    if( !CA8pre ) continue;

    int AK5preCount(0);
    for (int ij=0; ij < JetInfo.Size; ++ij) {
      if ( JetInfo.Pt[ij] > 30. ) AK5preCount++;
    }
    if( AK5preCount < 2 ) continue;

    evtwt_ *= puweight_ ; 


    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
        continue; //// apply jet pT cut
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.MassPruned[ifatjet] < 50. //eleni fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 
      retca8.set(false);
      if ( fatjetIDLoose(ifatjet,retca8) == 0 ) continue; //// apply loose jet ID

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);


      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dR = subjet1_p4.DeltaR(subjet2_p4);
      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      if( subjet_dyphi < (FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]) ) 
        continue; //// skip infrared unsafe configurations


       //// Selecting fat jets
      if (fatjet_p4.Mag() > 75.) { //elenifatJetMassMin_ 
//eleni          && fatjet_p4.Mag() < fatJetMassMax_ ) { 
//eleni          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_) {

	  fatJets.push_back(fatjet_p4) ;

         //// Higgs tagging
         if ( SubJetInfo.CombinedSVBJetTags[iSubJet1] >= 0. //eleni > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] >= 0. //eleni > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) {

           float DRgh(5.);
	   if ( genDecay_ ) {
             for (std::vector<TLorentzVector>::const_iterator igenH = GenHiggsbb.begin(); igenH != GenHiggsbb.end(); ++igenH) {
               if( fatjet_p4.DeltaR(*igenH) < DRgh )
                 DRgh = fatjet_p4.DeltaR(*igenH);
             }
             h_DR_GenHjet->Fill(DRgh ,evtwt_);
	   }

           if ( !genDecay_ || ( genDecay_ && DRgh < 0.5 ) ) {
             higgsJets.push_back(fatjet_p4) ;
	     higgsJetsIndex.push_back(ifatjet) ;
	   }

        } //// Higgs tagging

      } //// Selecting fat jets

    } //// Loop over fat jets 


    float AK5HT(0.);
    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      retak5.set(false);
      if ( jetIDTight(ijet,retak5) == 0 ) continue;

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);
      AK5HT += JetInfo.Pt[ijet];

      bool isJetNotHiggs(false) ; 
      for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) {
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (!isJetNotHiggs) continue ; //// Higgs-b jet disambiguation  

      jets.push_back(jet_p4) ; 

      if (JetInfo.Pt[ijet] > bjetPtMin_ && JetInfo.CombinedSVBJetTags[ijet] > 0.244 ) { //&& JetInfo.CombinedSVBJetTags[ijet] > 0.679) { //eleni

        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
            JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

	bJets.push_back(bjet_p4) ;
	bJetsIndex.push_back(ijet) ;

      } //// Select b-tagged AK5 jets 

    } //// Loop over AK5 jets 


    if (higgsJets.size() >= 1) {

      if (bJets.size() >= 2 ) { 

	float maxHiggspT(0.), maxbpT(0.), maxbEta(0.);
        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          HT += ihig->Pt() ; 
	  if ( ihig->Pt() > maxHiggspT ) maxHiggspT = ihig->Pt();
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
          HT += ib->Pt() ; 
	  if ( ib->Pt() > maxbpT ) {
	    maxbpT = ib->Pt();
	    maxbEta = ib->Eta();
	  }
        }

        //eleni if (HT < HTMin_ || HT > HTMax_) continue ; 


        //// Reconstruct b' candidates

        float maxHig(0.), Bpt(0.), maxHigBmass(0.), DR_bHiggs(0.);
        float maxHigEta(0.), maxHigPhi(0.), BEta(0.), BPhi(0.);
	float HJet_Mass(0.), HJet_PrunedMass(0.), HJet_Nsubjetiness(0.), bJet_CSV(0.);
	float Subjet_MaxPt(0.), Subjet_MinPt(0.), Subjet_MaxCSV(0.), Subjet_MinCSV(0.), Subjet_CSV_MaxPt(0.), Subjet_CSV_MinPt(0.), Subjet_DR(0.);
	int ihigi(0);
        int ibi(0);
        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          const TLorentzVector* closestB_p4 ;
          double deltaR(TMath::Pi()) ; 
          for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
            if ( ihig->DeltaR(*ib) < deltaR) {
              deltaR = ihig->DeltaR(*ib) ; 
              closestB_p4 = &(*ib) ; 
            }
          }
          if (deltaR < TMath::Pi()) {
            bprimes.push_back(*ihig + *closestB_p4) ; 

            // For tree filling ---- 
            if ( ihig->Pt() > maxHig ) {

              maxHig = ihig->Pt();
              maxHigEta = ihig->Eta();
              maxHigPhi = ihig->Phi();
              Bpt = closestB_p4->Pt();
              BEta = closestB_p4->Eta();
              BPhi = closestB_p4->Phi();
	      DR_bHiggs = deltaR;
	      maxHigBmass = (*ihig + *closestB_p4).M();

	      int ihigindex = higgsJetsIndex[ihigi];
	      HJet_Mass = FatJetInfo.Mass[ihigindex];
	      HJet_PrunedMass = FatJetInfo.MassPruned[ihigindex];
	      HJet_Nsubjetiness = FatJetInfo.tau2[ihigindex]/FatJetInfo.tau1[ihigindex];
	      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ihigindex];
              int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ihigindex];
	      int pt1 = SubJetInfo.Pt[iSubJet1];
              int pt2 = SubJetInfo.Pt[iSubJet2];
              float csv1 = SubJetInfo.CombinedSVBJetTags[iSubJet1];
              float csv2 = SubJetInfo.CombinedSVBJetTags[iSubJet2];
              Subjet_MaxPt = ( pt1 > pt2 ? pt1 : pt2 );
              Subjet_MinPt = ( pt1 > pt2 ? pt2 : pt1 );
              Subjet_MaxCSV = ( csv1 > csv2 ? csv1 : csv2 );
              Subjet_MinCSV = ( csv1 > csv2 ? csv2 : csv1 );
              Subjet_CSV_MaxPt = ( pt1 > pt2 ? csv1 : csv2 );
              Subjet_CSV_MinPt = ( pt1 > pt2 ? csv2 : csv1 );
              //Subjet_DR = 
              int ibindex = bJetsIndex[ibi];
	      bJet_CSV = JetInfo.CombinedSVBJetTags[ibi];

            } // -------------------

          }
	  ++ihigi;
        } //// Reconstruct b' candidates

        //// Fill the tree
        final_nHiggsJets = higgsJets.size();
        final_nBJets = bJets.size();
        final_HiggsJets_pt = maxHig;
        final_BJets_pt = Bpt;
        final_HiggsJets_eta = maxHigEta;
        final_BJets_eta = BEta;
        final_HT = HT;
        final_evtwt = evtwt_;
	final_AK5HT = AK5HT;
	final_HiggsJets_phi = maxHigPhi;
	final_BJets_phi = BPhi;
	final_DR_bHiggs = DR_bHiggs;
	final_Bprime_mass = maxHigBmass;
	final_HJet_Mass = HJet_Mass;
	final_HJet_PrunedMass = HJet_PrunedMass;
	final_HJet_Nsubjetiness = HJet_Nsubjetiness;
	final_Subjet_MaxPt = Subjet_MaxPt;
	final_Subjet_MinPt = Subjet_MinPt;
	final_Subjet_CSV_MaxPt = Subjet_CSV_MaxPt;
	final_Subjet_CSV_MinPt = Subjet_CSV_MinPt;
	final_Subjet_MaxCSV = Subjet_MaxCSV;
	final_Subjet_MinCSV = Subjet_MinCSV;
	final_Subjet_DR = Subjet_DR;
	final_bJet_CSV = bJet_CSV;

        // NOTE: Events with maxHig=0 and Bpt=0 have not passed the B' reconstruction.
        // NOTE: Use the final_evtwt to weight events in the TMVA macro.

        tree->Fill();

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  } //// entry loop 

 } // doTree_

}


// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobHDevelopment::endJob() { 

  tree->Write();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimeTobHDevelopment::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobHDevelopment);
