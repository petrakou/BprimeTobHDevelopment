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
#include <TH1D.h>
#include <TH1I.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Bprime_kit/BprimeTobH/interface/format.h"
#include "Bprime_kit/BprimeTobH/interface/TriggerBooking.h"
#include "Bprime_kit/BprimeTobH/interface/Njettiness.hh"
#include "Bprime_kit/BprimeTobH/interface/Nsubjettiness.hh"

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

    virtual void CreateHistos(const TString&) ; 
    virtual void AddHisto(const TString&, const TString& ,const int&, const double&, const double&) ; 
    template <class Type>
      void FillHisto(const TString& name, const Type value, const double weight);

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;

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
    const double fatJetTau2ByTau1Min_ ; 
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
//eleni
    bool TriggerStudyOn_ ;

    TH1I*              h_FatJets_Index;
    TH1I*              h_FatJets_NTracks;
    TH1D*              h_FatJets_Et;
    TH1D*              h_FatJets_Pt;
    TH1D*              h_FatJets_Eta;
    TH1D*              h_FatJets_Phi;
    TH1D*              h_FatJets_Energy;
    TH1D*              h_FatJets_Px;
    TH1D*              h_FatJets_Py;
    TH1D*              h_FatJets_Pz;
    TH1D*              h_FatJets_Mass;
    TH1D*              h_FatJets_Area;
    TH1I*              h_FatJets_JetCharge;
    TH1I*              h_FatJets_NConstituents;
    TH1D*              h_FatJets_NCH;
    TH1D*              h_FatJets_CEF;
    TH1D*              h_FatJets_NHF;
    TH1D*              h_FatJets_NEF;
    TH1D*              h_FatJets_CHF;
    TH1D*              h_FatJets_QGTagsMLP;
    TH1D*              h_FatJets_QGTagsLikelihood;
    TH1I*              h_FatJets_JetIDLOOSE;
    TH1I*              h_FatJets_JetIDTIGHT;
    TH1D*              h_FatJets_PtCorrRaw;
    TH1D*              h_FatJets_PtCorrL2;
    TH1D*              h_FatJets_PtCorrL3;
    TH1D*              h_FatJets_PtCorrL7g;
    TH1D*              h_FatJets_PtCorrL7uds;
    TH1D*              h_FatJets_PtCorrL7c;
    TH1D*              h_FatJets_PtCorrL7b;
    TH1D*              h_FatJets_JetBProbBJetTags;
    TH1D*              h_FatJets_JetProbBJetTags; 
    TH1D*              h_FatJets_TrackCountHiPurBJetTags; 
    TH1D*              h_FatJets_CombinedSVBJetTags; 
    TH1D*              h_FatJets_CombinedSVMVABJetTags; 
    TH1D*              h_FatJets_SoftElecByIP3dBJetTags; 
    TH1D*              h_FatJets_SoftElecByPtBJetTags;
    TH1D*              h_FatJets_SoftMuonBJetTags;
    TH1D*              h_FatJets_SoftMuonByIP3dBJetTags; 
    TH1D*              h_FatJets_SoftMuonByPtBJetTags;
    TH1D*              h_FatJets_DoubleSVHighEffBJetTags; 
    TH1D*              h_FatJets_GenJetPt; 
    TH1D*              h_FatJets_GenJetEta;
    TH1D*              h_FatJets_GenJetPhi;
    TH1I*              h_FatJets_GenMCTag;
    TH1D*              h_FatJets_GenPt; // for partons 
    TH1D*              h_FatJets_GenEta;
    TH1D*              h_FatJets_GenPhi;
    TH1D*              h_FatJets_GenPdgID;
    TH1D*              h_FatJets_GenFlavor;
    // Specific to Fat Jets: 
    TH1D*              h_FatJets_EtPruned;
    TH1D*              h_FatJets_PtPruned;
    TH1D*              h_FatJets_EtaPruned;
    TH1D*              h_FatJets_PhiPruned;
    TH1D*              h_FatJets_EnergyPruned;
    TH1D*              h_FatJets_PxPruned; 
    TH1D*              h_FatJets_PyPruned; 
    TH1D*              h_FatJets_PzPruned; 
    TH1D*              h_FatJets_MassPruned; 
    TH1D*              h_FatJets_AreaPruned; 
    TH1I*              h_FatJets_Jet_SubJet1Idx;
    TH1I*              h_FatJets_Jet_SubJet2Idx;
    TH1D*              h_FatJets_tau1;
    TH1D*              h_FatJets_tau2;
    TH1D*              h_FatJets_tau3;
    TH1D*              h_FatJets_tau2ByTau1; 
    TH1D*              h_FatJets_tau3ByTau2; 
    TH1D*              h_FatJets_tau3ByTau1; 

    TH1D*              h_SubJets_Pt;
    // Specific to Subjets:
    TH1I*              h_SubJets_Jet_FatJetIdx;

    TH1D* h_SubJet1_Pt ; 
    TH1D* h_SubJet1_Eta ;
    TH1D* h_SubJet1_Mass ; 
    TH1D* h_SubJet1_CombinedSVBJetTags ; 

    TH1D* h_SubJet2_Pt ; 
    TH1D* h_SubJet2_Eta ;
    TH1D* h_SubJet2_Mass ; 
    TH1D* h_SubJet2_CombinedSVBJetTags ; 

    TH1D* h_HiggsJet_Pt ; 
    TH1D* h_HiggsJet_Eta ;
    TH1D* h_HiggsJet_Mass ; 

    TH1D* h_nJets ; 
    TH1D* h_nBJets ; 
    TH1D* h_nHJets ; 

    TH1D* h_HT ; 
    TH1D* h_bprimePt ; 
    TH1D* h_bprimeMass ; 

    TH1D* h_HiggsPt ; 
    TH1D* h_HiggsPtMatchedJet ; 
    TEfficiency* teff_HiggsJetMatch ; 

    TH1D* h_cutflow ; 

    std::map<TString, TH1D*> hmap_1d ;  

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

/*    TH1D* h_nHiggsJets_NoCuts ;
    TH1D* h_nBJets_NoCuts ;
    TH1D* h_HT_NoCuts ;
    TH1D* h_pt_HiggsJets_NoCuts ;
    TH1D* h_pt_BJets_NoCuts ;
*/
};

//
// constructors and destructor
//
BprimeTobHDevelopment::BprimeTobHDevelopment(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
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
  fatJetTau2ByTau1Min_(iConfig.getParameter<double>("FatJetTau2ByTau1Min")),
  subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
  subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
  subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
  subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
  HTMin_(iConfig.getParameter<double>("HTMin")), 
  HTMax_(iConfig.getParameter<double>("HTMax")),
  isData_(0),
  evtwt_(1),
  TriggerStudyOn_(0), //eleni
  jetPtPreselection_(30.)
{ 

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

  // example from ExcitedQuark 
  //eSelector = new eventSelector(SelectionParameters,EvtInfo,LepInfo,JetInfo,VtxInfo);
  //cutLevels = eSelector->getCutLevels();

  h_FatJets_Pt                 = fs->make<TH1D>("h_FatJets_Pt"                 ,"Fat jet p_{T} [GeV]"       ,100 ,0.  ,2000.);
  h_FatJets_Eta                = fs->make<TH1D>("h_FatJets_Eta"                ,"Fat jet #eta"              ,50  ,-4. ,4.   );
  h_FatJets_Mass               = fs->make<TH1D>("h_FatJets_Mass"               ,"Fat jet mass [GeV]"        ,100 ,0.  ,2000.); 
  h_FatJets_MassPruned         = fs->make<TH1D>("h_FatJets_MassPruned"         ,"Fat jet pruned mass [GeV]" ,100 ,0.  ,2000.); 
  h_FatJets_tau2ByTau1         = fs->make<TH1D>("h_FatJets_tau2ByTau1"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
  h_FatJets_tau3ByTau2         = fs->make<TH1D>("h_FatJets_tau3ByTau2"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
  h_FatJets_tau3ByTau1         = fs->make<TH1D>("h_FatJets_tau3ByTau1"         ,"Fat jet #tau2/#tau1"       ,20  ,0.  ,1.   ); 
  h_FatJets_CombinedSVBJetTags = fs->make<TH1D>("h_FatJets_CombinedSVBJetTags" ,"Fat jet CSV discriminator" ,20  ,0.  ,1.   ); 

  h_SubJet1_Pt                 = fs->make<TH1D>("h_SubJet1_Pt"                 ,"SubJet1 p_{T} [GeV]"       ,100 ,0.  ,2000.);
  h_SubJet1_Eta                = fs->make<TH1D>("h_SubJet1_Eta"                ,"SubJet1 #eta"              ,50  ,-4. ,4.   );
  h_SubJet1_Mass               = fs->make<TH1D>("h_SubJet1_Mass"               ,"SubJet1 mass [GeV]"        ,100 ,0.  ,2000.); 
  h_SubJet1_CombinedSVBJetTags = fs->make<TH1D>("h_SubJet1_CombinedSVBJetTags" ,"SubJet1 CSV discriminator" ,20  ,0.  ,1.   ); 

  h_SubJet2_Pt                 = fs->make<TH1D>("h_SubJet2_Pt"                 ,"SubJet2 p_{T} [GeV]"       ,100 ,0.  ,2000.);
  h_SubJet2_Eta                = fs->make<TH1D>("h_SubJet2_Eta"                ,"SubJet2 #eta"              ,50  ,-4. ,4.   );
  h_SubJet2_Mass               = fs->make<TH1D>("h_SubJet2_Mass"               ,"SubJet2 mass [GeV]"        ,100 ,0.  ,2000.); 
  h_SubJet2_CombinedSVBJetTags = fs->make<TH1D>("h_SubJet2_CombinedSVBJetTags" ,"SubJet2 CSV discriminator" ,20  ,0.  ,1.   ); 

  h_HiggsJet_Pt                = fs->make<TH1D>("h_HiggsJet_Pt"                ,"Higgs jet p_{T} [GeV]"     ,100 ,0.  ,2000.);
  h_HiggsJet_Eta               = fs->make<TH1D>("h_HiggsJet_Eta"               ,"Higgs jet #eta"            ,50  ,-4. ,4.   );
  h_HiggsJet_Mass              = fs->make<TH1D>("h_HiggsJet_Mass"              ,"Higgs jet mass [GeV]"      ,100 ,0.  ,200. ); 

  h_HiggsPt                    = fs->make<TH1D>("h_HiggsPt"                   ,"Higgs p_{T} [GeV]"          ,100 ,0.  ,2000.);
  h_HiggsPtMatchedJet          = fs->make<TH1D>("h_HiggsPtMatchedJet"         ,"Matched Higgs p_{T} [GeV]"  ,100 ,0.  ,2000.);

  h_nJets                      = fs->make<TH1D>("h_nJets"                     ,"N_{jets}"                   ,51  ,-0.5,50.5 ); 
  h_nBJets                     = fs->make<TH1D>("h_nBJets"                    ,"N_{b jets}"                 ,51  ,-0.5,50.5 ); 
  h_nHJets                     = fs->make<TH1D>("h_nHJets"                    ,"N_{Higgs jets}"             ,51  ,-0.5,50.5 ); 

  h_HT                         = fs->make<TH1D>("h_HT"                        ,"H_{T}[GeV]"                 ,200 ,0.  ,4000.);
  h_bprimePt                   = fs->make<TH1D>("h_bprimePt"                  ,"b' p_{T} [GeV]"             ,100 ,0.  ,2000.);
  h_bprimeMass                 = fs->make<TH1D>("h_bprimeMass"                ,"b' mass [GeV]"              ,40  ,0.  ,2000.);

  h_cutflow                    = fs->make<TH1D>("h_cutflow"                   ,"Cut flow"                   ,20  ,0.  ,20.  ); 

  teff_HiggsJetMatch           = fs->make<TEfficiency>("teff_HiggsJetMatch" ,"Higgs-Jet matching efficiency" ,100 ,0.  ,2000.) ; 


  h_nHiggsJets_NoHiggsJetpTCut 	= fs->make<TH1D>("h_nHJets_NoHiggsJetpTCut" 	    ,"N_{Higgs jets}"             ,11  ,-0.5,10.5 );
  h_nBJets_NoHiggsJetpTCut 	= fs->make<TH1D>("h_nBJets_NoHiggsJetpTCut"         ,"N_{b jets}"		  ,21  ,-0.5,20.5 );
  h_HT_NoHiggsJetpTCut 		= fs->make<TH1D>("h_HT_NoHiggsJetpTCut"             ,"H_{T}[GeV]"		  ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoHiggsJetpTCut = fs->make<TH1D>("h_HiggsJet_Pt_NoHiggsJetpTCut"   ,"Higgs jet p_{T} [GeV]"	  ,100 ,0.  ,2000.);
  h_pt_BJets_NoHiggsJetpTCut 	= fs->make<TH1D>("h_BJet_Pt_NoHiggsJetpTCut"        ,"b jet p_{T} [GeV]"     	  ,100 ,0.  ,2000.);

  h_nHiggsJets_NoBJetpTCut 	= fs->make<TH1D>("h_nHJets_NoBJetpTCut"             ,"N_{Higgs jets}"		  ,11  ,-0.5,10.5 );
  h_nBJets_NoBJetpTCut 		= fs->make<TH1D>("h_nBJets_NoBJetpTCut"             ,"N_{b jets}"		  ,21  ,-0.5,20.5 );
  h_HT_NoBJetpTCut 		= fs->make<TH1D>("h_HT_NoBJetpTCut"                 ,"H_{T}[GeV]"		  ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoBJetpTCut 	= fs->make<TH1D>("h_HiggsJet_Pt_NoBJetpTCut"        ,"Higgs jet p_{T} [GeV]"      ,100 ,0.  ,2000.);
  h_pt_BJets_NoBJetpTCut 	= fs->make<TH1D>("h_BJet_Pt_NoBJetpTCut"            ,"b jet p_{T} [GeV]"     	  ,100 ,0.  ,2000.);

  h_nHiggsJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_nHJets_NoNumberHiggsJetsCut"        ,"N_{Higgs jets}"		  ,11  ,-0.5,10.5 );
  h_nBJets_NoNumberHiggsJetsCut	= fs->make<TH1D>("h_nBJets_NoNumberHiggsJetsCut"            ,"N_{b jets}"		  ,21  ,-0.5,20.5 );
  h_HT_NoNumberHiggsJetsCut 	= fs->make<TH1D>("h_HT_NoNumberHiggsJetsCut"                ,"H_{T}[GeV]"		  ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_HiggsJet_Pt_NoNumberHiggsJetsCut" ,"Higgs jet p_{T} [GeV]"	  ,100 ,0.  ,2000.);
  h_pt_BJets_NoNumberHiggsJetsCut = fs->make<TH1D>("h_BJet_Pt_NoNumberHiggsJetsCut"         ,"b jet p_{T} [GeV]"     	  ,100 ,0.  ,2000.);

  h_nHiggsJets_NoNumberBJetsCut	= fs->make<TH1D>("h_nHJets_NoNumberBJetsCut"                ,"N_{Higgs jets}"		  ,11  ,-0.5,10.5 );
  h_nBJets_NoNumberBJetsCut 	= fs->make<TH1D>("h_nBJets_NoNumberBJetsCut"                ,"N_{b jets}"		  ,21  ,-0.5,20.5 );
  h_HT_NoNumberBJetsCut 	= fs->make<TH1D>("h_HT_NoNumberBJetsCut"                    ,"H_{T}[GeV]"		  ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoNumberBJetsCut = fs->make<TH1D>("h_HiggsJet_Pt_NoNumberBJetsCut"         ,"Higgs jet p_{T} [GeV]"	  ,100 ,0.  ,2000.);
  h_pt_BJets_NoNumberBJetsCut 	= fs->make<TH1D>("h_BJet_Pt_NoNumberBJetsCut"               ,"b jet p_{T} [GeV]"     	  ,100 ,0.  ,2000.);

  h_nHiggsJets_NoHTCut 		= fs->make<TH1D>("h_nHJets_NoHTCut"                  ,"N_{Higgs jets}"  	   ,11  ,-0.5,10.5 );
  h_nBJets_NoHTCut 		= fs->make<TH1D>("h_nBJets_NoHTCut"                  ,"N_{b jets}"		   ,21  ,-0.5,20.5 );
  h_HT_NoHTCut 			= fs->make<TH1D>("h_HT_NoHTCut"                      ,"H_{T}[GeV]"		   ,200 ,0.  ,4000.);
  h_pt_HiggsJets_NoHTCut 	= fs->make<TH1D>("h_HiggsJet_Pt_NoHTCut"             ,"Higgs jet p_{T} [GeV]"	   ,100 ,0.  ,2000.);
  h_pt_BJets_NoHTCut 		= fs->make<TH1D>("h_BJet_Pt_NoHTCut"                 ,"b jet p_{T} [GeV]"     	   ,100 ,0.  ,2000.);


  h_FatJets_Pt                 -> Sumw2() ; 
  h_FatJets_Eta                -> Sumw2() ; 
  h_FatJets_Mass               -> Sumw2() ; 
  h_FatJets_MassPruned         -> Sumw2() ; 
  h_FatJets_tau2ByTau1         -> Sumw2() ; 
  h_FatJets_tau3ByTau2         -> Sumw2() ; 
  h_FatJets_tau3ByTau1         -> Sumw2() ; 
  h_FatJets_CombinedSVBJetTags -> Sumw2() ; 

  h_SubJet1_Pt                 -> Sumw2() ; 
  h_SubJet1_Eta                -> Sumw2() ; 
  h_SubJet1_Mass               -> Sumw2() ; 
  h_SubJet1_CombinedSVBJetTags -> Sumw2() ; 

  h_SubJet2_Pt                 -> Sumw2() ; 
  h_SubJet2_Eta                -> Sumw2() ; 
  h_SubJet2_Mass               -> Sumw2() ; 
  h_SubJet2_CombinedSVBJetTags -> Sumw2() ; 

  h_HiggsJet_Pt                -> Sumw2() ;
  h_HiggsJet_Eta               -> Sumw2() ;
  h_HiggsJet_Mass              -> Sumw2() ;

  h_HiggsPt                    -> Sumw2() ; 
  h_HiggsPtMatchedJet          -> Sumw2() ; 

  h_nJets  -> Sumw2() ; 
  h_nBJets -> Sumw2() ; 
  h_nHJets -> Sumw2() ; 

  h_HT         -> Sumw2() ; 
  h_bprimePt   -> Sumw2() ; 
  h_bprimeMass -> Sumw2() ; 

  h_cutflow    -> Sumw2() ; 

  h_cutflow -> GetXaxis() -> SetBinLabel(1,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2,"TriggerSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(3,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5,"HTSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(6,"TriggerStudy") ;

  for (int ii = 1; ii <= 6; ++ii) 
    CreateHistos(h_cutflow->GetXaxis()->GetBinLabel(ii)) ; 

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

  return ;  

}

void BprimeTobHDevelopment::CreateHistos(const TString& cutname) {

  AddHisto(cutname ,"_nJets"  ,20 ,-0.5 ,19.5) ; 
  AddHisto(cutname ,"_nBJets" ,20 ,-0.5 ,19.5) ; 
  AddHisto(cutname ,"_nHJets" ,20 ,-0.5 ,19.5) ; 

  return ; 
}

void BprimeTobHDevelopment::AddHisto(const TString& cutname, const TString& histname, const int& nbins, const double& min, const double& max) { 

  TH1D* h1d_mc ; 
  h1d_mc = fs->make<TH1D>(cutname+histname+"_mc", cutname+histname+"_mc", nbins, max, min);  
  h1d_mc -> Sumw2() ; 
  hmap_1d[cutname+histname+"_mc"] = h1d_mc ; 

  TH1D* h1d_data ; 
  h1d_data = fs->make<TH1D>(cutname+histname+"_data", cutname+histname+"_data", nbins, max, min);  
  h1d_data -> Sumw2() ; 
  hmap_1d[cutname+histname+"_data"] = h1d_data ; 

  return ; 

}

template <class Type>
void BprimeTobHDevelopment::FillHisto(const TString& name, const Type value, const double weight = 1.){
  if (!isData_) hmap_1d[name+"_mc"]->Fill(double(value),weight);
  else hmap_1d[name+"_data"]->Fill(double(value),weight); 

  return ; 

}


// ------------ method called for each event  ------------
void BprimeTobHDevelopment::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {

    //// Event variables 
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    int njets(0) ; 
    double HT(0) ; 
    int  nGoodVtxs(0) ;

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    nGoodVtxs = 0 ;
    //// Select good vertices 
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { cout << endl << "Vtx!" << endl; continue ; }
/*
    isData_ = EvtInfo.McFlag ? 0 : 1; 
//eleni    evtwt_  = GenInfo.Weight ; 

    h_cutflow -> Fill("AllEvents", 1) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_) ; 
*/
    if (EvtInfo.TrgBook[3225]==1||EvtInfo.TrgBook[4893]==1) {
//      h_cutflow -> Fill("TriggerSel", 1) ;
    } else {
      if (!TriggerStudyOn_)
	continue;
    }


    // 1. NO HIGGS JETS PT CUT ---------------------------------------- 

    double HT1(0) ;
    std::vector<TLorentzVector>higgsJets1 ;
    std::vector<TLorentzVector>bJets1 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) { }
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets1.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

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

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  if ( !(HT1 < HTMin_ || HT1 > HTMax_) && higgsJets1.size()!=0 && bJets1.size()!=0) {  
    h_nHiggsJets_NoHiggsJetpTCut->Fill(higgsJets1.size()) ;
    h_nBJets_NoHiggsJetpTCut->Fill(bJets1.size()) ;
    h_HT_NoHiggsJetpTCut->Fill(HT1) ;

    for (std::vector<TLorentzVector>::const_iterator ihig_ = higgsJets1.begin(); ihig_ != higgsJets1.end(); ++ihig_) {
        h_pt_HiggsJets_NoHiggsJetpTCut->Fill(ihig_->Pt()) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib_ = bJets1.begin(); ib_ != bJets1.end(); ++ib_) {
        h_pt_BJets_NoHiggsJetpTCut->Fill(ib_->Pt()) ;
    }
  }
    njets = 0;


    // 2. NO B-JETS PT CUT ---------------------------------------- 

    double HT2(0) ;
    std::vector<TLorentzVector>higgsJets2 ;
    std::vector<TLorentzVector>bJets2 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ )
	continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets2.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) {}
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

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

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  if ( !(HT1 < HTMin_ || HT1 > HTMax_) && higgsJets1.size()!=0 && bJets1.size()!=0) {
    h_nHiggsJets_NoBJetpTCut->Fill(higgsJets2.size()) ;
    h_nBJets_NoBJetpTCut->Fill(bJets2.size()) ;
    h_HT_NoBJetpTCut->Fill(HT2) ;

    for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets2.begin(); ihig != higgsJets2.end(); ++ihig) {
        h_pt_HiggsJets_NoBJetpTCut->Fill(ihig->Pt()) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib = bJets2.begin(); ib != bJets2.end(); ++ib) {
        h_pt_BJets_NoBJetpTCut->Fill(ib->Pt()) ;
    }
  }
    njets = 0;


    // 3. NO CUT ON NUMBER OF HIGGS JETS ---------------------------------------- 

    double HT3(0) ;
    std::vector<TLorentzVector>higgsJets3 ;
    std::vector<TLorentzVector>bJets3 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
	continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 

        higgsJets3.push_back(fatjet_p4) ; 

      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtPreselection_ ) continue;
      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

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

      } //// If at least two b-jets 

  if ( !(HT1 < HTMin_ || HT1 > HTMax_) && higgsJets1.size()!=0 && bJets1.size()!=0) {
    h_nHiggsJets_NoNumberHiggsJetsCut->Fill(higgsJets3.size()) ;
    h_nBJets_NoNumberHiggsJetsCut->Fill(bJets3.size()) ;
    h_HT_NoNumberHiggsJetsCut->Fill(HT3) ;

    for (std::vector<TLorentzVector>::const_iterator ihig__ = higgsJets3.begin(); ihig__ != higgsJets3.end(); ++ihig__) {
        h_pt_HiggsJets_NoNumberHiggsJetsCut->Fill(ihig__->Pt()) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib__ = bJets3.begin(); ib__ != bJets3.end(); ++ib__) {
        h_pt_BJets_NoNumberHiggsJetsCut->Fill(ib__->Pt()) ;
    }
  }
    njets = 0;


    // 4. NO CUT ON NUMBER OF B-JETS ---------------------------------------- 

    double HT4(0) ;
    std::vector<TLorentzVector>higgsJets4 ;
    std::vector<TLorentzVector>bJets4 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
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
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

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

    } //// If at least one Higgs jet 

  if ( !(HT1 < HTMin_ || HT1 > HTMax_) && higgsJets1.size()!=0 && bJets1.size()!=0) {
    h_nHiggsJets_NoNumberBJetsCut->Fill(higgsJets4.size()) ;
    h_nBJets_NoNumberBJetsCut->Fill(bJets4.size()) ;
    h_HT_NoNumberBJetsCut->Fill(HT4) ;

    for (std::vector<TLorentzVector>::const_iterator ihig_ = higgsJets4.begin(); ihig_ != higgsJets4.end(); ++ihig_) {
        h_pt_HiggsJets_NoNumberBJetsCut->Fill(ihig_->Pt()) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib_ = bJets4.begin(); ib_ != bJets4.end(); ++ib_) {
        h_pt_BJets_NoNumberBJetsCut->Fill(ib_->Pt()) ;
    }
  }
    njets = 0;


    // 5. NO HT CUT ---------------------------------------- 

    double HT5(0) ;
    std::vector<TLorentzVector>higgsJets5 ;
    std::vector<TLorentzVector>bJets5 ;

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < jetPtPreselection_ ) continue;
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue;
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
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
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

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

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 


  if ( HT1 > HTMin_ && higgsJets1.size()!=0 && bJets1.size()!=0) {
    h_nHiggsJets_NoHTCut->Fill(higgsJets5.size()) ;
    h_nBJets_NoHTCut->Fill(bJets5.size()) ;
    h_HT_NoHTCut->Fill(HT5) ;

    for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets5.begin(); ihig != higgsJets5.end(); ++ihig) {
        h_pt_HiggsJets_NoHTCut->Fill(ihig->Pt()) ;
    }
    for (std::vector<TLorentzVector>::const_iterator ib = bJets5.begin(); ib != bJets5.end(); ++ib) {
        h_pt_BJets_NoHTCut->Fill(ib->Pt()) ;
    }
  }

  } //// entry loop 


//==============================================================================================  

  for(int entry=0; entry<maxEvents_; entry++) {

    //// Event variables 
    std::vector<TLorentzVector>higgsJets ; 
    std::vector<TLorentzVector>jets ; 
    std::vector<TLorentzVector>bJets ; 
    //std::vector<std::pair<int,TLorentzVector> > higgsJets ; 
    //std::vector<std::pair<int,TLorentzVector> > jets ; 
    //std::vector<std::pair<int,TLorentzVector> > bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    int njets(0) ; 
    double HT(0) ; 
    int  nGoodVtxs(0) ;

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    nGoodVtxs = 0 ;
    //// Select good vertices 
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) {
      if (   VtxInfo.Type[iVtx]==1
          && VtxInfo.isFake[iVtx]==false
          && VtxInfo.Ndof[iVtx]>4
          && VtxInfo.Rho[iVtx]<2.
          && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }
    }
    if (nGoodVtxs < 1)  { cout << endl << "Vtx!" << endl; continue ; }

    isData_ = EvtInfo.McFlag ? 0 : 1; 
    evtwt_  = EvtInfo.Weight ; 

    h_cutflow -> Fill("AllEvents", 1) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_) ; 

    if (EvtInfo.TrgBook[3225]==1||EvtInfo.TrgBook[4893]==1) {
      h_cutflow -> Fill("TriggerSel", 1) ;
    } else {
      if (!TriggerStudyOn_)
	continue;
    }


    TLorentzVector higgs_p4 ; 
    for (int igen=0; igen < GenInfo.Size; ++igen) {

      if ( GenInfo.Status[igen] == 3 
          && TMath::Abs(GenInfo.PdgID[igen])==25 
          //&& GenInfo.nDa[igen] == 2 
          //&& TMath::Abs(GenInfo.PdgID[GenInfo.Da1[igen]])==5 
          //&& TMath::Abs(GenInfo.PdgID[GenInfo.Da2[igen]])==5 
         ) { 

        higgs_p4.SetPtEtaPhiM(GenInfo.Pt[igen], GenInfo.Eta[igen], GenInfo.Phi[igen], GenInfo.Mass[igen]) ; 
        TLorentzVector fatjet_p4;
        bool matched ; 

        h_HiggsPt ->Fill(higgs_p4.Pt()) ; 

        for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) { 

//eleni 1.5??? 
          if ( fabs(FatJetInfo.Eta[ifatjet]) > 1.5 ) continue ; 

          fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
              FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

          if (higgs_p4.DeltaR(fatjet_p4) < 0.5) {
            matched = true ; 
            h_HiggsPtMatchedJet ->Fill(higgs_p4.Pt()) ; 
            break ; 
          }
          else 
            matched = false ; 

        } //// Loop over fat jets 

        teff_HiggsJetMatch -> Fill(matched, GenInfo.Pt[igen]) ; 

      } //// Get Higgs boson 

    } //// Loop over all gen particles 


    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {

      //Fix h_FatJets_Pt->Fill(FatJetInfo.Pt[ifatjet]);

      //// Fat jet selection
      if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) 
        continue; //// apply jet pT cut
      if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) 
        continue; //// apply jet eta cut
      if ( FatJetInfo.JetIDLOOSE[ifatjet]==0 ) 
        continue; //// apply loose jet ID
      if ( FatJetInfo.MassPruned[ifatjet] < fatJetPrunedMassMin_ 
          || FatJetInfo.MassPruned[ifatjet] > fatJetPrunedMassMax_ ) 
        continue; //// apply pruned jet mass cut 

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

      //// Higgs tagging
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Min_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
          && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 
	higgsJets.push_back(fatjet_p4) ;
      }

    } //// Loop over fat jets 


    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) { 

      if ( JetInfo.Pt[ijet] < jetPtMin_ || JetInfo.Pt[ijet] > jetPtMax_ ) continue ; 
      if ( fabs(JetInfo.Eta[ijet]) > jetAbsEtaMax_ ) continue ; 
      if ( JetInfo.JetIDTIGHT[ijet]==0 ) continue ; 

      ++njets ; 

      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], 
          JetInfo.Phi[ijet], JetInfo.Mass[ijet]);

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


    if (higgsJets.size() >= 1) {

      h_cutflow -> Fill("HiggsJetSel", 1) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nJets"), njets, 1) ; 

      if (bJets.size() >= 2 ) { 

        h_cutflow -> Fill("BJetsSel", 1) ;
        FillHisto(TString("BJetsSel")+TString("_nJets"), JetInfo.Size, evtwt_) ; 

        for (std::vector<TLorentzVector>::const_iterator ihig = higgsJets.begin(); ihig != higgsJets.end(); ++ihig) { 
          HT += ihig->Pt() ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ib = bJets.begin(); ib != bJets.end(); ++ib) { 
          HT += ib->Pt() ; 
        }

        if (HT < HTMin_ || HT > HTMax_) continue ; 

        h_cutflow -> Fill("HTSel", 1) ; 
        FillHisto(TString("HTSel")+TString("_nJets"), njets, 1) ; 
        FillHisto(TString("HTSel")+TString("_nBJets"), bJets.size(), 1) ; 
        FillHisto(TString("HTSel")+TString("_nHJets"), higgsJets.size(), 1) ; 

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
            h_bprimePt   -> Fill((*ihig + *closestB_p4).Pt()) ; 
            h_bprimeMass-> Fill((*ihig + *closestB_p4).Mag()) ;
          }
        }

        h_nJets  -> Fill(njets) ; 
        h_nBJets -> Fill(bJets.size()) ; 
        h_nHJets -> Fill(higgsJets.size()) ; 

        h_HT -> Fill(HT) ; 

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  } //// entry loop 
//*/
//==============================================================================================  
  
}


// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobHDevelopment::endJob() { 
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
