#include <TH1D.h>
#include <TH1I.h>
#include <TEfficiency.h> 

class Histograms {

  public: 

    Histograms () ; 
    ~Histograms () {} ; 

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

};

Histograms::Histograms () {
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

  teff_HiggsJetMatch           = fs->make<TEfficiency>("teff_HiggsJetMatch" ,"Higgs-Jet matching efficiency" ,100 ,0.  ,2000.) ; 

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

}

