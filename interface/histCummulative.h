#include <TH1F.h>
#include <TH1D.h>
#include <TString.h> 
#include <TGraphAsymmErrors.h>
#include <TFile.h>

#include <cmath> 
#include <iostream>

using namespace std ; 

template <class Type> 
TGraphAsymmErrors* getCummulative (const Type histp, const char* option) {

  int nbins(histp->GetNbinsX());

  double val1(histp->GetBinContent(nbins+1));
  double val2(histp->GetBinContent(nbins));
  double err1(histp->GetBinError(nbins+1));
  double err2(histp->GetBinError(nbins));
  histp->SetBinContent(nbins+1, 0);
  histp->SetBinContent(nbins, val1 + val2);
  histp->SetBinError(nbins, sqrt(err1*err1 + err2*err2));

/*eleni
  double Val1(histp->GetBinContent(nbins+1));
  double Val2(histp->GetBinContent(nbins));
  double Err1(histp->GetBinError(nbins+1));
  double Err2(histp->GetBinError(nbins));
*/
  double Val1(histp->GetBinContent(0));
  double Val2(histp->GetBinContent(1));
  double Err1(histp->GetBinError(0));
  double Err2(histp->GetBinError(1));
  histp->SetBinContent(0, 0);
  histp->SetBinContent(1, Val1 + Val2);
  histp->SetBinError(1, sqrt(Err1*Err1 + Err2*Err2));

  Type hcummulp = (Type)histp->Clone("hcummul") ; 
//  hcummulp->Sumw2() ; 
  hcummulp->Reset(0) ; 

  TString mode(option) ; 
  double val(0) ; 
  double err(0) ; 
  for (int ibin = 1; ibin <= nbins; ++ibin) {
    if (mode.Contains("Up")) {
      val += histp->GetBinContent(ibin) ; 
      err = (err*err) + (histp->GetBinError(ibin)*histp->GetBinError(ibin)) ; 
      err = sqrt(err) ; 
      hcummulp->SetBinContent(ibin, val) ; 
      hcummulp->SetBinError(ibin, err) ; 
    }
    else if (mode.Contains("Down")) {
      val += histp->GetBinContent(nbins + 1 - ibin) ; 
      err = (err*err) + (histp->GetBinError(nbins + 1 - ibin)*histp->GetBinError(nbins + 1 - ibin)) ; 
      err = sqrt(err) ; 
      hcummulp->SetBinContent(nbins + 1 - ibin, val) ; 
      hcummulp->SetBinError(nbins + 1 - ibin, err) ; 
    }
    else {
      std::cout << " Error: Option " << option << " not recognised\n" ; 
    }
  }

  Type htot = (Type)histp->Clone("htot") ; 
//  htot->Sumw2() ; 
  for (int ibin = 1; ibin <= nbins; ++ibin) {
    if (mode.Contains("Up")) {
      htot->SetBinContent(ibin, hcummulp->GetBinContent(nbins)) ; 
      htot->SetBinError(ibin, hcummulp->GetBinError(nbins)) ; 
    }
    else if (mode.Contains("Down")) {
      htot->SetBinContent(nbins + 1 - ibin, hcummulp->GetBinContent(1)) ; 
      htot->SetBinError(nbins + 1 - ibin, hcummulp->GetBinError(1)) ; 
    }
    else {
      std::cout << " Error: Option " << option << " not recognised\n" ; 
    }
  }

  TGraphAsymmErrors* geff = new TGraphAsymmErrors(hcummulp, htot, "cl=0.683 b(1,1) mode") ; 
  geff->SetName("geff"+mode) ; 

  return geff ; 
}


