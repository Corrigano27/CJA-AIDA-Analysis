#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TMath.h"
#include "TFile.h"
//#include "Math/WrappedMultiTF1.h"
#include "TF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <functional>

#include "TAxis.h"
#include "TArrayD.h"

Double_t ScaleX(Double_t x)
{
  Double_t v;
  v = 101/100 * x + 1.639; // "linear scaling" function example
  return v;
}
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin()), // new Xmin
              Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}

void ScaleXaxis(TH1 *h, Double_t (*Scale)(Double_t))
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), Scale);
  return;
}

int betapfit(std::string iName){

    std::vector<TH1D *> inputHistogram[6];
    std::vector<TH1D *> inputHistogram_11[6];

    TFile * ifile = TFile::Open( iName.c_str(), "read");

    TFile * ifile_11 = TFile::Open("EnergyHistos_11.root", "read");

    TCanvas *c1 = new TCanvas("c1","Energy Spectra", 1024, 768);
    c1->cd();
    c1->SetLeftMargin(0.19);
    c1->SetBottomMargin(0.17);
    c1->SetTopMargin(0.04);
    c1->SetRightMargin(0.04);

    gStyle->SetOptStat(0);

    char element[][6] = {"Sn","In","Cd","Ag","Pd","Rh"};

    int isotopeStart[6] = {99,97,95,94,93,92};
    int isotopeEnd[6] = {101,100,99,96,95,92}; 

    for (int i = 0; i < 6; i++){
        for (int j = 0; j <= isotopeEnd[i] - isotopeStart[i]; j++){

            char histoDir[100];

            int isotope = isotopeStart[i] + j;

            std::string nucleus = element[i] + std::to_string(isotope);

            //printf("%s",histoDir);
            //if (nucleus == "Sn99" || nucleus == "In97" || nucleus == "Cd95" || nucleus == "Sn100"){
            if (nucleus == "Sn99" || nucleus == "In97" || nucleus == "Cd95"){
                sprintf(histoDir, "/%s/%s%d/%s%ddelayed1pEnergyTotal_AllDSSD_ExMaxSumCorr0",element[i],element[i],isotope,element[i],isotope);
            }
            else{
                sprintf(histoDir, "/%s/%s%d/%s%ddelayed1pEnergyTotal_AllDSSD_ExSumCorr0",element[i],element[i],isotope,element[i],isotope);
            }

            inputHistogram[i].push_back((TH1D*) ifile->Get(histoDir));
            inputHistogram_11[i].push_back((TH1D*) ifile_11->Get(histoDir));

            inputHistogram[i].at(j)->GetXaxis()->SetTitleSize(0.08);
            inputHistogram[i].at(j)->GetXaxis()->SetLabelSize(0.06);
            inputHistogram[i].at(j)->GetXaxis()->SetTitleOffset(0.9);
            inputHistogram[i].at(j)->GetXaxis()->CenterTitle(true);
            inputHistogram[i].at(j)->GetXaxis()->SetTitle("E_{p} (MeV)");
            inputHistogram[i].at(j)->GetYaxis()->CenterTitle(true);
            inputHistogram[i].at(j)->GetYaxis()->SetTitleSize(0.08);
            inputHistogram[i].at(j)->GetYaxis()->SetLabelSize(0.06);
            inputHistogram[i].at(j)->GetYaxis()->SetMaxDigits(6);
            inputHistogram[i].at(j)->GetYaxis()->SetTitleOffset(1.0);
            inputHistogram[i].at(j)->SetMinimum(0);

            inputHistogram_11[i].at(j)->GetXaxis()->SetTitleSize(0.08);
            inputHistogram_11[i].at(j)->GetXaxis()->SetLabelSize(0.06);
            inputHistogram_11[i].at(j)->GetXaxis()->SetTitleOffset(0.9);
            inputHistogram_11[i].at(j)->GetXaxis()->CenterTitle(true);
            inputHistogram_11[i].at(j)->GetXaxis()->SetTitle("E_{p} (MeV)");
            inputHistogram_11[i].at(j)->GetYaxis()->CenterTitle(true);
            inputHistogram_11[i].at(j)->GetYaxis()->SetTitleSize(0.08);
            inputHistogram_11[i].at(j)->GetYaxis()->SetLabelSize(0.06);
            inputHistogram_11[i].at(j)->GetYaxis()->SetMaxDigits(6);
            inputHistogram_11[i].at(j)->GetYaxis()->SetTitleOffset(1.0);
            inputHistogram_11[i].at(j)->SetMinimum(0);

            if (nucleus == "In100" || nucleus == "In98" || nucleus == "Ag94" || nucleus == "Sn101"){
                inputHistogram[i].at(j)->Rebin(2);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.04 MeV");
            }
            if (nucleus == "Rh92"){
                inputHistogram[i].at(j)->Rebin(4);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.08 MeV");
            }
            if (nucleus == "Pd94" || nucleus == "Pd95"){
                inputHistogram[i].at(j)->Rebin(6);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.12 MeV");
            }
            if (nucleus == "In97"){
                inputHistogram[i].at(j)->Rebin(10);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.2 MeV");
            }
            if (nucleus == "Sn100"){
                inputHistogram[i].at(j)->Rebin(10);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.2 MeV");
                inputHistogram_11[i].at(j)->Rebin(10);
                inputHistogram_11[i].at(j)->GetYaxis()->SetTitle("Counts / 0.2 MeV");
                //ScaleXaxis(inputHistogram[i].at(j),ScaleX);
                //inputHistogram[i].at(j)->GetXaxis()->SetTitle("100In Excitation Energy (MeV)");

            }
            if (nucleus == "Sn101"){
              inputHistogram[i].at(j)->Print("all");
              ScaleXaxis(inputHistogram[i].at(j),ScaleX);
              inputHistogram[i].at(j)->GetXaxis()->SetTitle("^{101}In Excitation Energy (MeV)");
            }


            if (nucleus == "In99"){
                inputHistogram[i].at(j)->Rebin(15);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.3 MeV");
            }
            if (nucleus == "Cd98" || nucleus == "Cd99"){
                inputHistogram[i].at(j)->Rebin(20);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.4 MeV");
            }
            if (nucleus == "Cd95" || nucleus == "Sn99"){
                inputHistogram[i].at(j)->Rebin(25);
                inputHistogram[i].at(j)->GetYaxis()->SetTitle("Counts / 0.5 MeV");
            }
            inputHistogram[i].at(j)->Draw();
            inputHistogram_11[i].at(j)->SetLineColor(kRed);
            //inputHistogram_11[i].at(j)->Draw("same");
            char printDir[100];

            sprintf(printDir, "./EnergySpectra/%s%d.pdf",element[i],isotope);

            c1->Print(printDir);

        }
    
    }

   //->Draw();


    return 0;
}

