#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TString.h"
#include "TLine.h"
#include <vector>
#include <iostream>
using namespace std;

template <typename T>
string NumberToString ( T Number )
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

void Unmatched(){

  TH1::SetDefaultSumw2 (kTRUE);

  float ptRatiobins[9]={30,40,50,60,80,100,120,250,500};
  float etaRatiobins[5]={0.0,2.0,2.5,3.0,4.7};
  //float fractionbins[11]={0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TH2F *h_map = new TH2F("h_map","", 8, ptRatiobins, 4, etaRatiobins);

  TFile *file= new TFile("res/QCD_cat2_singleH.root");

  cout << "opened files" << endl;
  for(int i=0; i<8; i++){
    for(int j=0; j<4; j++){
      TH1::SetDefaultSumw2 (kTRUE);
            
      TH1F*  histo_undef= (TH1F*) file->Get("h_nPart_undef_"+NumberToString(i)+(TString)"_"+NumberToString(j));
      TH1F*  histo_gluon= (TH1F*) file->Get("h_nPart_gluon_"+NumberToString(i)+(TString)"_"+NumberToString(j));
      TH1F*  histo_quark= (TH1F*) file->Get("h_nPart_quark_"+NumberToString(i)+(TString)"_"+NumberToString(j));

      float fraction = histo_undef->GetEntries()/(histo_undef->GetEntries() + histo_gluon->GetEntries() + histo_quark->GetEntries());

      cout << i << "\t" << j << "\t" << fraction << endl;

      h_map->Fill(ptRatiobins[i], etaRatiobins[j], fraction);
    }
  }
      
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  gStyle->SetOptStat(0);
  h_map->Draw("COLZTEXT45");
  h_map->SetTitle("");
  h_map->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_map->GetYaxis()->SetTitle("#eta");
  //  h_map->GetZaxis()->SetTitle("unmatched [%]");
  
  canvas->cd();
  canvas->Print("plots/pythia/map_dijets.pdf");
}
