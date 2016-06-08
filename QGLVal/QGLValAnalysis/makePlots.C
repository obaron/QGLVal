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

TString histos_JetHT[9] = {"h_nPart_0_0", "h_nPart_1_0", "h_nPart_2_0", "h_nPart_3_0","h_nPart_4_0","h_nPart_5_0","h_nPart_6_0","h_nPart_7_0","h_nPart_8_0"};
TString histos_QCD_pythia[27] = {"h_nPart_quark_0_0", "h_nPart_gluon_0_0", "h_nPart_undef_0_0", "h_nPart_quark_1_0", "h_nPart_gluon_1_0", "h_nPart_undef_1_0", "h_nPart_quark_2_0", "h_nPart_gluon_2_0", "h_nPart_undef_2_0", "h_nPart_quark_3_0", "h_nPart_gluon_3_0", "h_nPart_undef_3_0", "h_nPart_quark_4_0", "h_nPart_gluon_4_0", "h_nPart_undef_4_0","h_nPart_quark_5_0", "h_nPart_gluon_5_0", "h_nPart_undef_5_0", "h_nPart_quark_6_0", "h_nPart_gluon_6_0", "h_nPart_undef_6_0", "h_nPart_quark_7_0", "h_nPart_gluon_7_0", "h_nPart_undef_7_0", "h_nPart_quark_8_0", "h_nPart_gluon_8_0", "h_nPart_undef_8_0"     };

void makePlots(){

  TH1::SetDefaultSumw2 (kTRUE);
  for(int i=0; i<8; i++){
    TH1::SetDefaultSumw2 (kTRUE);
    TFile *JetHT = new TFile("res/ZeroBias_cat2_singleH.root"); 
    TFile *QCD_pythia= new TFile("res/QCD_cat2_singleH.root");

    float lumi = 2.6;

    TH1F*  h_JetHT= (TH1F*) JetHT->Get(histos_JetHT[i]);
    TH1F*  h_QCD_pythia_quark= (TH1F*) QCD_pythia->Get(histos_QCD_pythia[3*i]);
    TH1F*  h_QCD_pythia_gluon= (TH1F*) QCD_pythia->Get(histos_QCD_pythia[3*i+1]);
    TH1F*  h_QCD_pythia_undef= (TH1F*) QCD_pythia->Get(histos_QCD_pythia[3*i+2]);
   
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);
    gStyle->SetOptStat(0);

    TPad *pad1 = new TPad("pad1", "pad1", 0., 0.35, 1.,1.);
    pad1->SetBottomMargin(0.);
    pad1->Draw();
    pad1->cd();
    h_QCD_pythia_quark->SetFillColor(kAzure-8);
    h_QCD_pythia_gluon->SetFillColor(kRed-6);
    h_QCD_pythia_undef->SetFillColor(kGray+1);
    h_QCD_pythia_quark->SetLineColor(kBlack);
    h_QCD_pythia_gluon->SetLineColor(kBlack);
    h_QCD_pythia_undef->SetLineColor(kBlack);
    
    THStack *hs = new THStack("hs","");
    hs->Add(h_QCD_pythia_undef);
    hs->Add(h_QCD_pythia_gluon);
    hs->Add(h_QCD_pythia_quark);

    float nData = h_JetHT->Integral();
    TH1F* h_hs_last= (TH1F*)hs->GetStack()->Last();
    float nBkg = h_hs_last->Integral();
    float sf_norm = nData/nBkg;

    h_QCD_pythia_quark->Scale(sf_norm);
    h_QCD_pythia_gluon->Scale(sf_norm);
    h_QCD_pythia_undef->Scale(sf_norm);
    h_hs_last->Scale(sf_norm); 
    h_JetHT->Sumw2();
    h_hs_last->Sumw2();

    hs->Modified();
  
    hs->SetMaximum(1.2*max(h_JetHT->GetMaximum(),hs->GetMaximum()));
    hs->Draw("histo");
    hs->SetMinimum(0.1);
    hs->GetYaxis()->SetTitle(Form("Events/%.2f",(h_JetHT->GetBinWidth(2))));
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetTitleOffset(0.8);
    
    canvas->Update();
    
    h_JetHT->SetMarkerStyle(20);
    h_JetHT->SetMarkerSize(1); 
    h_JetHT->SetLineColor(kBlack);
    h_JetHT->SetMarkerColor(kBlack);
    h_JetHT->DrawCopy("same ep");   
    hs->SetMaximum(1.2*max(h_JetHT->GetMaximum(),hs->GetMaximum()));

    TLatex latex;
    TString event ="Dijets";
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.65, 0.82,event);

    TLegend *leg= new TLegend(0.63,0.65,0.78,0.8); 
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_JetHT,"Data","lep");
    leg->AddEntry(h_QCD_pythia_quark,"quark","f");
    leg->AddEntry(h_QCD_pythia_gluon,"gluon","f");
    leg->AddEntry(h_QCD_pythia_undef,"undefined","f");
    leg->Draw();

    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold                                                               
    bool writeExtraText = true;
    TString extraText   = "Work in progress";
    float extraTextFont = 52;  // default is helvetica-italics                                                           
    // text sizes and text offsets with respect to the top frame in unit of the top margin size                          
    float lumiTextSize     = 0.6;
    float lumiTextOffset   = 0.2;
    float cmsTextSize      = 0.55;
    float cmsTextOffset    = 0.1;  // only used in outOfFrame version                                                    
    float relPosX    = 0.045;
    float relPosY    = 0.035;
    float relExtraQCD_pythia = 1.2;
    // ratio of "CMS" and extra text size                                                                                
    float extraOverCmsTextSize  = 0.76;
    TString lumi_13TeV = "2.6 fb^{-1}, #sqrt{s} = 13 TeV ";
    TString lumi_8TeV  = "19.7 fb^{-1}";
    TString lumi_7TeV  = "5.1 fb^{-1}";
    TString lumiText;
    lumiText = "#sqrt{s} = 13 TeV ";
    float t = pad1->GetTopMargin();
    float b = pad1->GetBottomMargin();
    float r = pad1->GetRightMargin();
    float l = pad1->GetLeftMargin();
    //TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(lumiTextSize*t*1.3);
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumi_13TeV);
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11);
    latex.SetTextSize(cmsTextSize*t);
    latex.DrawLatex(l+0.03,1-t+lumiTextOffset*t-0.09,cmsText);
    latex.SetTextFont(extraTextFont);
    latex.SetTextSize(extraTextSize*t);
    latex.DrawLatex(l+0.03, 1-t+lumiTextOffset*t-0.09-0.06, extraText);

    float ptRatioBins[10]={30,40,50,60,80,100,120,250,500,800};
    float etaRatioBins[5]={0.0,2.0,2.5,3.0,4.7};

    int position_lu = histos_JetHT[i].Last('_');
    TString temp( histos_JetHT[i](0,position_lu-1) );
    int position_pu = temp.Last('_');
    
    int value_pt = 0;
    if(position_pu+1==position_lu-1){ 
      int position_pt= (int) ((histos_JetHT[i].Length()))-3;
      value_pt = (char) (histos_JetHT[i])(position_pt) - '0';
    } else {
      temp = histos_JetHT[i](position_pu+1,position_lu-1);
      value_pt = temp.Atoi();
      }
    
    TString binningpt = Form("%.1f GeV < p_{T} < %.1f GeV", ptRatioBins[value_pt], ptRatioBins[value_pt+1]);
    int position_eta = (int) ((histos_JetHT[i].Length()))-1;
    int value_eta = (char) (histos_JetHT[i])(position_eta) - '0'; 
    TString binningeta = Form("%.1f < |#eta| < %.1f", etaRatioBins[value_eta], etaRatioBins[value_eta+1]);
    
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.6, 0.6,binningpt);
    latex.DrawLatex(0.6, 0.55,binningeta);
     
    h_JetHT->SetMaximum(1.3*max(h_JetHT->GetMaximum(),hs->GetMaximum()));
    h_JetHT->SetMinimum(0.1);
    canvas->Update();
    canvas->cd();
    
    TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.3);
    pad2->SetBottomMargin(0.4);
    pad2->SetTopMargin(0.);
    pad2->Draw();
    pad2->cd();
    
    h_JetHT->SetStats(0);
    
    THStack *hs_1 = new THStack("hs_0","");
    hs_1->Add(h_QCD_pythia_undef);
    hs_1->Add(h_QCD_pythia_gluon);
    hs_1->Add(h_QCD_pythia_quark);

    float nData_1 = h_JetHT->Integral();
    TH1F* h_hs_1_last= (TH1F*)hs_1->GetStack()->Last();
    float nBkg_1 = h_hs_1_last->Integral();
    float sf_norm_1 = nData_1/nBkg_1;

    h_QCD_pythia_quark->Scale(sf_norm);
    h_QCD_pythia_gluon->Scale(sf_norm);
    h_QCD_pythia_undef->Scale(sf_norm);
    h_JetHT->Sumw2();
    h_hs_1_last->Sumw2();

    canvas->Update();
 
    h_JetHT->Divide(h_hs_1_last);
    
    canvas->Update();
    h_JetHT->SetMarkerStyle(20);
    h_JetHT->GetYaxis()->SetTitle("#frac{DATA}{MC}");
    h_JetHT->GetYaxis()->SetTitleSize(0.04);
    h_JetHT->GetYaxis()->SetTitleOffset(1.);
    h_JetHT->GetYaxis()->CenterTitle();

    h_JetHT->SetMaximum(2.1);
    h_JetHT->SetMinimum(0);
    h_JetHT->Draw("ep");

    h_JetHT->SetTitle("");
    h_JetHT->GetXaxis()->SetTitle("#sigma_{2}");
    h_JetHT->GetYaxis()->SetTitleSize(0.3);
    h_JetHT->GetYaxis()->SetTitleOffset(0.5);
    h_JetHT->GetXaxis()->SetLabelSize(0.08);
    h_JetHT->GetXaxis()->SetTitleSize(0.10);
    h_JetHT->GetYaxis()->SetTitleSize(0.08);
    h_JetHT->GetYaxis()->SetLabelSize(0.06);
    h_JetHT->GetYaxis()->SetTitleOffset(0.40);

    canvas->Update();
    
    TLine *line=new TLine(canvas->GetUxmin(),1.0,7.,1.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->SetLineStyle(1);
    line->Draw();

    TLine *line2=new TLine(canvas->GetUxmin(),0.5,7.,0.5);
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(1);
    line2->SetLineStyle(2);
    line2->Draw();

    TLine *line3=new TLine(canvas->GetUxmin(),1.5,7.,1.5);
    line3->SetLineColor(kBlack);
    line3->SetLineWidth(1);
    line3->SetLineStyle(2);
    line3->Draw();

    canvas->cd();
    canvas->Print("plots/pythia/"+histos_JetHT[i]+(TString)"_dijets.pdf");

  }
}
