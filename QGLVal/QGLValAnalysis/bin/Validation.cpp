#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include <stdio.h>
#include <string.h>
#include "ttDM/TopTagResolved/interface/KinematicFitter.hh"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "QGLVal/QGLValAnalysis/interface/Weights.h"
#include "QGLVal/QGLValAnalysis/interface/MT2Utility.h"
#include "QGLVal/QGLValAnalysis/interface/mt2w_bisect.h"
#include "QGLVal/QGLValAnalysis/interface/mt2bl_bisect.h"
#include "QGLVal/QGLValAnalysis/interface/Mt2Com_bisect.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "ttDM/localQGLikelihoodCalculator/localQGLikelihoodCalculator.h"
#include "ttDM/localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

void callme(){
  std::cout<<" NaN value"<<std::endl;
}

int main(int argc, char **argv) {

  std::cout<<"Let's start"<<endl;

  string sample(argv[1]) ;
  std::cout<<"sample: "<<sample<<endl;

   string path(argv[2]);
  std::cout<<"File list to open: "<<path<<endl;

  string channel(argv[3]);
  std::cout<<"channel: "<<channel<<endl;

  string cat(argv[4]);
  std::cout<<"category:" <<cat<<endl;

  string sys(argv[5]);
  std::cout<<"systematics: "<<sys<<endl;

  string sync(argv[6]);
  std::cout<<"synchro: "<<sync<<endl;

  string isData(argv[7]);
  std::cout<<"isData: "<<isData<<endl;

  //  TString path_ = path + "/trees*.root";                                                                                 
  TString path_ = path ;
 std::cout<<"File to open: "<<path_<<endl;

  std::cout << "Loading file collection from " << path << std::endl;
  TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
  std::cout << "Files found : " << fc.GetNFiles() << std::endl;
  cout << "=========" << sample << endl;

  string reportName = "SelectedEvents_"+channel+"_"+cat+"_"+sample+".txt";
  ofstream fileout;
  fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
  fileout<<"RunNumber EvtNumber Lumi "<<std::endl;

  //  TString outfile = "res/"+sample + "_" +cat+"_"+channel+".root";
  //TFile fout(outfile, "recreate");

  std::cout<<"File to open: "<<path_<<endl;
  TString treePath = "DMTreesDumper/ttDM__noSyst";
                                                                                                       
  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());

  Int_t nEvents = (Int_t)chain.GetEntries();
  std::cout<<"Number of Events: "<<nEvents<< endl;
  //  if(nEvents>5000000) nEvents=5000000;
  nEvents=10000000;
                                                                              
  int sizeMax=50;
  Int_t jetSize, genPartSize, muonSize;
  float passTrig300(0.), passTrig600(0.), passTrig800(0.), passTrigBias(0.);
  Float_t nPV, nGoodPV, nTruePV;

  float rho(0.),w(1.), runNumber(0.), lumiSec(0.), Ht(0.);
  double evtNumber(0.);
  int n_trig(0), n_isBackToBack(0), n_radiation(0);

  float genpartpt[sizeMax], genparteta[sizeMax], genpartphi[sizeMax], genparte[sizeMax], genpartid[sizeMax], genpartstatus[sizeMax], genpartmomid[sizeMax], genpartmass[sizeMax];
  float jetflavour[sizeMax], jete[sizeMax], jetpt[sizeMax], jetphi[sizeMax], jeteta[sizeMax], jetcsv[sizeMax], jetIsTight[sizeMax],jetmult[sizeMax], jetptd[sizeMax], jetaxis2[sizeMax], pileupJetIdptD[sizeMax],pileupJetIdnParticles[sizeMax],pileupJetIdminW[sizeMax], jetqgl[sizeMax],pileupJetIdRMS[sizeMax],pileupJetIdbeta[sizeMax], pileupJetIdbetaClassic[sizeMax], pileupJetIdbetaStar[sizeMax],pileupJetIdbetaStarClassic[sizeMax],pileupJetIddR2Mean[sizeMax],pileupJetIddRMean[sizeMax], pileupJetIddZ[sizeMax],pileupJetIdfrac01[sizeMax],pileupJetIdfrac02[sizeMax],pileupJetIdfrac03[sizeMax],pileupJetIdfrac04[sizeMax],pileupJetIdjetR[sizeMax],pileupJetIdjetRchg[sizeMax],pileupJetIdmajW[sizeMax],pileupJetIdnCharged[sizeMax],pileupJetIdnNeutrals[sizeMax],pileupJetIdpull[sizeMax] ;
   float muone[sizeMax], muonpt[sizeMax], muonphi[sizeMax], muoneta[sizeMax], muonIsTight[sizeMax], muoncharge[sizeMax];
   
   chain.SetBranchAddress("genpart_Pt", genpartpt);
   chain.SetBranchAddress("genpart_Eta", genparteta);
   chain.SetBranchAddress("genpart_Phi", genpartphi);
   chain.SetBranchAddress("genpart_E", genparte);
   chain.SetBranchAddress("genpart_Mass", genpartmass);
   chain.SetBranchAddress("genpart_ID", genpartid);
   chain.SetBranchAddress("genpart_Status", genpartstatus);
   chain.SetBranchAddress("genpart_Mom0ID", genpartmomid);
   chain.SetBranchAddress("genpart_size", &genPartSize);

   chain.SetBranchAddress("jetsAK4Tight_PartonFlavour", jetflavour);
   chain.SetBranchAddress("jetsAK4Tight_E", jete);
   chain.SetBranchAddress("jetsAK4Tight_Pt", jetpt);
   chain.SetBranchAddress("jetsAK4Tight_Phi", jetphi);
   chain.SetBranchAddress("jetsAK4Tight_Eta", jeteta);
   chain.SetBranchAddress("jetsAK4Tight_IsTight", jetIsTight);
   chain.SetBranchAddress("jetsAK4Tight_CSVv2", jetcsv);
   chain.SetBranchAddress("jetsAK4Tight_QGL", jetqgl);
   chain.SetBranchAddress("jetsAK4Tight_mult", jetmult);
   chain.SetBranchAddress("jetsAK4Tight_ptD", jetptd);
   chain.SetBranchAddress("jetsAK4Tight_axis2", jetaxis2);
   chain.SetBranchAddress("jetsAK4Tight_size", &jetSize);

   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdptD", pileupJetIdptD);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdnParticles",pileupJetIdnParticles);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdminW",pileupJetIdminW);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdRMS",pileupJetIdRMS);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdbeta", pileupJetIdbeta);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdbetaClassic", pileupJetIdbetaClassic);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdbetaStar", pileupJetIdbetaStar);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdbetaStarClassic",pileupJetIdbetaStarClassic);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIddR2Mean", pileupJetIddR2Mean);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIddRMean", pileupJetIddRMean);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIddZ", pileupJetIddZ);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdfrac01", pileupJetIdfrac01);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdfrac02",pileupJetIdfrac02);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdfrac03",pileupJetIdfrac03);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdfrac04",pileupJetIdfrac04);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdjetR",pileupJetIdjetR);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdjetRchg", pileupJetIdjetRchg);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdmajW", pileupJetIdmajW);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdnCharged", pileupJetIdnCharged);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdnNeutrals",pileupJetIdnNeutrals);
   chain.SetBranchAddress("jetsAK4Tight_pileupJetIdpull", pileupJetIdpull);

   chain.SetBranchAddress("muonsMedium_E", muone);
   chain.SetBranchAddress("muonsMedium_Pt", muonpt);
   chain.SetBranchAddress("muonsMedium_Phi", muonphi);
   chain.SetBranchAddress("muonsMedium_Eta", muoneta);
   chain.SetBranchAddress("muonsMedium_Charge", muoncharge);
   chain.SetBranchAddress("muonsMedium_size", &muonSize);
   chain.SetBranchAddress("muonsMedium_IsTightMuon", &muonIsTight);
   
   chain.SetBranchAddress("Event_passesHLT_PFHT600_v3", &passTrig600);
   chain.SetBranchAddress("Event_passesHLT_PFHT300_v2", &passTrig300);
   chain.SetBranchAddress("Event_passesHLT_PFHT800_v2", &passTrig800);
   chain.SetBranchAddress("Event_passesHLT_ZeroBias_v2", &passTrigBias);
   chain.SetBranchAddress("Event_Ht", &Ht);
   
   chain.SetBranchAddress("Event_Rho", &rho);
   chain.SetBranchAddress("Event_nPV", &nPV);
   chain.SetBranchAddress("Event_nGoodPV", &nGoodPV);
   chain.SetBranchAddress("Event_nTruePV", &nTruePV);
   chain.SetBranchAddress("Event_RunNumber", &runNumber);
   chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
   chain.SetBranchAddress("Event_EventNumber", &evtNumber);

   TH1F *h_nPV = new TH1F("h_nPV","nPV", 30, 0,30);
   TH1F *h_nTruePV = new TH1F("h_nTruePV","nTruePV", 30, 0,30);
   TH1F *h_nGoodPV = new TH1F("h_nGoodPV","nGoodPV", 30, 0,30);
   TH1F *h_deltaR = new TH1F("h_deltaR", "deltaR", 40, 0, 2);
   
   int a,b;
   int N = 10, M=5;

   TH1F *h_qgl[N][M];
   TH1F *h_ptD[N][M] ;
   TH1F *h_minW[N][M] ;
   TH1F *h_nPart[N][M];
   TH1F *h_ptD_quark[N][M] ;
   TH1F *h_minW_quark[N][M] ;
   TH1F *h_nPart_quark[N][M];
   TH1F *h_ptD_gluon[N][M] ;
   TH1F *h_minW_gluon[N][M] ;
   TH1F *h_nPart_gluon[N][M];
   TH1F *h_ptD_undef[N][M] ;
   TH1F *h_minW_undef[N][M] ;
   TH1F *h_nPart_undef[N][M];
   TH1F *h_qgl_quark[N][M];
   TH1F *h_qgl_gluon[N][M];
   TH1F *h_qgl_undef[N][M];
   TH1F *h_chg[N][M];
   TH1F *h_chg_quark[N][M];
   TH1F *h_chg_gluon[N][M];
   TH1F *h_chg_undef[N][M];
   TH1F *h_neu[N][M];
   TH1F *h_neu_quark[N][M];
   TH1F *h_neu_gluon[N][M];
   TH1F *h_neu_undef[N][M];

   TH1F *h_pt[M];
   TH1F *h_pt_quark[M];
   TH1F *h_pt_gluon[M];
   TH1F *h_pt_undef[M];
   
   char pt[100], pt_quark[100], pt_gluon[100], pt_undef[100];
   char chg[100], chg_quark[100], chg_gluon[100], chg_undef[100];
   char neu[100], neu_quark[100], neu_gluon[100], neu_undef[100];
   char qgl[100], qgl_quark[100], qgl_gluon[100], qgl_undef[100];
   char ptD[100], minW[100], nPart[100];
   char ptD_quark[100], minW_quark[100], nPart_quark[100];
   char ptD_gluon[100], minW_gluon[100], nPart_gluon[100];
   char ptD_undef[100], minW_undef[100], nPart_undef[100];

   for(b=0; b<M; b++){
       sprintf(pt, "h_pt_%d",b);
       h_pt[b] = new TH1F(pt, pt, 400, 0, 800);
       sprintf(pt_gluon, "h_pt_gluon_%d", b);
       h_pt_gluon[b] = new TH1F(pt_gluon, pt_gluon,400, 0, 800);
       sprintf(pt_quark, "h_pt_quark_%d",b);
       h_pt_quark[b] = new TH1F(pt_quark, pt_quark,400, 0, 800);
       sprintf(pt_undef, "h_pt_undef_%d",b);
       h_pt_undef[b] = new TH1F(pt_undef, pt_undef,400, 0, 800);
   }
   for(a=0; a<N; a++){
     for(b=0; b<M; b++){
       
       sprintf(qgl, "h_qgl_%d_%d", a,b);
       h_qgl[a][b] = new TH1F(qgl, qgl, 25, 0, 1);
       sprintf(qgl_gluon, "h_qgl_gluon_%d_%d", a,b);
       h_qgl_gluon[a][b] = new TH1F(qgl_gluon, qgl_gluon,25, 0, 1);
       sprintf(qgl_quark, "h_qgl_quark_%d_%d", a,b);
       h_qgl_quark[a][b] = new TH1F(qgl_quark, qgl_quark,25, 0, 1);
       sprintf(qgl_undef, "h_qgl_undef_%d_%d", a,b);
       h_qgl_undef[a][b] = new TH1F(qgl_undef, qgl_undef,25, 0, 1);

        sprintf(chg, "h_chg_%d_%d", a,b);
       h_chg[a][b] = new TH1F(chg, chg, 40 , 0-.5,40-.5);
       sprintf(chg_gluon, "h_chg_gluon_%d_%d", a,b);
       h_chg_gluon[a][b] = new TH1F(chg_gluon, chg_gluon,40 , 0-.5,40-.5);
       sprintf(chg_quark, "h_chg_quark_%d_%d", a,b);
       h_chg_quark[a][b] = new TH1F(chg_quark, chg_quark,40 , 0-.5,40-.5);
       sprintf(chg_undef, "h_chg_undef_%d_%d", a,b);
       h_chg_undef[a][b] = new TH1F(chg_undef, chg_undef,40 , 0-.5,40-.5);

       sprintf(neu, "h_neu_%d_%d", a,b);
       h_neu[a][b] = new TH1F(neu, neu, 40 , 0-.5,40-.5);
       sprintf(neu_gluon, "h_neu_gluon_%d_%d", a,b);
       h_neu_gluon[a][b] = new TH1F(neu_gluon, neu_gluon,40 , 0-.5,40-.5);
       sprintf(neu_quark, "h_neu_quark_%d_%d", a,b);
       h_neu_quark[a][b] = new TH1F(neu_quark, neu_quark,40 , 0-.5,40-.5);
       sprintf(neu_undef, "h_neu_undef_%d_%d", a,b);
       h_neu_undef[a][b] = new TH1F(neu_undef, neu_undef,40 , 0-.5,40-.5);


       sprintf(ptD, "h_ptD_%d_%d", a,b);
       h_ptD[a][b] = new TH1F(ptD, ptD,50,0,1);
       sprintf(minW, "h_minW_%d_%d", a,b);
       h_minW[a][b]  = new TH1F(minW, minW,50,0,7);
       sprintf(nPart, "h_nPart_%d_%d", a,b);
       h_nPart[a][b] = new TH1F(nPart, nPart,50 , 0-.5,50-.5);

       sprintf(ptD_quark, "h_ptD_quark_%d_%d", a,b);
       h_ptD_quark[a][b] = new TH1F(ptD_quark, ptD_quark,50,0,1);
       sprintf(minW_quark, "h_minW_quark_%d_%d", a,b);
       h_minW_quark[a][b]  = new TH1F(minW_quark, minW_quark,50,0,7);
       sprintf(nPart_quark, "h_nPart_quark_%d_%d", a,b);
       h_nPart_quark[a][b] = new TH1F(nPart_quark, nPart_quark,50 , 0-.5,50-.5);

       sprintf(ptD_gluon, "h_ptD_gluon_%d_%d", a,b);
       h_ptD_gluon[a][b] = new TH1F(ptD_gluon, ptD_gluon,50,0,1);
       sprintf(minW_gluon, "h_minW_gluon_%d_%d", a,b);
       h_minW_gluon[a][b]  = new TH1F(minW_gluon, minW_gluon,50,0,7);
       sprintf(nPart_gluon, "h_nPart_gluon_%d_%d", a,b);
       h_nPart_gluon[a][b] = new TH1F(nPart_gluon, nPart_gluon,50 , 0-.5,50-.5);

       sprintf(ptD_undef, "h_ptD_undef_%d_%d", a,b);
       h_ptD_undef[a][b] = new TH1F(ptD_undef, ptD_undef,50,0,1);
       sprintf(minW_undef, "h_minW_undef_%d_%d", a,b);
       h_minW_undef[a][b]  = new TH1F(minW_undef, minW_undef,50,0,7);
       sprintf(nPart_undef, "h_nPart_undef_%d_%d", a,b);
       h_nPart_undef[a][b] = new TH1F(nPart_undef, nPart_undef,50 , 0-.5,50-.5);
     }
   }

   float p=1;
   float x = (float)(1./nEvents);
   
   if(strcmp (sample.c_str(),"QCD_Pt15to30") == 0) p=(float)(x * (float)(2.6 * 1837410000 * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt30to50") == 0) p=(float)(x * (float)(2.6 * 140932000  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt50to80") == 0) p=(float)(x * (float)(2.6 *19204300.  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt85to120") == 0) p=(float)(x * (float)(2.6 * 2762530.  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt120to170") == 0) p=(float)(x * (float)(2.6 * 471100. * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt170to300") == 0) p=(float)(x * (float)(2.6 * 117276  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt300to470") == 0) p=(float)(x * (float)(2.6 * 7823 * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt470to600") == 0) p=(float)(x * (float)(2.6 * 648.2 * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt600to800") == 0) p=(float)(x * (float)(2.6 * 186.9  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt800to1000") == 0) p=(float)(x * (float)(2.6 * 32.293  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt1000to1400") == 0) p=(float)(x * (float)(2.6 * 9.4183 * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt1400to1800") == 0) p=(float)(x * (float)(2.6 * 0.84265  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt1800to2400") == 0) p=(float)(x * (float)(2.6 * 0.114943  * 1000));
   else if(strcmp (sample.c_str(),"QCD_Pt2400to3200") == 0) p=(float)(x * (float)(2.6 * 0.00682981  * 1000)); 
 
   cout << "weight is" << p << endl;
   //p=1;
   TH1F *h_cutFlow = new TH1F("h_cutFlow", "cutFlow", 5, -0.5, 4.5);

   QGLikelihoodCalculator localQG("/mnt/t3nfs01/data01/shome/grauco/JetMET/CMSSW_7_6_3_patch2/src/QGLVal/QGLValAnalysis/pdfQG_AK4chs_13TeV_v2_PU20bx25_QCD_AllPtBins.root");

   //   if(nEvents>500000) nEvents = 500000;
   
   for(Int_t i=0; i<nEvents; i++ )
     {
     
       if(i%100000==1 ){
	 cout<<"Running on event: "<<i<<endl;
       }
       chain.GetEntry(i);
       
       float ptRatiobins[10]={30,40,50,60,80,100,120,250,500,800};
       float etaRatiobins[5]={0.0,2.0,2.5,3.0,4.7};

       struct TightJets{
	 TLorentzVector vect;
	 float deltaphi;
	 float qgl;
	 float ptD;
	 int nPart;
	 float minW;
	 float chg;
	 float neu;
	 float pt;
	 int jetflavour;
       };
       std::vector<TightJets> SelectedJets;
    
       if((passTrigBias)>0.){	
	 n_trig+= w;
	 if(jetSize>1 && jetpt[0]>30 && jetpt[1]>20){
	   TLorentzVector jet0;
	   jet0.SetPtEtaPhiE(jetpt[0], jeteta[0], jetphi[0], jete[0]);
	   TLorentzVector jet1;
	   jet1.SetPtEtaPhiE(jetpt[1], jeteta[1], jetphi[1], jete[1]);
	   float deltaPhi = jet1.DeltaPhi(jet0);
	   bool areBackToBack = deltaPhi>2.5;
	   float ptaverage = (jet1.Pt()+jet0.Pt())/2;
	   bool thirdjet = 0;
	   if(jetSize>2 && jetpt[2]<(0.3*ptaverage)) thirdjet=1;
	   if(jetSize==2) thirdjet=1;
	   if(areBackToBack==1 && thirdjet==1){
	     TightJets b;
	     b.vect = jet0;
	     b.deltaphi = jet0.DeltaPhi(jet1);
	     b.qgl=localQG.computeQGLikelihood(jet0.Pt(), jet0.Eta(), rho, {(float) jetmult[0], jetptd[0], -log(jetaxis2[0])});
	     b.minW=jetaxis2[0];
	     b.ptD=jetptd[0];
	     b.nPart = jetmult[0];
	     b.chg = pileupJetIdnCharged[0] ;
	     b.neu = pileupJetIdnNeutrals[0];
	     b.pt = jet0.Pt();
	     b.jetflavour=jetflavour[0];
	     SelectedJets.push_back(b);

	     TightJets b2;
	     b2.vect = jet1;
	     b2.deltaphi = jet0.DeltaPhi(jet1);
	     b2.qgl=localQG.computeQGLikelihood(jet1.Pt(), jet1.Eta(), rho, {(float) jetmult[1], jetptd[1], -log(jetaxis2[1])});
	     b2.minW=jetaxis2[1];
	     b2.ptD=jetptd[1];
	     b2.nPart = jetmult[1];
	     b2.chg = pileupJetIdnCharged[1];
	     b2.neu = pileupJetIdnNeutrals[1];
	     b2.pt = jet1.Pt();
	     b2.jetflavour = jetflavour[1];
	     SelectedJets.push_back(b2);
	   }
	     
	   if(thirdjet){
	     n_radiation+=w;
	     if(areBackToBack){
	       n_isBackToBack+=w;
	       if(SelectedJets.size()>1){
		 
		 h_nPV->Fill(nPV,p);
		 h_nGoodPV->Fill(nGoodPV,p);
		 h_nTruePV->Fill(nTruePV,p);
		 
		 for(int j=0;j<4;j++){
		   if(abs((SelectedJets[0].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[0].vect).Eta())<etaRatiobins[j+1]){
		     h_pt[j]->Fill(SelectedJets[0].pt,p);
		   }
		   if(abs((SelectedJets[1].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[1].vect).Eta())<etaRatiobins[j+1]){
		     h_pt[j]->Fill(SelectedJets[1].pt,p);
                   }
		 }
		 
		 for(int i=0;i<9;i++){
		   for(int j=0;j<4;j++){	       
		     if(abs((SelectedJets[0].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[0].vect).Eta())<etaRatiobins[j+1] && ((SelectedJets[1].vect).Pt())>ptRatiobins[i] &&  ((SelectedJets[1].vect).Pt())<ptRatiobins[i+1]){	 
		       h_ptD[i][j]->Fill(SelectedJets[0].ptD,p);
		       h_nPart[i][j]->Fill(SelectedJets[0].nPart,p);
		       h_minW[i][j]->Fill(-log(SelectedJets[0].minW),p);
		       h_qgl[i][j]->Fill(SelectedJets[0].qgl,p);
		       h_chg[i][j]->Fill(SelectedJets[0].chg,p);
		       h_neu[i][j]->Fill(SelectedJets[0].neu,p);
		       
		     }
		     if(abs((SelectedJets[1].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[1].vect).Eta())<etaRatiobins[j+1] && ((SelectedJets[0].vect).Pt())>ptRatiobins[i] &&  ((SelectedJets[0].vect).Pt())<ptRatiobins[i+1]){
		       h_ptD[i][j]->Fill(SelectedJets[1].ptD,p);
		       h_nPart[i][j]->Fill(SelectedJets[1].nPart,p);
		       h_minW[i][j]->Fill(-log(SelectedJets[1].minW),p);
		       h_qgl[i][j]->Fill(SelectedJets[1].qgl,p);
		       h_chg[i][j]->Fill(SelectedJets[1].chg,p);
		       h_neu[i][j]->Fill(SelectedJets[1].neu,p);
		     }
		   } 
		 }
		 
		 float deltaR_jetallpartons[2];
		 float mindeltaR_jetallpartons[2];
		 int genindex[2];
		 float maxpt[2];
		 //		 int maxindex[2];
		 
		 
		 
		 //int index[2];
		 for(int n=0; n<2; n++){
		   deltaR_jetallpartons[n]= -999;
		   mindeltaR_jetallpartons[n]= 999;
		   genindex[n]=999;
		   maxpt[n]=-999;
		   //		   maxindex[n]=999;
		   //index[n]=999;
		   for(int i=0; i<genPartSize; i++){	   
		     if(genpartpt[i]>0.1){
		       if( (genpartid[i]<=5 && genpartid[i]>=-5) || (genpartid[i]<=21.2 && genpartid[i]>=20.9)){
			 if(deltaR((SelectedJets[n].vect).Eta(),(SelectedJets[n].vect).Phi(),genparteta[i],genpartphi[i])<0.4){
			   if(genpartpt[i]>maxpt[n]){
			     maxpt[n]=genpartpt[i];
			     //			     maxindex[n]=i;
			   }
			 }
		       }		       
		       if((genpartstatus[i]<=23.1 && genpartstatus[i]>=22.9)){
			 TLorentzVector genParticle;
			 genParticle.SetPtEtaPhiE(genpartpt[i], genparteta[i], genpartphi[i], genparte[i]);
			 deltaR_jetallpartons[n]=deltaR((SelectedJets[n].vect).Eta(),(SelectedJets[n].vect).Phi(),genparteta[i],genpartphi[i]);
			 if(deltaR_jetallpartons[n] < mindeltaR_jetallpartons[n]){                 
			   genindex[n]=genpartid[i];
			   mindeltaR_jetallpartons[n] = deltaR_jetallpartons[n];
			   //index[n]=i;  
			 }
		       }
		     }
		   }
		   //		   if(genpartid[maxindex[n]]>20.9) cout << "ieooo" << endl;
		 }
		 h_deltaR->Fill(mindeltaR_jetallpartons[0],p);
		 h_deltaR->Fill(mindeltaR_jetallpartons[1],p);
		 
		 for(int j=0;j<4;j++){
		   if(abs((SelectedJets[0].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[0].vect).Eta())<etaRatiobins[j+1]){
		     if(mindeltaR_jetallpartons[0]<0.8){
		       if((genindex[0]<=5 && genindex[0]>=-5)){
			 h_pt_quark[j]->Fill(SelectedJets[0].pt,p);
		       }
		       else if(genindex[0]<=21.1 && genindex[0]>=20.9){
			 h_pt_gluon[j]->Fill(SelectedJets[0].pt,p);
		       }
		     }
		     else{
		       h_pt_undef[j]->Fill(SelectedJets[0].pt,p);
		     }
		   }
		   if(abs((SelectedJets[1].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[1].vect).Eta())<etaRatiobins[j+1]){
		     if(mindeltaR_jetallpartons[1]<0.8){
		       if((genindex[1]<=5 && genindex[1]>=-5)){
			 h_pt_quark[j]->Fill(SelectedJets[1].pt,p);
		       }
		       else if(genindex[1]<=21.1 && genindex[1]>=20.9){
			 h_pt_gluon[j]->Fill(SelectedJets[1].pt,p);
		       }
		     }
		     else{
		       h_pt_undef[j]->Fill(SelectedJets[1].pt,p);
		     }
		   }
		   
		 }
		 
		 for(int i=0;i<9;i++){
		   for(int j=0;j<5;j++){
		     if(abs((SelectedJets[0].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[0].vect).Eta())<etaRatiobins[j+1] && ((SelectedJets[1].vect).Pt())>ptRatiobins[i] &&  ((SelectedJets[1].vect).Pt())<ptRatiobins[i+1]){
		       //if((genpartpt[index[0]]/(SelectedJets[0].vect).Pt())<2 && ((genpartpt[index[0]]/(SelectedJets[0].vect).Pt()) > 0.5)){
		       if(mindeltaR_jetallpartons[0]<0.4){
			 // if((( SelectedJets[0].jetflavour<=5 && SelectedJets[0].jetflavour>=-5 && SelectedJets[0].jetflavour!=0))){ 
			 if(((genindex[0]<=5 && genindex[0]>=-5))){
			 //if (genpartid[maxindex[0]]<=5 && genpartid[maxindex[0]]>=-5){
			   h_ptD_quark[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_quark[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_quark[i][j]->Fill(-log(SelectedJets[0].minW),p);                                   
			   h_qgl_quark[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_quark[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_quark[i][j]->Fill(SelectedJets[0].neu,p); 
			 }
			 
			 else if(genindex[0]<=21.1 && genindex[0]>=20.9){
			   //else if (genpartid[maxindex[0]]<=21.1 && genpartid[maxindex[0]]>=20.9){
			 //else if((SelectedJets[0].jetflavour<=21.1 && SelectedJets[0].jetflavour>=20.9)){
			   h_ptD_gluon[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_gluon[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_gluon[i][j]->Fill(-log(SelectedJets[0].minW),p);
			   h_qgl_gluon[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_gluon[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_gluon[i][j]->Fill(SelectedJets[0].neu,p); 
			 }
			 else{
			   h_ptD_undef[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_undef[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_undef[i][j]->Fill(-log(SelectedJets[0].minW),p);
			   h_qgl_undef[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_undef[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_undef[i][j]->Fill(SelectedJets[0].neu,p);
			   }
		       }
		       else{
			 /*if((SelectedJets[0].jetflavour<=5 && SelectedJets[0].jetflavour>=-5 && SelectedJets[0].jetflavour!=0)){
			   h_ptD_quark[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_quark[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_quark[i][j]->Fill(-log(SelectedJets[0].minW),p);
			   h_qgl_quark[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_quark[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_quark[i][j]->Fill(SelectedJets[0].neu,p);

			 }

			 else if((SelectedJets[0].jetflavour<=21.1 && SelectedJets[0].jetflavour>=20.9)){
			   h_ptD_gluon[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_gluon[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_gluon[i][j]->Fill(-log(SelectedJets[0].minW),p);
			   h_qgl_gluon[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_gluon[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_gluon[i][j]->Fill(SelectedJets[0].neu,p);

			 }
			 
			 else if(SelectedJets[0].jetflavour==0){
			   h_ptD_undef[i][j]->Fill(SelectedJets[0].ptD,p);
			   h_nPart_undef[i][j]->Fill(SelectedJets[0].nPart,p);
			   h_minW_undef[i][j]->Fill(-log(SelectedJets[0].minW),p);
			   h_qgl_undef[i][j]->Fill(SelectedJets[0].qgl,p);
			   h_chg_undef[i][j]->Fill(SelectedJets[0].chg,p);
			   h_neu_undef[i][j]->Fill(SelectedJets[0].neu,p); 
			   }
			 */
			 h_ptD_undef[i][j]->Fill(SelectedJets[0].ptD,p);
			 h_nPart_undef[i][j]->Fill(SelectedJets[0].nPart,p);
			 h_minW_undef[i][j]->Fill(-log(SelectedJets[0].minW),p);
			 h_qgl_undef[i][j]->Fill(SelectedJets[0].qgl,p);
			 h_chg_undef[i][j]->Fill(SelectedJets[0].chg,p);
			 h_neu_undef[i][j]->Fill(SelectedJets[0].neu,p);
			 
		        }
		       // }
		       // }
		     }
		     if(abs((SelectedJets[1].vect).Eta())>etaRatiobins[j] && abs((SelectedJets[1].vect).Eta())<etaRatiobins[j+1] && ((SelectedJets[0].vect).Pt())>ptRatiobins[i] &&  ((SelectedJets[0].vect).Pt())<ptRatiobins[i+1]){                                             
		       //if((genpartpt[index[1]]/(SelectedJets[1].vect).Pt())<2 && ((genpartpt[index[1]]/(SelectedJets[1].vect).Pt()) > 0.5)){ 
		       if(mindeltaR_jetallpartons[1]<0.4){
			 if((genindex[1]<=5 && genindex[1]>=-5)){
		       //if((SelectedJets[1].jetflavour<=5 && SelectedJets[1].jetflavour>=-5 && SelectedJets[1].jetflavour!=0)){ 
			   //if (genpartid[maxindex[1]]<=5 && genpartid[maxindex[1]]>=-5){
			   h_ptD_quark[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_quark[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_quark[i][j]->Fill(-log(SelectedJets[1].minW),p);                                     
			   h_qgl_quark[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_quark[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_quark[i][j]->Fill(SelectedJets[1].neu,p);
			 }
			 // else if((SelectedJets[1].jetflavour<=21.1 && SelectedJets[1].jetflavour>=20.9)){ 
		       else if((genindex[1]<=21.1 && genindex[1]>=20.9)){
			   //else if (genpartid[maxindex[1]]<=21.1 && genpartid[maxindex[1]]>=20.9){
			   h_ptD_gluon[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_gluon[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_gluon[i][j]->Fill(-log(SelectedJets[1].minW),p);
			   h_qgl_gluon[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_gluon[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_gluon[i][j]->Fill(SelectedJets[1].neu,p);
			 }
			 
			 else{
			   h_ptD_undef[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_undef[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_undef[i][j]->Fill(-log(SelectedJets[1].minW),p);
			   h_qgl_undef[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_undef[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_undef[i][j]->Fill(SelectedJets[1].neu,p);
			   }
		       }
		       else{
			 
			 /*		 if((SelectedJets[1].jetflavour<=5 && SelectedJets[1].jetflavour>=-5 && SelectedJets[1].jetflavour!=0)){
			   
			   h_ptD_quark[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_quark[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_quark[i][j]->Fill(-log(SelectedJets[1].minW),p);
			   h_qgl_quark[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_quark[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_quark[i][j]->Fill(SelectedJets[1].neu,p);
			 }

			 else if((SelectedJets[1].jetflavour<=21.1 && SelectedJets[1].jetflavour>=20.9)){
			   
			   h_ptD_gluon[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_gluon[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_gluon[i][j]->Fill(-log(SelectedJets[1].minW),p);
			   h_qgl_gluon[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_gluon[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_gluon[i][j]->Fill(SelectedJets[1].neu,p);
			 }

			 else if((SelectedJets[1].jetflavour==0)){
			   
			   h_ptD_undef[i][j]->Fill(SelectedJets[1].ptD,p);
			   h_nPart_undef[i][j]->Fill(SelectedJets[1].nPart,p);
			   h_minW_undef[i][j]->Fill(-log(SelectedJets[1].minW),p);
			   h_qgl_undef[i][j]->Fill(SelectedJets[1].qgl,p);
			   h_chg_undef[i][j]->Fill(SelectedJets[1].chg,p);
			   h_neu_undef[i][j]->Fill(SelectedJets[1].neu,p);
			   }
			 */
			h_ptD_undef[i][j]->Fill(SelectedJets[1].ptD,p);
			 h_nPart_undef[i][j]->Fill(SelectedJets[1].nPart,p);
			 h_minW_undef[i][j]->Fill(-log(SelectedJets[1].minW),p);
			 h_qgl_undef[i][j]->Fill(SelectedJets[1].qgl,p);
			 h_chg_undef[i][j]->Fill(SelectedJets[1].chg,p);
			 h_neu_undef[i][j]->Fill(SelectedJets[1].neu,p); 
			 
			 
			 // }
		     	 }
			 }		
		   }
		   }
	       }	 
	     }
	   }
	 }
       }
     }
   
   TString outfile = "res/"+sample + "_" +cat+"_"+channel+".root";
   TFile fout(outfile, "recreate");
   
   for(b=0; b<M; b++){
     
     h_pt[b]->Write();
     h_pt_gluon[b]->Write();
     h_pt_quark[b]->Write();
     h_pt_undef[b]->Write();
   }
   
   for(a=0; a<N; a++){
     for(b=0; b<M; b++){
       
       h_qgl[a][b]->Write();
       h_qgl_gluon[a][b]->Write();
       h_qgl_quark[a][b]->Write();
       h_qgl_undef[a][b]->Write();

        
       h_chg[a][b]->Write();
       h_chg_gluon[a][b]->Write();
       h_chg_quark[a][b]->Write();
       h_chg_undef[a][b]->Write();

        
       h_neu[a][b]->Write();
       h_neu_gluon[a][b]->Write();
       h_neu_quark[a][b]->Write();
       h_neu_undef[a][b]->Write();

       h_ptD_quark[a][b]->Write();
       h_minW_quark[a][b]->Write();
       h_nPart_quark[a][b]->Write();

       h_ptD_gluon[a][b]->Write();
       h_minW_gluon[a][b]->Write();
       h_nPart_gluon[a][b]->Write();

       h_ptD_undef[a][b]->Write();
       h_minW_undef[a][b]->Write();
       h_nPart_undef[a][b]->Write();

       h_ptD[a][b]->Write();
       h_minW[a][b]->Write();
       h_nPart[a][b]->Write();
      
     }
   }

   h_cutFlow->SetBinContent(1,nEvents);
   h_cutFlow->GetXaxis()->SetBinLabel(1,"no selection");
   h_cutFlow->SetBinContent(2,n_trig);
   h_cutFlow->GetXaxis()->SetBinLabel(2,"trigger");
   h_cutFlow->SetBinContent(3,n_radiation);
   h_cutFlow->GetXaxis()->SetBinLabel(3,"third jet");
   h_cutFlow->SetBinContent(4,n_isBackToBack);
   h_cutFlow->GetXaxis()->SetBinLabel(4,"DeltaPhi");
   h_cutFlow->Write();
   
   h_nPV->Write();
   h_nGoodPV->Write();
   h_nTruePV->Write();
   h_deltaR->Write();

   fout.Close();
   fileout.close();
   
   std::cout<< "---> "<<sample<<std::endl;
   std::cout<< "Number of events           : "<<nEvents<<std::endl;
   std::cout<< "Events after trigger cut   : "<<n_trig<<std::endl;
   std::cout<< "Events after third jet cut   : "<<n_radiation<<std::endl;
   std::cout<< "Events after deltaphi cut   : "<<n_isBackToBack<<std::endl;

}
