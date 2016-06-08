//================================================================================================
// Demo usage of resolved top tagger MVA
//________________________________________________________________________________________________

#include <TLorentzVector.h>
#include <TMVA/Reader.h>
#include <string>
#include <iostream>
#include "ttDM/TopTagResolved/interface/KinematicFitter.hh"


int main(int argc, char **argv)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialization
  //
  //*********************************

  //
  // Set up MVA reader
  //

  // spectator variables, not used for MVA evaluation
  int isSig, b_mis, w_mis, wb_mis;
  float mtop;
  // MVA input variables
  float bdt_qgid1, bdt_qgid2;
  float bdt_dphij1b, bdt_dphij2b, bdt_drj1b, bdt_drj2b;
  float bdt_bjcsv, bdt_jet1csv, bdt_jet2csv;
  float bdt_prob;
  TMVA::Reader res_topmvaReader("");
  res_topmvaReader.AddSpectator("isSig",  &isSig);
  res_topmvaReader.AddSpectator("b_mis",  &b_mis);
  res_topmvaReader.AddSpectator("w_mis",  &w_mis);
  res_topmvaReader.AddSpectator("wb_mis", &wb_mis);
  res_topmvaReader.AddSpectator("mtop",   &mtop);
  res_topmvaReader.AddVariable("qgid1",   &bdt_qgid1);    // QGL for one of the W-jet candidates
  res_topmvaReader.AddVariable("qgid2",   &bdt_qgid2);    // QGL for the other W-jet candidate
  res_topmvaReader.AddVariable("dphij1b", &bdt_dphij1b);  // |deltaPhi| between W-jet #1 and b-jet
  res_topmvaReader.AddVariable("dphij2b", &bdt_dphij2b);  // |deltaPhi| between W-jet #2 and b-jet
  res_topmvaReader.AddVariable("drj1b",   &bdt_drj1b);    // |deltaR| between W-jet #1 and b-jet
  res_topmvaReader.AddVariable("drj2b",   &bdt_drj2b);    // |deltaR| between W-jet #2 and b-jet
  res_topmvaReader.AddVariable("bjcsv",   &bdt_bjcsv);    // CSVv2+IVF value of b-jet
  res_topmvaReader.AddVariable("jet1csv", &bdt_jet1csv);  // CSVv2+IVF value of W-jet #1
  res_topmvaReader.AddVariable("jet2csv", &bdt_jet2csv);  // CSVv2+IVF value of W-jet #2
  res_topmvaReader.AddVariable("prob",    &bdt_prob);     // probability of kinematic fit

  const std::string cmssw_base = getenv("CMSSW_BASE");
  const std::string weightsfile = cmssw_base + std::string("/src/TopTagger/Resolved/data/toptrainingbits_prob.root_BDTG.weights.xml");
  res_topmvaReader.BookMVA("BDTG", weightsfile.c_str()); 

  //
  // Kinematic fitter
  //
  KinematicFitter fitter;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Stuff that goes into event loop
  //
  //*********************************

  //
  // Example tri-jet combination
  // The MVA training is done such that the "b-jet" is the jet in the triplet
  // that has the highest CSVv2+IVF value. I've called this one "jet3" in the example.
  //
  TLorentzVector jet1;
  jet1.SetPtEtaPhiM(62.698604, 0.9635512, -1.307790, 9.1356649);
  float jet1csv  = 0.4292013;
  float jet1qgid = 0.9191989;

  TLorentzVector jet2;
  jet2.SetPtEtaPhiM(39.438114, 1.5529232, -2.586280, 5.3947892);
  float jet2csv  = 0.5377913;
  float jet2qgid = 0.9061799;

  TLorentzVector jet3;
  jet3.SetPtEtaPhiM(61.190555, 2.2724466, -1.735949, 14.285264);
  float jet3csv = 0.5930559;

  //
  // Perform kinematic fit
  //

  // initialize object that stores fit results
  FitResults fitres;
  fitres.converged = false;
  fitres.prob      = 0.;
  fitres.chisq     = 999.;
  fitres.cost      = 999.;
  fitres.fitmass   = 0.;
  fitres.fitmassW  = 0.;

  TopCandidate::TopCandidateParticle wjet1(jet1, std::string("unmatched"), 3, 0);
  TopCandidate::TopCandidateParticle wjet2(jet2, std::string("unmatched"), 3, 0);
  TopCandidate::TopCandidateParticle bjet (jet3, std::string("unmatched"), 3, 0);

  // The constructor, TopCandidate(j1,j2,j3), assumes "j3" corresponds to the b-jet while "j1" and "j2" are the W-jets
  TopCandidate combo(wjet1, wjet2, bjet);
  combo.reset();
  fitter.fit(combo, fitres);

  if(fitres.converged) { std::cout << "Fit converged! Top quark probability = " << fitres.prob << std::endl; }
  else                 { std::cout << "Fit did not converge!" << std::endl; }


  //
  // Set the inputs for the MVA
  //
  bdt_qgid1   = jet1qgid;
  bdt_qgid2   = jet2qgid;
  bdt_dphij1b = fabs(jet1.DeltaPhi(jet3));
  bdt_dphij2b = fabs(jet2.DeltaPhi(jet3));
  bdt_drj1b   = jet1.DeltaR(jet3);
  bdt_drj2b   = jet2.DeltaR(jet3);
  bdt_bjcsv   = jet3csv;
  bdt_jet1csv = jet1csv;
  bdt_jet2csv = jet2csv;
  bdt_prob    = fitres.prob;

  //
  // Compute the MVA value
  //
  std::cout << "MVA value = " << res_topmvaReader.EvaluateMVA("BDTG") << std::endl;

  return 0;
}
