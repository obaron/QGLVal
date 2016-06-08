#ifndef KINEMATIC_FITTER
#define KINEMATIC_FITTER

/*

single top mass constrain to get teh right combinations

S term
-------
3 jet 4 vectors
deltas to those componenents
covariance matrix for measured params

constraint
-----------
mass constraint, here only depends on measured params f(y)
need derivative matrix 'A', 1x12 for 1 constraint eqn



 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "JetResolutions.hh"
#include "TopCandidate.hh"
#include "FitResults.hh"

class KinematicFitter { 

public:

  JetResolutions jetres;
  FitResults     fitstat;  

  const double MTop = 173.2;
  //const double MTop = 175.;
  const double MW   = 81.385;

  KinematicFitter () {
    y.ResizeTo(12);
    dy.ResizeTo(12);
    dy_last.ResizeTo(12);
    V.ResizeTo(12,12);
    B.ResizeTo(2,12);
  };
  ~KinematicFitter () {};

  bool fit( TopCandidate & topcand, FitResults & fitstat); 

  
  private:

  TVectorD y;
  TVectorD dy, dy_last;
  TMatrixD V;
  TMatrixD B;
    
  //  std::pair<unsigned, unsigned> bitrep;
  //  std::string label;
  std::vector<TopCandidate::TopCandidateParticle> the4vecs;

  void set_y( TVectorD & y, TopCandidate & combo );
  void zero_dy( TVectorD & dy );
  void setB( TMatrixD & B, TopCandidate & combo );
  void setV( TMatrixD & V, JetResolutions & myres, TopCandidate & combo );
  void updateVectors( TopCandidate & combo, TVectorD & dy );
};
  
  
#endif
