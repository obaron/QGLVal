#ifndef TOP_CANDIDATE
#define TOP_CANDIDATE

#include <vector>
#include <iostream>
#include "TLorentzVector.h"


class TopCandidate {
public:

  class TopCandidateParticle { 
  public:
    TopCandidateParticle () {;}
    TopCandidateParticle (TLorentzVector v, int num):vec(v),topnum(num) {;}
    TopCandidateParticle (TLorentzVector v, std::string l, int num):vec(v),label(l),topnum(num) {;}
    TopCandidateParticle (TLorentzVector v, std::string l, int num, unsigned b):vec(v),label(l),topnum(num),bitrep(b) {;}
    ~TopCandidateParticle () {;}

    static const unsigned BitB=0x1;
    static const unsigned BitWd1=0x2;
    static const unsigned BitWd2=0x4;
    
    
    TLorentzVector vec;
    std::string label;
    int topnum;
    unsigned bitrep;
  };
  

  TLorentzVector topvec;
  TLorentzVector Wvec;
  std::vector<TopCandidateParticle> particles;
  std::vector<TopCandidateParticle> particles_fit;
  std::string label;
  std::pair< unsigned, unsigned > bitrep; // 1st "W", 2nd "B"

  TopCandidate(TopCandidateParticle p0,TopCandidateParticle p1,TopCandidateParticle p2) {
    particles.push_back(p0);
    particles.push_back(p1);
    particles.push_back(p2);
    topvec = p0.vec + p1.vec + p2.vec;
    Wvec   = p0.vec + p1.vec;    
    label = p0.label + std::string(":") + p1.label + std::string(":") + p2.label;
    bitrep.first = p0.bitrep | p1.bitrep;
    bitrep.second = p2.bitrep;
  };
  ~TopCandidate() {};

  void updateLabeling ();
  bool matched();
  int  nunmatched();
  void reset();
  static bool compare(TopCandidateParticle first, TopCandidateParticle second );

};

#endif
