#include "../interface/TopCandidate.hh"

bool TopCandidate::matched() { 
  if( label == std::string("Wd2:Wd1:B") || 
      label == std::string("Wd1:Wd2:B") ) {
    return true;
  }
  return false;
};  

int TopCandidate::nunmatched() { 
  int n=0;
  if( particles[0].bitrep > 0x4 ) n++;
  if( particles[1].bitrep > 0x4 ) n++;
  if( particles[2].bitrep > 0x4 ) n++;
  return n;
};  


void TopCandidate::reset() { 
  particles_fit = particles;
  topvec = particles[0].vec + particles[1].vec + particles[2].vec;
  Wvec   = particles[0].vec + particles[1].vec;
};


bool TopCandidate::compare(TopCandidateParticle first, TopCandidateParticle second ) { 
  if (first.vec.Pt() > second.vec.Pt()) return true;
  return false;
}

