#ifndef FIT_RESULTS
#define FIT_RESULTS

class FitResults { 
public:

  FitResults():
    mass(0.), massW(0.),
    fitmass(0.), fitmassW(0.),
    prob(0.), chisq(999.), cost(999.),
    nsteps(0), nunmatched(0),
    matched(false), converged(false)
  {}
  ~FitResults(){}

  std::string label;
  std::pair<unsigned,unsigned> bitrep; // 1st "W", 2nd "B"
  double mass, massW;
  double fitmass, fitmassW;
  double prob, chisq, cost;
  int nsteps, nunmatched;
  bool matched, converged;

  int topnum;

};

#endif
