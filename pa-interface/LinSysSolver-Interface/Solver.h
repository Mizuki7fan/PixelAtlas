#pragma once
#include <vector>

class Solver {
private:
  /* data */
public:
  virtual ~Solver() {};
  virtual void PardisoInit() = 0;
  virtual bool Factorize() = 0;
  virtual void PardisoSolver() {};
  virtual void FreeNumericalFactorizationMemory() {};

  std::vector<double> result_;

  std::vector<int> ia_;
  std::vector<int> ja_;
  std::vector<double> a_;
  std::vector<double> rhs_;

  int nnz_;
  int num_;
};
