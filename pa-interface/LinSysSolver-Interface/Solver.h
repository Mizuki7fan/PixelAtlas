#pragma once
#include <vector>

class Solver {
private:
  /* data */
public:
  virtual ~Solver() {};
  virtual void pardiso_init() = 0;
  virtual bool factorize() = 0;
  virtual void pardiso_solver() {};
  virtual void free_numerical_factorization_memory() {};

  std::vector<double> result_;

  std::vector<int> ia_;
  std::vector<int> ja_;
  std::vector<double> a_;
  std::vector<double> rhs_;

  int nnz_;
  int num_;
};
