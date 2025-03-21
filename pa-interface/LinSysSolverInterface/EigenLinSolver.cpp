#include "EigenLinSolver.h"

EigenLinSolver::EigenLinSolver() {}

EigenLinSolver::~EigenLinSolver() {}

void EigenLinSolver::pardiso_init() {}

bool EigenLinSolver::factorize() {
  update_coef();
  simplicialLDLT.compute(coefMtr);
  return true;
}

void EigenLinSolver::pardiso_solver() {
  Eigen::Map<Eigen::VectorXd> b(rhs_.data(), rhs_.size(), 1);

  Eigen::VectorXd res_ = simplicialLDLT.solve(b);
  result_.resize(num_);
  for (std::size_t i = 0; i < rhs_.size(); i++) {
    result_[i] = res_[i];
  }
}

void EigenLinSolver::free_numerical_factorization_memory() {}

void EigenLinSolver::update_coef() {
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletvec;
  for (int i = 0; i < num_; i++) {
    int j = ia_[i];
    tripletvec.emplace_back(i, i, a_[j]);
    for (j = ia_[i] + 1; j < ia_[i + 1]; j++) {
      tripletvec.emplace_back(i, ja_[j], a_[j]);
      tripletvec.emplace_back(ja_[j], i, a_[j]);
    }
  }
  coefMtr.resize(num_, num_);
  coefMtr.setFromTriplets(tripletvec.begin(), tripletvec.end());
}
