#pragma once
#include "Solver.h"
#include <Eigen\Eigen>

class EigenLinSolver : public Solver {
public:
  EigenLinSolver();
  ~EigenLinSolver();

  void pardiso_init();
  bool factorize();
  void pardiso_solver();
  void free_numerical_factorization_memory();

private:
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;
  Eigen::SparseMatrix<double> coefMtr;
  void update_coef();
};