#pragma once
#include "Solver.h"
#include <Eigen\Eigen>

class EigenLinSolver : public Solver {
public:
  EigenLinSolver();
  ~EigenLinSolver();

  void PardisoInit();
  bool Factorize();
  void PardisoSolver();
  void FreeNumericalFactorizationMemory();

private:
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicial_LDLT_;
  Eigen::SparseMatrix<double> coef_matrix_;
  void UpdateCoef();
};