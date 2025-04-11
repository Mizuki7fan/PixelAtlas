#pragma once

#include "ShellData.h"
#include <CGAL-Interface/CGAL-Interface.h>
#include <LinSysSolver-Interface/Solver.h>

#include <Eigen/Core>

class ParaFun {
public:
  ParaFun(ShellData &data);
  ~ParaFun();

  void AfterMeshImprove();
  double ComputeEnergy(const Eigen::MatrixXd &x, bool whole = false);
  void AdjustShellWeight(double new_weight);
  double BPE(bool is_ip_convrate, bool is_slim_convrate);

private:
  void Init();
  void HandleMinTri();
  void InitArea();
  void SetVirtualTri();
  void PreCalculate();

  void CM(bool is_interp = false);
  void SLIM(bool is_interp = false);
  void UpdateSourceSameT();

  void FunGrid(const Eigen::VectorXd &x);

  void Energy(const Eigen::VectorXd &x, double &energy, bool is_interp = false,
              bool is_whole = true);
  void EnergySource();

  void MaxStep(const Eigen::VectorXd &xx, const Eigen::VectorXd &dd,
               double &step);
  void DetectTmax(const Eigen::VectorXd &x, const Eigen::VectorXd &d,
                  double &tmax);
  void BacktrackingLineSearch(const Eigen::VectorXd &x,
                              const Eigen::VectorXd &d,
                              const Eigen::VectorXd &negetive_grad,
                              double &alpha, bool is_interp = false);
  double GetSmallestPosQuadZero(double a, double b, double c);

  double GetDistance(double s0, double s1, double e0, double e1, double p0,
                     double p1);

  double NewtonEquation(const double &a, const double &b, const double &K);
  void LocalCoordinateInverse(int i, double &p00, double &p01, double &p10,
                              double &p11);
  void LocalCoordinateInverseScaf(int i, double &p00, double &p01, double &p10,
                                  double &p11);

  bool CheckIntersection(const Eigen::VectorXd &pos);

public:
  double energy_mesh_;
  int AV_F_N_H_;
  double time_1_ = 0, time_2_ = 0, time_3_ = 0;
  double density_ = 0;
  std::vector<double> area_;
  double energy_shell_, energy_barrier_;
  double barrer_coef_;

private:
  ShellData &shell_data_;
  int total_num_;
  int F_N_;
  int V_N_;
  static constexpr int kDim = 2;
  std::array<std::vector<int>, 3> F_;
  Eigen::VectorXd position_of_mesh_;
  std::vector<double> area_src_;
  double area_threshold_;

  bool is_first_;
  double Intp_T_Min_;
  double change_to_cm_flag_;
  double bound_distortion_K_;
  std::vector<double> source_p00_, source_p01_, source_p10_, source_p11_;
  std::vector<double> update_p00_, update_p01_, update_p10_, update_p11_;

  std::unique_ptr<Solver> solver_ = nullptr;
  std::vector<int> pardiso_ia_, pardiso_ja_;
  std::vector<double> pardiso_a_, pardiso_b_;

  std::vector<int> id_h00_, id_h01_, id_h02_, id_h03_, id_h04_, id_h05_;
  std::vector<int> id_h11_, id_h12_, id_h13_, id_h14_, id_h15_;
  std::vector<int> id_h22_, id_h23_, id_h24_, id_h25_;
  std::vector<int> id_h33_, id_h34_, id_h35_;
  std::vector<int> id_h44_, id_h45_;
  std::vector<int> id_h55_;

  double threhold;
  double average_length;

  int BE_N_;
  int V_F_N_;
  int AV_F_N_;
  std::vector<int> V_F0_;
  std::vector<int> V_F1_;
  std::vector<int> V_F2_;
  std::vector<int> V_F0_H_;
  std::vector<int> V_F1_H_;
  std::vector<int> V_F2_H_;
  std::vector<int> is_active_;
  std::vector<int> AV_ID_;
  std::vector<int> boundary_vertexID_;

  int cellx_num;
  int celly_num;
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double lengthgrid_x;
  double lengthgrid_y;
  std::vector<std::vector<int>> cell_points;
};
