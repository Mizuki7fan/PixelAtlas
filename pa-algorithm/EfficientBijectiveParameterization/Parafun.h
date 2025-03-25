#pragma once
#include "LinSysSolverInterface/Solver.h"
#include "ShellData.h"

class Parafun {
public:
  explicit Parafun(const std::unique_ptr<ShellData> &data);
  ~Parafun() {};

public:
  void AfterMeshImprove();
  double ComputeEnergy(const Eigen::MatrixXd &x, bool whole = false);
  void AdjustShellWeight(double new_weight);
  double BPE(bool is_ip_convrate, bool is_slim_convrate);
  double energy_mesh_, energy_barrier_, energy_shell_, energy_all_;
  int AV_F_N_H_;
  double time_1_ = 0, time_2_ = 0, time_3_ = 0;
  std::vector<double> area_;
  double barrier_coef_;

private:
  void Init();
  void InitArea();
  void HandleMinTri();
  void SetVirtualTri();
  void PreCalculate();
  void LocalCoordinateInverse(int i, double &p00, double &p01, double &p10,
                              double &p11);
  void LocalCoordinateInverseScaf(int i, double &p00, double &p01, double &p10,
                                  double &p11);
  void FunGrid(const Eigen::VectorXd &x);
  double GetDistance(double s0, double s1, double e0, double e1, double p0,
                     double p1);
  double NewtonEquation(const double &a, const double &b, const double &K);
  void MaxStep(const Eigen::VectorXd &xx, const Eigen::VectorXd &dd,
               double &step);
  void BacktrackingLineSearch(const Eigen::VectorXd &x,
                              const Eigen::VectorXd &d,
                              const Eigen::VectorXd &negetive_grad,
                              double &alpha, bool is_interp = false);
  double GetSmallestPosQuadZero(double a, double b, double c);
  void DetectTmax(const Eigen::VectorXd &x, const Eigen::VectorXd &d,
                  double &tmax);
  void Energy(const Eigen::VectorXd &x, double &energy, bool is_interp = false,
              bool is_whole = true);
  void Energysource();
  bool CheckIntersection(const Eigen::VectorXd &pos);
  void CM(bool is_interp = false);
  void SLIM(bool is_interp = false);
  void UpdateSourceSameT();

  const std::unique_ptr<ShellData> &shell_data_;
  int total_num_;
  int F_N_, V_N_;
  std::vector<int> F_0_, F_1_, F_2_;
  Eigen::VectorXd position_of_mesh_;
  const int kDim = 2;
  std::unique_ptr<Solver> solver_;
  std::vector<int> pardiso_ia_, pardiso_ja_;
  std::vector<double> pardiso_a_, pardiso_b_;
  bool is_first_;
  double bound_distortion_K_;

  std::vector<double> source_p00_, source_p01_, source_p10_, source_p11_;
  std::vector<double> update_p00_, update_p01_, update_p10_, update_p11_;
  std::vector<double> area_src_;
  double area_threshold_;

  double x_min_, x_max_, y_min_, y_max_;

  int cellx_num_, celly_num_;
  double average_length_, threshold_;
  std::vector<std::vector<int>> cell_points_;

  int AV_F_N_, BE_N_, V_F_N_;
  std::vector<int> V_F_0_, V_F_1_, V_F_2_;
  std::vector<int> V_F0_H_, V_F1_H_, V_F2_H_;
  std::vector<int> is_active_, AV_ID_, boundary_vertexID_;

  double Intp_T_Min_;
  double change_to_cm_flag_;

  std::vector<int> id_h00_, id_h01_, id_h02_, id_h03_, id_h04_, id_h05_;
  std::vector<int> id_h11_, id_h12_, id_h13_, id_h14_, id_h15_;
  std::vector<int> id_h22_, id_h23_, id_h24_, id_h25_;
  std::vector<int> id_h33_, id_h34_, id_h35_;
  std::vector<int> id_h44_, id_h45_;
  std::vector<int> id_h55_;

  double lengthgrid_x_, lengthgrid_y_;
  double threhold_;
};