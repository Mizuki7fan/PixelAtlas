#pragma once

#include "ShellData.h"
#include <CGAL-Interface/CGAL-Interface.h>
#include <LinSysSolver-Interface/Solver.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>

class ParaFun {
public:
  ParaFun(ShellData &data);
  ~ParaFun();

  void AfterMeshImprove();
  double ComputeEnergy(const Eigen::MatrixXd &x, bool whole = false);
  void adjust_shell_weight(double new_weight);
  double BPE(bool is_ip_convrate, bool is_slim_convrate);

private:
  void init();
  void handle_mintri();
  void init_area();
  void setvirtualtri();
  void Pre_calculate();

  void CM(bool is_interp = false);
  void SLIM(bool is_interp = false);
  void Update_source_same_t();

  void fungrid(const Eigen::VectorXd &x);

  void Energy(const Eigen::VectorXd &x, double &energy, bool is_interp = false,
              bool is_whole = true);
  void Energysource();

  void max_step(const Eigen::VectorXd &xx, const Eigen::VectorXd &dd,
                double &step);
  void tmaxdetect(const Eigen::VectorXd &x, const Eigen::VectorXd &d,
                  double &tmax);
  void backtracking_line_search(const Eigen::VectorXd &x,
                                const Eigen::VectorXd &d,
                                const Eigen::VectorXd &negetive_grad,
                                double &alpha, bool is_interp = false);
  double get_smallest_pos_quad_zero(double a, double b, double c);

  double distance(double s0, double s1, double e0, double e1, double p0,
                  double p1);

  double newton_equation(const double &a, const double &b, const double &K);
  void local_coordinate_inverse(int i, double &p00, double &p01, double &p10,
                                double &p11);
  void local_coordinate_inverse_scaf(int i, double &p00, double &p01,
                                     double &p10, double &p11);

  bool check_intersection(const Eigen::VectorXd &pos);

public:
  double energy_mesh;
  int AV_F_N_H;
  double time1 = 0, time2 = 0, time3 = 0;
  double density = 0;
  std::vector<double> area;
  double energy_shell, energy_barrier;
  double barrer_coef;

private:
  ShellData &d_;
  int total_num;
  int F_N;
  int V_N;
  int kDim = 2;
  std::vector<int> F0;
  std::vector<int> F1;
  std::vector<int> F2;
  Eigen::VectorXd position_of_mesh;
  std::vector<double> area_src;
  double area_threshold;

  bool is_first;
  double Intp_T_Min;
  double changetocm_flag;
  double bound_distortion_K;
  std::vector<double> source_p00;
  std::vector<double> source_p01;
  std::vector<double> source_p10;
  std::vector<double> source_p11;
  std::vector<double> update_p00;
  std::vector<double> update_p01;
  std::vector<double> update_p10;
  std::vector<double> update_p11;

  Solver *pardiso;
  std::vector<int> pardiso_ia;
  std::vector<int> pardiso_ja;
  std::vector<double> pardiso_a;
  std::vector<double> pardiso_b;

  double energy_all;

  std::vector<int> id_h00;
  std::vector<int> id_h01;
  std::vector<int> id_h02;
  std::vector<int> id_h03;
  std::vector<int> id_h04;
  std::vector<int> id_h05;
  std::vector<int> id_h11;
  std::vector<int> id_h12;
  std::vector<int> id_h13;
  std::vector<int> id_h14;
  std::vector<int> id_h15;
  std::vector<int> id_h22;
  std::vector<int> id_h23;
  std::vector<int> id_h24;
  std::vector<int> id_h25;
  std::vector<int> id_h33;
  std::vector<int> id_h34;
  std::vector<int> id_h35;
  std::vector<int> id_h44;
  std::vector<int> id_h45;
  std::vector<int> id_h55;

  double threhold;
  double average_length;

  int BE_N;
  int V_F_N;
  int AV_F_N;
  std::vector<int> V_F0;
  std::vector<int> V_F1;
  std::vector<int> V_F2;
  std::vector<int> V_F0_H;
  std::vector<int> V_F1_H;
  std::vector<int> V_F2_H;
  std::vector<int> is_active;
  std::vector<int> AV_ID;
  std::vector<int> boundary_vertexID;

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
