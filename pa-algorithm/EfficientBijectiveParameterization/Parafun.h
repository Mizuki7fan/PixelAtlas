#pragma once
#include "LinSysSolverInterface/Solver.h"
#include "ShellData.h"
#include <Eigen/Geometry>
#include <iostream>

class Parafun {

public:
  Parafun(ShellData &data) //
      : d_(data), solver(nullptr) {
    is_first = true;
    bound_distortion_K = 250;
    barrer_coef = d_.mesh_measure * 1e-8;
    std::cout << barrer_coef << std::endl;
  };
  ~Parafun();

  void after_mesh_improve();
  void init();
  void handle_mintri();
  void init_area();
  void setvirtualtri();
  void Pre_calculate();

  double BPE(bool is_ip_convrate, bool is_slim_convrate);
  void CM(bool is_interp = false);
  void SLIM(bool is_interp = false);
  void Update_source_same_t();

  void fungrid(const Eigen::VectorXd &x);

  void Energy(const Eigen::VectorXd &x, double &energy, bool is_interp = false,
              bool is_whole = true);
  void Energysource();
  double compute_energy(const Eigen::MatrixXd &x, bool whole = false);
  void adjust_shell_weight(double new_weight);

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

  ShellData &d_;
  int total_num;
  int F_N;
  int V_N;
  std::size_t dim = 2;
  std::vector<int> F0;
  std::vector<int> F1;
  std::vector<int> F2;
  Eigen::VectorXd position_of_mesh;
  std::vector<double> area;
  std::vector<double> area_src;
  double area_threshold;

  bool is_first;
  double Intp_T_Min;
  double changetocm_flag;
  double bound_distortion_K;
  std::vector<double> source_p00, source_p01, source_p10, source_p11;
  std::vector<double> update_p00, update_p01, update_p10, update_p11;

  std::unique_ptr<Solver> solver;
  std::vector<int> pardiso_ia, pardiso_ja;
  std::vector<double> pardiso_a, pardiso_b;

  double energy_mesh;
  double energy_barrier;
  double energy_all;
  double energy_shell;

  std::vector<int> id_h00, id_h01, id_h02, id_h03, id_h04, id_h05;
  std::vector<int> id_h11, id_h12, id_h13, id_h14, id_h15;
  std::vector<int> id_h22, id_h23, id_h24, id_h25;
  std::vector<int> id_h33, id_h34, id_h35;
  std::vector<int> id_h44, id_h45;
  std::vector<int> id_h55;

  double barrer_coef;
  double threhold;
  double average_length;
  double time1 = 0;
  double time2 = 0;
  double time3 = 0;
  double density = 0;

  int BE_N;
  int V_F_N;
  int AV_F_N;
  int AV_F_N_H;
  std::vector<int> V_F0, V_F1, V_F2;
  std::vector<int> V_F0_H, V_F1_H, V_F2_H;
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
