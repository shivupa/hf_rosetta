#include "Eigen/Dense"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

typedef Eigen::Matrix< double , Eigen::Dynamic , Eigen::Dynamic > MatrixXd;

void load_int_scalar(int &in, std::string filename) {
  std::ifstream input_file(filename);
  input_file >> in;
}

void load_double_scalar(double &in, std::string filename) {
  std::ifstream input_file(filename);
  input_file >> in;
}

void load_double(std::vector<double> &in, std::string filename) {
  std::ifstream input_file(filename);
  double temp_val;
  std::string line;
  while (std::getline(input_file, line)) {
    std::istringstream ss(line);
    while (ss >> temp_val) {
      in.push_back(temp_val);
    }
  }
}

size_t idx2(size_t &i, size_t &j) {
  if (i >= j) {
    return i * (i + 1) / 2 + j;
  } else {
    return j * (j + 1) / 2 + i;
  }
}

size_t idx4(size_t &i, size_t &j, size_t &k, size_t &l) {
  size_t temp_ij = idx2(i, j);
  size_t temp_kl = idx2(k, l);
  return idx2(temp_ij, temp_kl);
}

int main(int argc, char const *argv[]) {
  int iteration_max;
  int num_elec_alpha;
  int num_elec_beta;
  int num_ao;

  load_int_scalar(num_elec_alpha, "../../data/num_elec_alpha.txt");
  load_int_scalar(num_elec_beta, "../../data/num_elec_beta.txt");
  load_int_scalar(iteration_max, "../../data/iteration_max.txt");
  load_int_scalar(num_ao, "../../data/num_ao.txt");

  std::vector<double> S_vec;
  load_double(S_vec, "../../data/S.txt");
  Eigen::Map<MatrixXd> S(S_vec.data(), num_ao, num_ao);

  std::vector<double> T_vec;
  load_double(T_vec, "../../data/T.txt");
  Eigen::Map<MatrixXd> T(T_vec.data(), num_ao, num_ao);

  std::vector<double> V_vec;
  load_double(V_vec, "../../data/V.txt");
  Eigen::Map<MatrixXd> V(V_vec.data(), num_ao, num_ao);

  std::vector<double> eri_vec;
  load_double(eri_vec, "../../data/eri.txt");
  Eigen::Map<MatrixXd> eri(eri_vec.data(),45150,1);

  double convergence_DM;
  load_double_scalar(convergence_DM, "../../data/convergence_DM.txt");
  double convergence_E;
  load_double_scalar(convergence_DM, "../../data/convergence_E.txt");
  double E_nuc;
  load_double_scalar(E_nuc, "../../data/E_nuc.txt");
  
  MatrixXd D = MatrixXd::Zero(num_ao,num_ao);
  MatrixXd D_last = MatrixXd::Zero(num_ao,num_ao);
  // loop variables
  int iteration_num = 0;
  double E_total = 0.0;
  double E_elec = 0.0;
  double E_elec_last = 0.0;
  double iteration_E_diff = 0.0;
  double iteration_rmsc_dm = 0.0;
  bool converged = false;
  bool exceeded_iterations = false;
  Eigen::SelfAdjointEigenSolver<MatrixXd> eigen_solver(num_ao);
  eigen_solver.compute(S);
  auto s = eigen_solver.eigenvalues().transpose();
  auto L = eigen_solver.eigenvectors();
  MatrixXd X = MatrixXd::Zero(num_ao,num_ao);
  for (size_t i = 0; i < s.size(); i++) {
    X(i, i) = 1.0 / std::sqrt(s(i));
  }
  X = L*X*L.transpose();
  auto H = T + V;
  MatrixXd G = MatrixXd::Zero(num_ao,num_ao);
  MatrixXd F = MatrixXd::Zero(num_ao,num_ao);
  MatrixXd F_prime = MatrixXd::Zero(num_ao,num_ao);
  MatrixXd C = MatrixXd::Zero(num_ao,num_ao);
  while (!converged && !exceeded_iterations) {
    E_elec_last = E_elec;
    D_last = D;

    iteration_num++;
    // form G matrix
    MatrixXd G = MatrixXd::Zero(num_ao,num_ao);
    for (size_t i = 0; i < num_ao; i++) {
      for (size_t j = 0; j < num_ao; j++) {
        for (size_t k = 0; k < num_ao; k++) {
          for (size_t l = 0; l < num_ao; l++) {
            G(i, j) += D(k, l) * ((2.0 * (eri(idx4(i, j, k, l)))) - (eri(idx4(i, k, j, l))));
          }
        }
      }
    }

    F = H + G;
    F_prime = X*F*X;
    
    eigen_solver.compute(F_prime);
    auto E_orbitals = eigen_solver.eigenvalues().transpose();
    auto C_prime = eigen_solver.eigenvectors();

    C = X*C_prime;
    D.setZero();
    
    for (size_t i = 0; i < num_ao; i++) {
      for (size_t j = 0; j < num_ao; j++) {
        for (size_t k = 0; k < num_elec_alpha; k++) {
          D(i, j) += C(i, k) * C(j, k);
        }
      }
    }
    E_elec = D.cwiseProduct(H + F).sum();
    iteration_E_diff = std::abs(E_elec - E_elec_last);
    iteration_rmsc_dm = std::sqrt((D-D_last).cwiseProduct(D-D_last).sum());
    if (iteration_E_diff < convergence_E &&
        iteration_rmsc_dm < convergence_DM) {
      converged = true;
    }
    if (iteration_num == iteration_max) {
      exceeded_iterations = true;
    }
  }
  E_total = E_elec + E_nuc;
  printf("%20.15f\n", E_total);
  return 0;
}
