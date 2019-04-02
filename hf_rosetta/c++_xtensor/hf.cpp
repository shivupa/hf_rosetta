#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xmanipulation.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

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
  auto S = xt::adapt(S_vec, {num_ao, num_ao});

  std::vector<double> T_vec;
  load_double(T_vec, "../../data/T.txt");
  auto T = xt::adapt(T_vec, {num_ao, num_ao});

  std::vector<double> V_vec;
  load_double(V_vec, "../../data/V.txt");
  auto V = xt::adapt(V_vec, {num_ao, num_ao});

  std::vector<double> eri_vec;
  load_double(eri_vec, "../../data/eri.txt");
  auto eri = xt::adapt(eri_vec);

  double convergence_DM;
  load_double_scalar(convergence_DM, "../../data/convergence_DM.txt");
  double convergence_E;
  load_double_scalar(convergence_DM, "../../data/convergence_E.txt");
  double E_nuc;
  load_double_scalar(E_nuc, "../../data/E_nuc.txt");

  xt::xarray<double> D = xt::zeros<double>({num_ao, num_ao});
  xt::xarray<double> D_last = xt::zeros<double>({num_ao, num_ao});
  // loop variables
  int iteration_num = 0;
  double E_total = 0.0;
  xt::xarray<double> E_elec;
  xt::xarray<double> E_elec_last;
  xt::xarray<double> iteration_E_diff;
  xt::xarray<double> iteration_rmsc_dm;
  bool converged = false;
  bool exceeded_iterations = false;

  auto eigen_S = xt::linalg::eigh(S);
  auto s = std::get<0>(eigen_S);
  auto L = std::get<1>(eigen_S);
  xt::xarray<double> X = xt::zeros<double>({num_ao, num_ao});
  for (size_t i = 0; i < s.size(); i++) {
    X(i, i) = 1.0 / std::sqrt(s(i));
  }
  X = xt::linalg::dot(L, xt::linalg::dot(X, xt::transpose(L)));
  xt::xarray<double> H = T + V;
  xt::xarray<double> G = xt::zeros<double>({num_ao, num_ao});
  xt::xarray<double> F = xt::zeros<double>({num_ao, num_ao});
  xt::xarray<double> F_prime = xt::zeros<double>({num_ao, num_ao});
  xt::xarray<double> C = xt::zeros<double>({num_ao, num_ao});
  while (!converged && !exceeded_iterations) {
    E_elec_last = E_elec;
    D_last = D;

    iteration_num++;
    // form G matrix
    G = xt::zeros<double>({num_ao, num_ao});
    for (size_t i = 0; i < num_ao; i++) {
      for (size_t j = 0; j < num_ao; j++) {
        for (size_t k = 0; k < num_ao; k++) {
          for (size_t l = 0; l < num_ao; l++) {
            G(i, j) += D(k, l) * ((2.0 * (eri(idx4(i, j, k, l)))) -
                                  (eri(idx4(i, k, j, l))));
          }
        }
      }
    }

    F = H + G;
    F_prime = xt::linalg::dot(X, xt::linalg::dot(F, X));

    auto F_prime_eigen = xt::linalg::eigh(F_prime);
    auto E_orbitals = std::get<0>(F_prime_eigen);
    auto C_prime = std::get<1>(F_prime_eigen);

    C = xt::linalg::dot(X, C_prime);
    D = xt::zeros<double>({num_ao, num_ao});
    
    for (size_t i = 0; i < num_ao; i++) {
      for (size_t j = 0; j < num_ao; j++) {
        for (size_t k = 0; k < num_elec_alpha; k++) {
          D(i, j) += C(i, k) * C(j, k);
        }
      }
    }
    E_elec = xt::sum(D * (H + F));
    iteration_E_diff = xt::abs(E_elec - E_elec_last);
    iteration_rmsc_dm = xt::sqrt(xt::sum(xt::pow((D - D_last), 2)));
    if (iteration_E_diff(0) < convergence_E &&
        iteration_rmsc_dm(0) < convergence_DM) {
      converged = true;
    }
    if (iteration_num == iteration_max) {
      exceeded_iterations = true;
    }
  }
  E_total = E_elec(0) + E_nuc;
  printf("%20.15f\n", E_total);
  return 0;
}

