// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Function to sort data (replacement for quick_sort)
// [[Rcpp::export]]
arma::vec sort_vector(arma::vec v) {
  return arma::sort(v);
}

// Center a matrix and calculate variances
// [[Rcpp::export]]
List center_matrix_rcpp(arma::mat& x) {
  int n = x.n_rows;
  int p = x.n_cols;

  arma::vec var_x(p);
  arma::mat x_centered = x; // Create a copy to avoid modifying input directly

  for (int j = 0; j < p; ++j) {
    double mean_j = arma::mean(x_centered.col(j));
    x_centered.col(j) -= mean_j;
    var_x(j) = arma::accu(x_centered.col(j) % x_centered.col(j)) / n;
  }

  return List::create(Named("x") = x_centered, Named("var_x") = var_x);
}

// Calculate t-statistics
// [[Rcpp::export]]
arma::mat get_t_statistics_rcpp(const arma::mat& x, const arma::mat& y,
                                const arma::vec& var_x, const arma::vec& var_y) {
  int n = x.n_rows;
  int p = x.n_cols;
  int q = y.n_cols;

  arma::mat t = (x.t() * y) / n;

  for (int k = 0; k < q; ++k) {
    for (int j = 0; j < p; ++j) {
      double theta = t(j, k) / var_x(j);
      double theta_var = var_y(k) / var_x(j) - theta * theta;
      if (theta_var < 0.0) theta_var = 0.0;
      t(j, k) = theta / sqrt(theta_var / n);
    }
  }

  return t;
}

// Aggregate by cumulative sum
// [[Rcpp::export]]
arma::vec aggregate_by_cum_sum_rcpp(const arma::vec& t) {
  int p = t.n_elem;
  arma::vec temp = -t % t;
  arma::vec sorted_temp = arma::sort(temp);

  arma::vec g(p);
  double s = 0.0;
  for (int j = 0; j < p; ++j) {
    s += sorted_temp(j);
    g(j) = -s;
  }

  return g;
}

// Aggregate marginals
// [[Rcpp::export]]
arma::mat aggregate_marginals_rcpp(const arma::mat& t) {
  int p = t.n_rows;
  int q = t.n_cols;

  arma::mat g(q, p);

  for (int k = 0; k < q; ++k) {
    g.row(k) = aggregate_by_cum_sum_rcpp(t.col(k)).t();
  }

  return g;
}

// Get p-value
// [[Rcpp::export]]
double get_p_value_rcpp(double obs, const arma::vec& null) {
  int n = null.n_elem;

  if (obs <= null(0))
    return 1.0;
  else if (obs >= null(n - 1))
    return 1.0 / n;

  for (int i = 0; i < n - 1; i++) {
    if (obs > null(i) && obs <= null(i + 1)) {
      double adj = (obs - null(i)) / (null(i + 1) - null(i));
      return 1.0 - ((double)i + adj) / n;
    }
  }

  return 0.0; // Should not reach here
}

// Simulate null distribution
// [[Rcpp::export]]
List simulate_null_rcpp(const arma::mat& x, const arma::vec& var_x, int q, int n_sim) {
  int n = x.n_rows;
  int p = x.n_cols;

  // Create containers for all simulations
  arma::cube all_t(p, q, n_sim);
  arma::cube all_g(q, p, n_sim);

  for (int sim = 0; sim < n_sim; ++sim) {
    // Generate random normal data
    arma::mat epsilon(n, q, arma::fill::randn);

    // Center epsilon and get variance
    arma::vec var_epsilon(q);
    for (int j = 0; j < q; ++j) {
      double mean_j = arma::mean(epsilon.col(j));
      epsilon.col(j) -= mean_j;
      var_epsilon(j) = arma::accu(epsilon.col(j) % epsilon.col(j)) / n;
    }

    // Calculate t-statistics
    arma::mat t = get_t_statistics_rcpp(x, epsilon, var_x, var_epsilon);
    all_t.slice(sim) = t;

    // Calculate g statistics
    arma::mat g = aggregate_marginals_rcpp(t);
    all_g.slice(sim) = g;
  }

  return List::create(
    Named("t") = all_t,
    Named("g") = all_g
  );
}

// Main function for effect detection
// [[Rcpp::export]]
NumericVector detect_effect_rcpp(arma::mat x, arma::vec y, arma::vec alpha, int num_sim) {
  int n = x.n_rows;
  int p = x.n_cols;
  int n_alpha = alpha.n_elem;

  // Center x and get variance
  List centered_x_result = center_matrix_rcpp(x);
  arma::vec var_x = as<arma::vec>(centered_x_result["var_x"]);
  arma::mat x_centered = as<arma::mat>(centered_x_result["x"]);

  // Center y and get variance
  double mean_y = arma::mean(y);
  arma::vec y_centered = y - mean_y;
  double var_y = arma::accu(y_centered % y_centered) / n;
  arma::vec var_y_vec(1);
  var_y_vec(0) = var_y;

  // Convert y to matrix for processing
  arma::mat y_mat = arma::mat(y_centered);

  // Get observed t-statistics
  arma::mat obs_t = get_t_statistics_rcpp(x_centered, y_mat, var_x, var_y_vec);

  // Get observed g statistics
  arma::mat obs_g = aggregate_marginals_rcpp(obs_t);

  // Simulate null distribution for comparison
  List null_sim = simulate_null_rcpp(x_centered, var_x, 1, num_sim);
  arma::cube null_g_cube = as<arma::cube>(null_sim["g"]);

  // Find minimum p-value for observed statistics
  double obs_g_min = 1.0;
  arma::vec obs_p_values(n_alpha);

  for (int k = 0; k < n_alpha; ++k) {
    // Extract and sort null distribution for this index
    arma::vec null_g_k(num_sim);
    for (int i = 0; i < num_sim; ++i) {
      null_g_k(i) = null_g_cube(0, k, i);
    }
    null_g_k = arma::sort(null_g_k);

    // Calculate p-value
    double p_value = get_p_value_rcpp(obs_g(0, k), null_g_k);
    obs_p_values(k) = p_value;
    if (p_value < obs_g_min) obs_g_min = p_value;
  }

  // Simulate for comparison
  List sim = simulate_null_rcpp(x_centered, var_x, 1, num_sim);
  arma::cube sim_g_cube = as<arma::cube>(sim["g"]);

  // Calculate minimum p-values across all simulations
  arma::vec sim_g_min(num_sim, arma::fill::ones);
  arma::mat sim_p_values(num_sim, n_alpha);

  for (int k = 0; k < n_alpha; ++k) {
    // Extract null distribution for this index
    arma::vec null_g_k(num_sim);
    for (int i = 0; i < num_sim; ++i) {
      null_g_k(i) = null_g_cube(0, k, i);
    }
    null_g_k = arma::sort(null_g_k);

    // Calculate p-values for each simulation
    for (int i = 0; i < num_sim; ++i) {
      double sim_g_value = sim_g_cube(0, k, i);
      double p_value = get_p_value_rcpp(sim_g_value, null_g_k);
      sim_p_values(i, k) = p_value;
      if (p_value < sim_g_min(i)) sim_g_min(i) = p_value;
    }
  }

  // Compare observed p-values to simulated p-values
  NumericVector extreme(n_alpha + 1);
  for (int k = 0; k < n_alpha; ++k) {
    int count = 0;
    for (int i = 0; i < num_sim; ++i) {
      if (sim_p_values(i, k) < obs_p_values(k)) ++count;
    }
    extreme[k] = (double)count / num_sim;
  }

  // Adaptive test
  int count = 0;
  for (int i = 0; i < num_sim; ++i) {
    if (sim_g_min(i) < obs_g_min) ++count;
  }
  extreme[n_alpha] = (double)count / num_sim;

  return extreme;
}
