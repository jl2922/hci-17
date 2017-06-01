#include "davidson.h"

void Davidson::diagonalize(const std::vector<double>& initial_vector) {
  const double TOLERANCE = 2.0e-7;
  const std::size_t MAX_ITERATIONS = 5;  // Good choice for direct evaluation of the hamiltonian.

  if (n == 1) {
    lowest_eigenvalue = diagonal[0];
    lowest_eigenvector = std::vector<double>(1, 1.0);
    diagonalized = true;
    return;
  }

  const std::size_t iterations = std::min(n, MAX_ITERATIONS);
  double lowest_eigenvalue = 0.0;
  double lowest_eigenvalue_prev = 0.0;
  double residual_norm = 0.0;

  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(n, iterations);

  if (initial_vector.size() != n) {
    v(0, 0) = 1.0;  // Start from HF.
  } else {
    for (std::size_t i = 0; i < n; i++) v(i, 0) = initial_vector[i];
    v.col(0).normalize();
  }

  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(n, iterations);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(n);  // Lowest eigenvector so far.
  Eigen::VectorXd Hw = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd h_krylov = Eigen::MatrixXd::Zero(iterations, iterations);
  Eigen::MatrixXd h_overwrite;
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(iterations);
  std::size_t len_work = 3 * iterations - 1;
  Eigen::VectorXd work(len_work);
  bool converged = false;
  std::vector<double> tmp_v(n);

  // Get diagonal elements.
  Eigen::VectorXd diag_elems(n);
  for (std::size_t i = 0; i < n; i++) diag_elems[i] = diagonal[i];

  // First iteration.
  for (std::size_t i = 0; i < n; i++) tmp_v[i] = v(i, 0);
  const auto& tmp_Hv = apply_hamiltonian(tmp_v);
  for (std::size_t i = 0; i < n; i++) Hv(i, 0) = tmp_Hv[i];
  lowest_eigenvalue = v.col(0).dot(Hv.col(0));
  h_krylov(0, 0) = lowest_eigenvalue;
  w = v.col(0);
  Hw = Hv.col(0);
  if (verbose) printf("Davidson Iteration #1. Eigenvalue: %#.15g\n", lowest_eigenvalue);

  residual_norm = 1.0;  // So at least one iteration is done.
  std::size_t n_iter = std::min(n, iterations);
  int n_diagonalize = 1;  // For print.

  for (std::size_t it = 1; it < n_iter; it++) {
    // Compute residual.
    for (std::size_t j = 0; j < n; j++) {
      v(j, it) = (Hw(j, 0) - lowest_eigenvalue * w(j, 0)) / (lowest_eigenvalue - diag_elems(j));
      if (fabs(lowest_eigenvalue - diag_elems[j]) < 1.0e-8) v(j, it) = -1.0;
    }

    // If residual is small, converge.
    residual_norm = v.col(it).norm();
    if (residual_norm < 1.0e-6) converged = true;

    // Orthogonalize and normalize.
    for (std::size_t i = 0; i < it; i++) {
      double norm = v.col(it).dot(v.col(i));
      v.col(it) -= norm * v.col(i);
    }
    v.col(it).normalize();

    // Apply H once.
    for (std::size_t i = 0; i < n; i++) tmp_v[i] = v(i, it);
    const auto& tmp_Hv2 = apply_hamiltonian(tmp_v);
    for (std::size_t i = 0; i < n; i++) Hv(i, it) = tmp_Hv2[i];

    // Construct Krylow matrix and diagonalize.
    for (std::size_t i = 0; i <= it; i++) {
      h_krylov(i, it) = v.col(i).dot(Hv.col(it));
      h_krylov(it, i) = h_krylov(i, it);
    }

    len_work = 3 * it + 2;
    h_overwrite = h_krylov.leftCols(it + 1).topRows(it + 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(
        h_krylov.leftCols(it + 1).topRows(it + 1));
    const auto& eigenvalues = eigenSolver.eigenvalues();
    const auto& eigenvectors = eigenSolver.eigenvectors();
    lowest_eigenvalue = eigenvalues[0];
    std::size_t lowest_id = 0;
    for (std::size_t i = 1; i <= it; i++) {
      if (eigenvalues[i] < lowest_eigenvalue) {
        lowest_eigenvalue = eigenvalues[i];
        lowest_id = i;
      }
    }
    w = v.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
    Hw = Hv.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);

    if (it > 1 && fabs(lowest_eigenvalue - lowest_eigenvalue_prev) < TOLERANCE) {
      converged = true;
      break;
    } else {
      lowest_eigenvalue_prev = lowest_eigenvalue;
      n_diagonalize++;
      if (verbose)
        printf("Davidson Iteration #%d. Eigenvalue: %#.15g\n", n_diagonalize, lowest_eigenvalue);
    }

    if (converged) break;
  }

  this->lowest_eigenvalue = lowest_eigenvalue;
  lowest_eigenvector.resize(n);
  for (std::size_t i = 0; i < n; i++) lowest_eigenvector[i] = w(i);
  diagonalized = true;
}

void Davidson::diagonalize(
    const std::vector<std::vector<double>>& initial_vectors, const std::size_t n_states) {
  const double TOLERANCE = 2.0e-7;
  const std::size_t MAX_ITERATIONS = 5;  // Good choice for direct evaluation of the hamiltonian.

  const std::size_t iterations = std::min(n, MAX_ITERATIONS);

  // Orthonormalize initial vectors.
  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(n, n_states * iterations);
  for (std::size_t i = 0; i < n_states; i++) {
    for (std::size_t j = 0; j < n; j++) v(j, i) = initial_vectors[i][j];
    v.col(i).normalize();
    if (i > 0) {
      for (std::size_t j = 0; j < i; j++) {
        const double cosine = v.col(i).dot(v.col(j));
        v.col(i) -= cosine * v.col(j);
      }
      v.col(i).normalize();
    }
  }

  if (n == 1) {
    lowest_eigenvalues.assign(1, diagonal[0]);
    lowest_eigenvectors.resize(1);
    lowest_eigenvectors[0].assign(1, 1.0);
    diagonalized = true;
    return;
  }

  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(n, n_states * iterations);
  Eigen::MatrixXd w = Eigen::MatrixXd::Zero(n, n_states);  // Lowest eigenvector so far.
  Eigen::MatrixXd Hw = Eigen::MatrixXd::Zero(n, n_states);
  Eigen::MatrixXd h_krylov = Eigen::MatrixXd::Zero(n_states * iterations, n_states * iterations);

  // Get diagonal elements.
  Eigen::VectorXd diag_elems(n);
  for (std::size_t i = 0; i < n; i++) diag_elems[i] = diagonal[i];

  // First iteration.
  std::vector<double> tmp_v(n);
  for (std::size_t i = 0; i < n_states; i++) {
    for (std::size_t j = 0; j < n; j++) tmp_v[j] = v(j, i);
    const auto& tmp_Hv = apply_hamiltonian(tmp_v);
    for (std::size_t j = 0; j < n; j++) Hv(j, i) = tmp_Hv[j];
  }

  lowest_eigenvalues.resize(n_states);
  for (std::size_t i = 0; i < n_states; i++) {
    lowest_eigenvalues[i] = v.col(i).dot(Hv.col(i));
    h_krylov(i, i) = lowest_eigenvalues[i];
    for (std::size_t j = i + 1; j < n_states; j++) {
      h_krylov(i, j) = v.col(i).dot(Hv.col(j));
      h_krylov(j, i) = h_krylov(i, j);
    }
  }

  if (verbose) {
    printf("Davidson Iteration #1. Eigenvalues:\n");
    for (std::size_t i = 0; i < n_states; i++) printf("%#.15g, ", lowest_eigenvalues[i]);
    printf("\n");
  }

  w = v.leftCols(n_states);
  Hw = Hv.leftCols(n_states);
  Eigen::VectorXd residual_norm = Eigen::VectorXd::Ones(n_states);
  std::vector<double> lowest_eigenvalues_prev(n_states);
  std::size_t n_iter = std::min(n, n_states * iterations);
  int n_diagonalize = 1;  // For print.
  bool converged = false;

  for (std::size_t it = n_states; it < n_iter; it++) {
    // Compute residual.
    std::size_t i = it % n_states;
    for (std::size_t j = 0; j < n; j++) {
      v(j, it) =
          (Hw(j, i) - lowest_eigenvalues[i] * w(j, i)) / (lowest_eigenvalues[i] - diag_elems(j));
      if (fabs(lowest_eigenvalues[i] - diag_elems[j]) < 1.0e-8) v(j, it) = -1.0;
    }

    // If residual is small, converge.
    residual_norm[i] = v.col(it).squaredNorm();
    if (residual_norm.sum() < 1.0e-12) converged = true;

    // Orthogonalize and normalize.
    for (std::size_t i = 0; i < it; i++) {
      double norm = v.col(it).dot(v.col(i));
      v.col(it) -= norm * v.col(i);
    }
    v.col(it).normalize();

    // Apply H once.
    for (std::size_t i = 0; i < n; i++) tmp_v[i] = v(i, it);
    const auto& tmp_Hv2 = apply_hamiltonian(tmp_v);
    for (std::size_t i = 0; i < n; i++) Hv(i, it) = tmp_Hv2[i];

    // Construct Krylov matrix and diagonalize.
    for (std::size_t i = 0; i <= it; i++) {
      h_krylov(i, it) = v.col(i).dot(Hv.col(it));
      h_krylov(it, i) = h_krylov(i, it);
    }

    if ((it + 1) % n_states == 0) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(h_krylov.leftCols(it).topRows(it));
      const auto& eigenvalues = eigenSolver.eigenvalues();
      const auto& eigenvectors = eigenSolver.eigenvectors();
      std::vector<std::size_t> lowest_ids(it);
      for (std::size_t i = 0; i < it; i++) lowest_ids[i] = i;
      std::stable_sort(
          lowest_ids.begin(),
          lowest_ids.end(),
          [eigenvalues](const std::size_t a, const std::size_t b) -> bool {
            return fabs(eigenvalues[a]) > fabs(eigenvalues[b]);
          });
      for (std::size_t i = 0; i < n_states; i++) {
        std::size_t lowest_id = lowest_ids[i];
        lowest_eigenvalues[i] = eigenvalues[lowest_id];
        w.col(i) = v.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
        Hw.col(i) = Hv.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
      }

      if (it > 0) {
        converged = true;
        for (std::size_t i = 0; i < n_states; i++) {
          if (fabs(lowest_eigenvalues[i] - lowest_eigenvalues_prev[i]) > TOLERANCE) {
            converged = false;
            break;
          }
        }
        if (converged) break;
      }

      for (std::size_t i = 0; i < n_states; i++) lowest_eigenvalues_prev[i] = lowest_eigenvalues[i];
      n_diagonalize++;
      if (verbose) {
        printf("Davidson Iteration #%d. Eigenvalues:\n", n_diagonalize);
        for (std::size_t i = 0; i < n_states; i++) printf("%#.15g, ", lowest_eigenvalues[i]);
        printf("\n");
      }
    }
  }

  lowest_eigenvectors.resize(n_states);
  for (std::size_t i = 0; i < n_states; i++) {
    lowest_eigenvectors[i].resize(n);
    for (std::size_t j = 0; j < n; j++) {
      lowest_eigenvectors[i][j] = w(j, i);
    }
  }
  diagonalized = true;
}
