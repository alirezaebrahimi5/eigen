#include "tensor_ops.h"
#include <symengine/matrix.h>
#include <symengine/printers.h>
#include <stdexcept>


// Perform Gauss-Jordan elimination to invert a square matrix `mat`.
// Throws if matrix is singular.
DenseMatrix invert_matrix(const DenseMatrix &mat) {
    size_t n = mat.nrows();
    if (n != mat.ncols()) {
        throw std::runtime_error("Matrix must be square for inversion");
    }

    DenseMatrix augmented(n, 2*n);

    // Initialize augmented matrix [mat | I]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented.set(i, j, mat.get(i,j));
        }
        for (size_t j = n; j < 2*n; ++j) {
            augmented.set(i, j, (i == (j-n)) ? one : zero);
        }
    }

    // Gauss-Jordan elimination
    for (size_t col = 0; col < n; ++col) {
        // Find pivot row
        size_t pivot = col;
        for (size_t row = col+1; row < n; ++row) {
            // Here just pick first non-zero pivot (symbolic comparison is complex)
            if (!eq(*augmented.get(row,col), *zero)) {
                pivot = row;
                break;
            }
        }
        if (eq(*augmented.get(pivot,col), *zero)) {
            throw std::runtime_error("Matrix is singular, cannot invert");
        }

        // Swap rows if needed
        if (pivot != col) {
            for (size_t j = 0; j < 2*n; ++j) {
                auto temp = augmented.get(col,j);
                augmented.set(col, j, augmented.get(pivot, j));
                augmented.set(pivot, j, temp);
            }
        }

        // Normalize pivot row
        auto pivot_val = augmented.get(col, col);
        for (size_t j = 0; j < 2*n; ++j) {
            augmented.set(col, j, div(augmented.get(col, j), pivot_val));
        }

        // Eliminate other rows
        for (size_t row = 0; row < n; ++row) {
            if (row == col)
                continue;
            auto factor = augmented.get(row, col);
            for (size_t j = 0; j < 2*n; ++j) {
                auto val = sub(augmented.get(row, j), mul(factor, augmented.get(col, j)));
                augmented.set(row, j, val);
            }
        }
    }

    // Extract inverse matrix (right half)
    DenseMatrix inv_mat(n, n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            inv_mat.set(i, j, augmented.get(i, j + n));

    return inv_mat;
}



// Compute inverse metric g_inv from metric g
Tensor compute_inverse_metric(const Tensor& g) {
    size_t dim = g.get_shape()[0];
    if (g.get_shape().size() != 2 || g.get_shape()[1] != dim) {
        throw std::runtime_error("Metric tensor must be square.");
    }

    SymEngine::DenseMatrix g_mat(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            g_mat.set(i, j, g.get_element({static_cast<int>(i), static_cast<int>(j)}));
        }
    }

    // Use your invert_matrix function
    SymEngine::DenseMatrix g_inv_mat = invert_matrix(g_mat);

    std::vector<int> shape = {static_cast<int>(dim), static_cast<int>(dim)};
    std::vector<std::string> indices = {"i", "j"};
    Tensor g_inv(shape, indices);

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            g_inv.set_element({static_cast<int>(i), static_cast<int>(j)}, g_inv_mat.get(i, j));
        }
    }

    return g_inv;
}



// Compute Christoffel symbols Γ^k_{ij}
Tensor compute_christoffel(const Tensor& g, const Tensor& g_inv, const std::vector<RCP<const Symbol>>& coords) {
    size_t dim = g.get_shape()[0];
    std::vector<int> shape = {static_cast<int>(dim), static_cast<int>(dim), static_cast<int>(dim)};
    std::vector<std::string> indices = {"k", "i", "j"};
    Tensor Gamma(shape, indices);

    for (size_t k = 0; k < dim; ++k) {
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                RCP<const Basic> sum = zero;
                for (size_t l = 0; l < dim; ++l) {
                    auto dg_jl_di = diff(g.get_element({(int)j, (int)l}), coords[i], true);
                    auto dg_il_dj = diff(g.get_element({(int)i, (int)l}), coords[j], true);
                    auto dg_ij_dl = diff(g.get_element({(int)i, (int)j}), coords[l], true);

                    auto term = mul(g_inv.get_element({(int)k, (int)l}),
                                    add(add(dg_jl_di, dg_il_dj), mul(integer(-1), dg_ij_dl)));

                    sum = add(sum, term);
                }
                sum = mul(div(integer(1), integer(2)), sum);
                Gamma.set_element({(int)k, (int)i, (int)j}, sum);
            }
        }
    }

    return Gamma;
}

// Compute Riemann curvature tensor R^m_{nij}
Tensor compute_riemann(const Tensor& Gamma, const std::vector<RCP<const Symbol>>& coords) {
    size_t dim = Gamma.get_shape()[0];
    std::vector<int> shape = {static_cast<int>(dim), static_cast<int>(dim), static_cast<int>(dim), static_cast<int>(dim)};
    std::vector<std::string> indices = {"m", "n", "i", "j"};
    Tensor R(shape, indices);

    for (size_t m = 0; m < dim; ++m) {
        for (size_t n = 0; n < dim; ++n) {
            for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    auto d_Gamma_mjn_di = diff(Gamma.get_element({(int)m, (int)j, (int)n}), coords[i], true);
                    auto d_Gamma_min_dj = diff(Gamma.get_element({(int)m, (int)i, (int)n}), coords[j], true);

                    RCP<const Basic> sum = zero;
                    for (size_t p = 0; p < dim; ++p) {
                        auto term1 = mul(Gamma.get_element({(int)m, (int)i, (int)p}),
                                         Gamma.get_element({(int)p, (int)j, (int)n}));
                        auto term2 = mul(Gamma.get_element({(int)m, (int)j, (int)p}),
                                         Gamma.get_element({(int)p, (int)i, (int)n}));
                        sum = add(sum, add(term1, mul(integer(-1), term2)));
                    }
                    RCP<const Basic> val = add(sub(d_Gamma_mjn_di, d_Gamma_min_dj), sum);
                    R.set_element({(int)m, (int)n, (int)i, (int)j}, val);
                }
            }
        }
    }
    return R;
}

// Contract Riemann tensor to get Ricci tensor: R_{mu,nu} = R^{alpha}_{mu,alpha,nu}
Tensor compute_ricci(const Tensor& Riemann) {
    auto shape = Riemann.get_shape(); // Expect {dim, dim, dim, dim}
    int dim = shape[0];

    Tensor Ricci({dim, dim}, {"mu", "nu"});

    for (int mu = 0; mu < dim; ++mu) {
        for (int nu = 0; nu < dim; ++nu) {
            RCP<const Basic> sum = zero;
            for (int alpha = 0; alpha < dim; ++alpha) {
                sum = add(sum, Riemann.get_element({alpha, mu, alpha, nu}));
            }
            Ricci.set_element({mu, nu}, sum);
        }
    }
    return Ricci;
}

// Compute Ricci scalar: R = g^{mu,nu} R_{mu,nu}
RCP<const Basic> compute_ricci_scalar(const Tensor& Ricci, const Tensor& g_inv) {
    auto shape = Ricci.get_shape();
    int dim = shape[0];

    RCP<const Basic> R = zero;

    for (int mu = 0; mu < dim; ++mu) {
        for (int nu = 0; nu < dim; ++nu) {
            auto term = mul(g_inv.get_element({mu, nu}), Ricci.get_element({mu, nu}));
            R = add(R, term);
        }
    }
    return R;
}

// Compute Einstein tensor: G_{mu,nu} = R_{mu,nu} - 1/2 * R * g_{mu,nu}
Tensor compute_einstein_tensor(const Tensor& Ricci, const RCP<const Basic>& RicciScalar, const Tensor& g) {
    auto shape = Ricci.get_shape();
    int dim = shape[0];

    Tensor Einstein({dim, dim}, {"mu", "nu"});

    auto half = div(integer(1), integer(2));

    for (int mu = 0; mu < dim; ++mu) {
        for (int nu = 0; nu < dim; ++nu) {
            auto val = sub(Ricci.get_element({mu, nu}), mul(half, mul(RicciScalar, g.get_element({mu, nu}))));
            Einstein.set_element({mu, nu}, val);
        }
    }
    return Einstein;
}
