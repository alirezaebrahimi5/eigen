#include "tensor.h"
#include "tensor_ops.h"
#include <symengine/functions.h>
#include <symengine/printers.h>
#include <iostream>

using namespace SymEngine;

void print_tensor(const Tensor& tensor, const std::string& name) {
    std::cout << "\n" << name << ":\n";
    const auto& shape = tensor.get_shape();

    if (shape.size() == 2) {
        // For rank 2 tensors, print with 4 repeated indices for formatting:
        // e.g. G_mu_nu_{ijij}
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                auto elem = tensor.get_element({i, j});
                std::cout << name << "_{" << i << j << " " << i << j << "} = "
                          << SymEngine::str(*elem) << "\n";
            }
        }
    }
    else if (shape.size() == 4) {
        // For rank 4 tensors, print normally
        for (int m = 0; m < shape[0]; ++m)
            for (int n = 0; n < shape[1]; ++n)
                for (int i = 0; i < shape[2]; ++i)
                    for (int j = 0; j < shape[3]; ++j)
                        std::cout << name << "^" << m << "_{" << n << " " << i << " " << j << "} = "
                                  << SymEngine::str(*tensor.get_element({m, n, i, j})) << "\n";
    }
    else {
        std::cout << "Tensor rank not supported for printing.\n";
    }
}

int main() {
    // Define 4D coordinates
    auto t = symbol("t");
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    std::vector<RCP<const Symbol>> coords = {t, x, y, z};

    int dim = 4;

    // Define scale factors a(t), theta(t)
    auto a = function_symbol("a", t);
    auto theta = function_symbol("theta", t);

    // Initialize metric tensor g_ij
    Tensor g({dim, dim}, {"i", "j"});
    g.set_element({0, 0}, integer(-1));  // g_tt = -1
    for (int i = 1; i < dim; ++i) {
        g.set_element({i, i}, add(mul(a, a), mul(theta, theta)));  // spatial parts
    }
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            if (i != j) {
                g.set_element({i, j}, zero);
            }
        }
    }

    // Compute inverse metric
    Tensor g_inv = compute_inverse_metric(g);

    // Compute Christoffel symbols Γ^k_ij
    Tensor Gamma = compute_christoffel(g, g_inv, coords);

    // Compute Riemann curvature tensor R^m_nij
    Tensor Riemann = compute_riemann(Gamma, coords);

    // Limit powers if needed (limit=2)
    int limit = 2;
    for (int m = 0; m < dim; ++m)
        for (int n = 0; n < dim; ++n)
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j) {
                    auto elem = Riemann.get_element({m, n, i, j});
                    Riemann.set_element({m, n, i, j}, enforce_limit(elem, t, limit));
                }

    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) {
            auto elem = g.get_element({i, j});
            g.set_element({i, j}, enforce_limit(elem, t, limit));
        }

    // Print tensors
    print_tensor(g, "FRW Metric Tensor g_ij");
    print_tensor(g_inv, "Inverse Metric Tensor g^ij");
    print_tensor(Gamma, "Christoffel Symbols Γ^k_ij");
    print_tensor(Riemann, "Riemann Tensor R^m_nij");

    // Compute Ricci, Ricci scalar, Einstein tensors
    Tensor Ricci = compute_ricci(Riemann);
    auto RicciScalar = compute_ricci_scalar(Ricci, g_inv);
    Tensor Einstein = compute_einstein_tensor(Ricci, RicciScalar, g);

    // Limit powers for Ricci and Einstein
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) {
            Ricci.set_element({i, j}, enforce_limit(Ricci.get_element({i, j}), t, limit));
            Einstein.set_element({i, j}, enforce_limit(Einstein.get_element({i, j}), t, limit));
        }

    print_tensor(Ricci, "Ricci Tensor R_mu_nu");
    print_tensor(Einstein, "Einstein Tensor G_mu_nu");

    return 0;
}