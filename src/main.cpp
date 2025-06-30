#include "tensor.h"
#include "tensor_ops.h"
#include <symengine/functions.h>
#include <symengine/printers.h>
#include <iostream>

using namespace SymEngine;

// Helper function to print tensors
void print_tensor(const Tensor& tensor, const std::string& name) {
    std::cout << "\n" << name << ":\n";
    const auto& shape = tensor.get_shape();

    if (shape.size() == 2) {
        // 2D tensor (Ricci, Einstein)
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                std::cout << name << "_{" << i << j << "} = "
                          << SymEngine::str(*tensor.get_element({i, j})) << "\n";
            }
        }
    } else if (shape.size() == 3) {
        // 3D tensor (Christoffel)
        for (int k = 0; k < shape[0]; ++k) {
            for (int i = 0; i < shape[1]; ++i) {
                for (int j = 0; j < shape[2]; ++j) {
                    std::cout << "Γ^" << k << "_{" << i << j << "} = "
                              << SymEngine::str(*tensor.get_element({k, i, j})) << "\n";
                }
            }
        }
    } else if (shape.size() == 4) {
        // 4D tensor (Riemann)
        for (int m = 0; m < shape[0]; ++m) {
            for (int n = 0; n < shape[1]; ++n) {
                for (int i = 0; i < shape[2]; ++i) {
                    for (int j = 0; j < shape[3]; ++j) {
                        std::cout << "R^" << m << "_{" << n << i << j << "} = "
                                  << SymEngine::str(*tensor.get_element({m, n, i, j})) << "\n";
                    }
                }
            }
        }
    } else {
        std::cout << "Tensor rank not supported for printing.\n";
    }
}


int main() {
    // Define coordinates: t (time) and x (space)
    auto t = symbol("t");
    auto x = symbol("x");
    std::vector<RCP<const Symbol>> coords = {t, x};

    int dim = static_cast<int>(coords.size());

    // Define the scale factor a(t) and its derivatives
    auto a = function_symbol("a", t);
    auto adot = function_symbol("adot", t);   // da/dt
    auto addot = function_symbol("addot", t); // d²a/dt²

    // Define the new scalar function Theta(t) and derivatives
    auto theta = function_symbol("theta", t);
    auto thetadot = function_symbol("thetadot", t);     // dTheta/dt
    auto thetaddot = function_symbol("thetaddot", t);   // d²Theta/dt²

    // Define the FRW metric g_ij with added Theta(t)^2 term
    Tensor g({dim, dim}, {"i", "j"});
    g.set_element({0, 0}, integer(-1));       // g_tt = -1
    g.set_element({1, 1}, add(mul(a, a), mul(theta, theta)));  // g_xx = a(t)^2 + theta(t)^2
    g.set_element({0, 1}, zero);             // g_tx = 0
    g.set_element({1, 0}, zero);             // g_xt = 0

    // Compute inverse metric g^ij
    Tensor g_inv = compute_inverse_metric(g);

    // Compute Christoffel symbols Γ^k_ij
    Tensor Gamma = compute_christoffel(g, g_inv, coords);

    // Compute Riemann curvature tensor R^m_nij
    Tensor Riemann = compute_riemann(Gamma, coords);

    // Create substitution map
    map_basic_basic subs_map;
    subs_map[diff(a, t, false)] = adot;                            // da/dt -> adot
    subs_map[diff(diff(a, t, false), t, false)] = addot;           // d²a/dt² -> addot

    subs_map[diff(theta, t, false)] = thetadot;                    // dTheta/dt -> thetadot
    subs_map[diff(diff(theta, t, false), t, false)] = thetaddot;   // d²Theta/dt² -> thetaddot

    // Substitute derivatives in Christoffel symbols and Riemann tensor
    for (int k = 0; k < dim; ++k) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Gamma.set_element({k, i, j}, Gamma.get_element({k, i, j})->subs(subs_map));
            }
        }
    }

    for (int m = 0; m < dim; ++m) {
        for (int n = 0; n < dim; ++n) {
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    Riemann.set_element({m, n, i, j}, Riemann.get_element({m, n, i, j})->subs(subs_map));
                }
            }
        }
    }

    // Compute Ricci tensor, scalar, and Einstein tensor
    Tensor Ricci = compute_ricci(Riemann);
    auto RicciScalar = compute_ricci_scalar(Ricci, g_inv);
    Tensor Einstein = compute_einstein_tensor(Ricci, RicciScalar, g);

    // Print tensors
    print_tensor(g, "FRW Metric Tensor g_ij");
    print_tensor(g_inv, "Inverse Metric Tensor g^ij");
    print_tensor(Gamma, "Christoffel Symbols Γ^k_ij");
    print_tensor(Riemann, "Riemann Tensor R^m_nij");
    print_tensor(Ricci, "Ricci Tensor R_mu_nu");
    print_tensor(Einstein, "Einstein Tensor G_mu_nu");

    return 0;
}