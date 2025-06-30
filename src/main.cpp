#include "tensor.h"
#include "tensor_ops.h"
#include <symengine/functions.h>
#include <symengine/printers.h>
#include <symengine/subs.h>
#include <iostream>

using namespace SymEngine;

int main() {
    // Define coordinates: t (time) and x (space)
    auto t = symbol("t");
    auto x = symbol("x");
    std::vector<RCP<const Symbol>> coords = {t, x};

    int dim = static_cast<int>(coords.size());

    // Define the scale factor a(t) and its derivatives
    auto a = function_symbol("a", t);
    auto adot = function_symbol("adot", t);   // First derivative da/dt
    auto addot = function_symbol("addot", t); // Second derivative d^2a/dt^2

    // Define the FRW metric g_ij
    Tensor g({dim, dim}, {"i", "j"});
    g.set_element({0, 0}, integer(-1));       // g_tt = -1
    g.set_element({1, 1}, mul(a, a));        // g_xx = a(t)^2
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
    subs_map[diff(a, t, false)] = adot;                          // da/dt -> adot
    subs_map[diff(diff(a, t, false), t, false)] = addot;         // d²a/dt² -> addot

    // Substitute in Christoffel symbols and Riemann tensor manually
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

    // Print results for all components
    std::cout << "FRW Metric Tensor g_ij:\n";
    g.print();

    std::cout << "\nInverse Metric Tensor g^ij:\n";
    g_inv.print();

    std::cout << "\nChristoffel Symbols Γ^k_ij (all components):\n";
    for (int k = 0; k < dim; ++k) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                std::cout << "Γ^" << k << "_{" << i << j << "} = "
                          << SymEngine::str(*Gamma.get_element({k, i, j})) << "\n";
            }
        }
    }

    std::cout << "\nRiemann Tensor R^m_nij (all components):\n";
    for (int m = 0; m < dim; ++m) {
        for (int n = 0; n < dim; ++n) {
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    std::cout << "R^" << m << "_{" << n << i << j << "} = "
                              << SymEngine::str(*Riemann.get_element({m, n, i, j})) << "\n";
                }
            }
        }
    }

    Tensor Ricci = compute_ricci(Riemann);
    auto RicciScalar = compute_ricci_scalar(Ricci, g_inv);
    Tensor Einstein = compute_einstein_tensor(Ricci, RicciScalar, g);

    std::cout << "\nRicci Tensor R_mu_nu:\n";
    Ricci.print();

    std::cout << "\nRicci Scalar R:\n";
    std::cout << *RicciScalar << "\n";

    std::cout << "\nEinstein Tensor G_mu_nu:\n";
    Einstein.print();

    return 0;
}
