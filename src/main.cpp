#include "tensor.h"
#include "tensor_ops.h"
#include <symengine/functions.h>
#include <symengine/printers.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <symengine/parser.h>

using namespace SymEngine;

void print_tensor_to_file(const Tensor& tensor, const std::string& name, std::ofstream& output_file) {
    output_file << "\n" << name << ":\n";
    const auto& shape = tensor.get_shape();
    size_t rank = shape.size();

    if (rank == 2) {
        // Handle metric tensors (rank-2 tensors)
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                auto elem = tensor.get_element({i, j});
                output_file << name << "_{" << i << " " << j << "} = " << SymEngine::str(*elem) << "\n";
            }
        }
    } 
    else if (rank == 3) {
        // Handle 3D tensors
        for (int k = 0; k < shape[0]; ++k) {
            for (int i = 0; i < shape[1]; ++i) {
                for (int j = 0; j < shape[2]; ++j) {
                    auto elem = tensor.get_element({k, i, j});
                    output_file << name << "^" << k << "_{" << i << " " << j << "} = " << SymEngine::str(*elem) << "\n";
                }
            }
        }
    } 
    else if (rank == 4) {
        // Handle 4D tensors
        for (int m = 0; m < shape[0]; ++m) {
            for (int n = 0; n < shape[1]; ++n) {
                for (int i = 0; i < shape[2]; ++i) {
                    for (int j = 0; j < shape[3]; ++j) {
                        auto elem = tensor.get_element({m, n, i, j});
                        output_file << name << "^" << m << "_{" << n << " " << i << " " << j << "} = " 
                                    << SymEngine::str(*elem) << "\n";
                    }
                }
            }
        }
    } 
    else {
        // Fallback for unsupported ranks
        output_file << "Tensor rank not supported for printing.\n";
    }
}


void parse_input_file(
    const std::string& file_path,
    int& dim,
    std::vector<RCP<const Symbol>>& coords,
    std::vector<RCP<const Basic>>& functions,
    Tensor& metric,
    Tensor& stress_energy_tensor, // Add a parameter for T_mu_nu
    RCP<const Basic>& c_val,
    RCP<const Basic>& G_val,
    RCP<const Basic>& pi_val,
    RCP<const Basic>& Lambda_val
) {
    std::ifstream input_file(file_path);
    if (!input_file) {
        throw std::runtime_error("Error opening input file.");
    }

    std::string line;
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            if (key == "dim") {
                dim = std::stoi(value);

            } else if (key == "variables") {
                std::istringstream vars(value);
                std::string var;
                while (std::getline(vars, var, ',')) {
                    coords.push_back(symbol(var));
                }

            } else if (key == "functions") {
                std::istringstream funcs(value);
                std::string func;
                while (std::getline(funcs, func, ',')) {
                    size_t open_paren = func.find('(');
                    size_t close_paren = func.find(')');
                    std::string name = func.substr(0, open_paren);
                    std::string arg = func.substr(open_paren + 1, close_paren - open_paren - 1);
                    functions.push_back(function_symbol(name, symbol(arg)));
                }

            } else if (key == "metric") {
                metric = Tensor({dim, dim}, {"i", "j"});
                std::istringstream elements(value);
                std::string elem;
                while (std::getline(elements, elem, ';')) {
                    size_t colon = elem.find(':');
                    std::string indices = elem.substr(0, colon);
                    std::string val = elem.substr(colon + 1);

                    int i = indices[0] - '0';
                    int j = indices[2] - '0';

                    metric.set_element({i, j}, parse(val));
                }

            } else if (key == "T") { // Add a new key for T_mu_nu
                stress_energy_tensor = Tensor({dim, dim}, {"mu", "nu"});
                std::istringstream elements(value);
                std::string elem;
                while (std::getline(elements, elem, ';')) {
                    size_t colon = elem.find(':');
                    std::string indices = elem.substr(0, colon);
                    std::string val = elem.substr(colon + 1);

                    int mu = indices[0] - '0';
                    int nu = indices[2] - '0';

                    stress_energy_tensor.set_element({mu, nu}, parse(val));
                }

            } else if (key == "c") {
                c_val = parse(value);

            } else if (key == "G") {
                G_val = parse(value);

            } else if (key == "pi") {
                pi_val = parse(value);

            } else if (key == "Lambda") {
                Lambda_val = parse(value);
            }
        }
    }
}



int main() {
    int dim;
    std::vector<RCP<const Symbol>> coords;
    std::vector<RCP<const Basic>> functions;
    Tensor metric({4, 4}, {"i", "j"});

    Tensor T({4, 4}, {"mu", "nu"}); 
    RCP<const Basic> c_val;
    RCP<const Basic> G_val;
    RCP<const Basic> pi_val;
    RCP<const Basic> Lambda_val;

    try {
        parse_input_file("input.txt", dim, coords, functions, metric, T, c_val, G_val, pi_val, Lambda_val);

        Tensor g_inv = compute_inverse_metric(metric);
        Tensor Gamma = compute_christoffel(metric, g_inv, coords);
        Tensor Riemann = compute_riemann(Gamma, coords);

        int limit = 2;
        for (int m = 0; m < dim; ++m)
            for (int n = 0; n < dim; ++n)
                for (int i = 0; i < dim; ++i)
                    for (int j = 0; j < dim; ++j) {
                        auto elem = Riemann.get_element({m, n, i, j});
                        Riemann.set_element({m, n, i, j}, enforce_limit(elem, coords[0], limit));
                    }

        Tensor Ricci = compute_ricci(Riemann);
        auto RicciScalar = compute_ricci_scalar(Ricci, g_inv);
        Tensor Einstein = compute_einstein_tensor(Ricci, RicciScalar, metric);

        std::ofstream output_file("output.txt");
        if (!output_file) {
            throw std::runtime_error("Error opening output file.");
        }


        Tensor EFE_left_side({dim, dim}, {"mu", "nu"});
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                auto G_mn = Einstein.get_element({i, j});
                auto g_mn = metric.get_element({i, j});
                // EFE left side: G_{mu nu} + Lambda * g_{mu nu}
                EFE_left_side.set_element({i, j}, add(G_mn, mul(Lambda_val, g_mn)));
            }
        }

        // Compute RHS: 8 * pi * G * T
        Tensor EFE_right_side({dim, dim}, {"mu", "nu"});
        auto eight_pi_G = mul(parse("8"), mul(pi_val, G_val));
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                auto T_mn = T.get_element({i, j});
                EFE_right_side.set_element({i, j}, mul(eight_pi_G, T_mn));
            }
        }

        // Compute the difference between LHS and RHS
        Tensor EFE_difference({dim, dim}, {"mu", "nu"});
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                auto lhs = EFE_left_side.get_element({i, j});
                auto rhs = EFE_right_side.get_element({i, j});
                EFE_difference.set_element({i, j}, sub(lhs, rhs));
            }
        }

        print_tensor_to_file(metric, "FRW Metric Tensor g_ij", output_file);
        print_tensor_to_file(g_inv, "Inverse Metric Tensor g^ij", output_file);
        print_tensor_to_file(Gamma, "Christoffel Symbols Γ^k_ij", output_file);
        print_tensor_to_file(Riemann, "Riemann Tensor R^m_nij", output_file);
        print_tensor_to_file(Ricci, "Ricci Tensor R_mu_nu", output_file);
        print_tensor_to_file(Einstein, "Einstein Tensor G_mu_nu", output_file);
        print_tensor_to_file(EFE_left_side, "EFE Left Side (G_mu_nu + Lambda * g_mu_nu)", output_file);
        print_tensor_to_file(EFE_right_side, "EFE Right Side (8 pi G T_mu_nu)", output_file);
        print_tensor_to_file(EFE_difference, "EFE Difference (LHS - RHS)", output_file);

        output_file.close();
        std::cout << "Calculation completed. Results saved to output.txt.\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
