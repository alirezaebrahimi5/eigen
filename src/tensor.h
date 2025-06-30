#ifndef TENSOR_H
#define TENSOR_H

#include <symengine/expression.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/symengine_rcp.h>
#include <symengine/solve.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace SymEngine;

class Tensor {
private:
    std::vector<int> shape;                // Tensor shape
    std::vector<std::string> indices;      // Indices (e.g., "i", "j")
    std::vector<SymEngine::RCP<const SymEngine::Basic>> elements; 
    int get_flat_index(const std::vector<int>& indices) const;

public:
    Tensor(const std::vector<int>& shape, const std::vector<std::string>& indices);
    void set_element(const std::vector<int>& indices, const RCP<const Basic>& value);
    RCP<const Basic> get_element(const std::vector<int>& indices) const;
    void simplify();
    void print() const;

    const std::vector<int>& get_shape() const;

    RCP<const Basic> sum_over_index(const std::string& index);
    Tensor contract(const std::string& index1, const std::string& index2);
    void differentiate(const RCP<const Symbol>& symbol);
    void substitute(const std::string& var, const RCP<const Basic>& value); // Substitute a variable
    RCP<const Basic> solve_equation(const RCP<const Basic>& equation, const RCP<const Symbol>& variable); // Solve equation

    // Christoffel symbols tensor: 3D tensor of symbolic expressions
    std::shared_ptr<Tensor> christoffel; // shape should be (dim, dim, dim)

    // Method to set Christoffel symbols
    void set_christoffel(const Tensor& Gamma);

    // Covariant derivative function
    Tensor covariant_derivative(const std::string& diff_index);
};

#endif // TENSOR_H
