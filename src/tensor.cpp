#include "tensor.h"

// Constructor
Tensor::Tensor(const std::vector<int>& shape, const std::vector<std::string>& indices)
    : shape(shape), indices(indices) {
    int total_size = 1;
    for (int dim : shape)
        total_size *= dim;
    elements.resize(total_size, zero); // Initialize with zeros
}

// Helper to calculate flat index
int Tensor::get_flat_index(const std::vector<int>& indices) const {
    int flat_index = 0;
    int stride = 1;
    for (int i = indices.size() - 1; i >= 0; --i) {
        flat_index += indices[i] * stride;
        stride *= shape[i];
    }
    return flat_index;
}

// Set an element
void Tensor::set_element(const std::vector<int>& indices, const RCP<const Basic>& value) {
    int flat_index = get_flat_index(indices);
    elements[flat_index] = value;
}

// Get an element
RCP<const Basic> Tensor::get_element(const std::vector<int>& indices) const {
    int flat_index = get_flat_index(indices);
    return elements[flat_index];
}

// Simplify tensor elements
void Tensor::simplify() {
    for (auto& element : elements) {
        element = expand(element); // Expand expressions
    }
}

// Substitute a variable with a value in all elements
void Tensor::substitute(const std::string& var, const RCP<const Basic>& value) {
    for (auto& element : elements) {
        element = element->subs({{symbol(var), value}});
    }
}



// Print the tensor
void Tensor::print() const {
    int rows = shape[0];
    int cols = shape[1];
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << *get_element({i, j}) << " ";
        }
        std::cout << std::endl;
    }
}


// Summation over an index
RCP<const Basic> Tensor::sum_over_index(const std::string& index) {
    int index_position = -1;
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] == index) {
            index_position = i;
            break;
        }
    }
    if (index_position == -1) {
        throw std::runtime_error("Index not found in tensor.");
    }

    RCP<const Basic> result = zero;
    for (int i = 0; i < shape[index_position]; ++i) {
        std::vector<int> coords(indices.size(), 0);
        coords[index_position] = i;
        result = add(result, get_element(coords));
    }
    return result;
}

// Contract two indices
Tensor Tensor::contract(const std::string& index1, const std::string& index2) {
    // Ensure index1 and index2 exist
    auto pos1 = std::find(indices.begin(), indices.end(), index1) - indices.begin();
    auto pos2 = std::find(indices.begin(), indices.end(), index2) - indices.begin();
    if (pos1 >= indices.size() || pos2 >= indices.size() || pos1 == pos2) {
        throw std::invalid_argument("Invalid indices for contraction");
    }

    std::vector<int> new_shape = shape;
    new_shape.erase(new_shape.begin() + std::max(pos1, pos2));
    new_shape.erase(new_shape.begin() + std::min(pos1, pos2));

    Tensor result(new_shape, indices);
    for (int i = 0; i < shape[pos1]; ++i) {
        for (int j = 0; j < shape[pos2]; ++j) {
            std::vector<int> coords(shape.size(), 0);
            coords[pos1] = i;
            coords[pos2] = j;

            auto new_value = add(result.get_element({i}), get_element(coords));
            result.set_element({i}, new_value);
        }
    }

    return result;
}


// Differentiate each tensor element with respect to a symbol
void Tensor::differentiate(const RCP<const Symbol>& symbol) {
    for (auto& element : elements) {
        element = diff(element, symbol, true); // Provide the third argument for proper differentiation
    }
}


RCP<const Basic> Tensor::solve_equation(const RCP<const Basic>& equation, const RCP<const Symbol>& variable) {
    // Use SymEngine's solve functionality
    return SymEngine::solve(equation, variable);
}


void Tensor::set_christoffel(const Tensor& Gamma) {
    christoffel = std::make_shared<Tensor>(Gamma);
}

// Assume diff_index corresponds to one of the tensor indices or a coordinate variable
Tensor Tensor::covariant_derivative(const std::string& diff_index) {
    // Step 1: Create new shape with one extra dimension for the derivative index
    std::vector<int> new_shape = shape; // Use int if your Tensor::shape uses int
    new_shape.insert(new_shape.begin(), 1); // Insert derivative index dimension at front

    // Step 2: Create new indices vector with diff_index prepended
    std::vector<std::string> new_indices;
    new_indices.push_back(diff_index);
    new_indices.insert(new_indices.end(), indices.begin(), indices.end());

    // Step 3: Initialize result tensor with new shape and indices
    Tensor result(new_shape, new_indices);

    // Step 4: Compute total elements in original tensor
    int total_elements = 1;
    for (auto dim : shape) {
        total_elements *= dim;
    }

    // Step 5: Compute partial derivatives element-wise
    for (int i = 0; i < total_elements; ++i) {
        auto elem = elements[i];
        // Differentiate with respect to the symbol corresponding to diff_index
        auto d_elem = SymEngine::diff(elem, symbol(diff_index), true);
        // Store in result elements at position with offset for new index dimension
        // Since new dimension size is 1, the index mapping here is straightforward
        result.elements[i] = d_elem;
    }

    // Step 6: TODO - Add connection terms involving Christoffel symbols and contractions
    // This requires detailed index bookkeeping and tensor algebra operations.

    return result;
}

// Getter for shape
const std::vector<int>& Tensor::get_shape() const {
    return shape;
}