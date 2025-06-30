#ifndef TENSOR_OPS_H
#define TENSOR_OPS_H

#include "tensor.h"
#include <symengine/symengine_rcp.h>
#include <vector>

using namespace SymEngine;

// Declare functions:
Tensor compute_inverse_metric(const Tensor& g);
Tensor compute_christoffel(const Tensor& g, const Tensor& g_inv, const std::vector<RCP<const Symbol>>& coords);
Tensor compute_riemann(const Tensor& Gamma, const std::vector<RCP<const Symbol>>& coords);

#endif // TENSOR_OPS_H
