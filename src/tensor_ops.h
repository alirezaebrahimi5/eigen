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

// Add declarations for Ricci, Ricci scalar, and Einstein tensors:
Tensor compute_ricci(const Tensor& Riemann);
RCP<const Basic> compute_ricci_scalar(const Tensor& Ricci, const Tensor& g_inv);
Tensor compute_einstein_tensor(const Tensor& Ricci, const RCP<const Basic>& RicciScalar, const Tensor& g);


#endif // TENSOR_OPS_H
