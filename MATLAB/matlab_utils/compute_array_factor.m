function array_factor = compute_array_factor(weights_complex, element_patterns_stacked)
% COMPUTE_ARRAY_FACTOR  Coherent array factor over the full angular grid.
%
%   array_factor = COMPUTE_ARRAY_FACTOR(weights_complex, element_patterns_stacked)
%
%   Implements AF(theta, phi) = sum_n  w_n * E_n(theta, phi), the coherent
%   superposition of all weighted element patterns.
%
%   Inputs:
%       weights_complex          : complex weights, length N_elements.
%                                  Units: dimensionless (V/V).
%       element_patterns_stacked : complex element patterns,
%                                  size (N_elements, N_theta, N_phi). Units: V/m.
%
%   Outputs:
%       array_factor : complex array factor, (N_theta x N_phi). Units: V/m.
%
%   Ported from src/cost/cost_function.py:compute_array_factor.
%   Layout note: the stack keeps the Python axis order (element, theta, phi);
%   the sum is taken over dimension 1 (elements).
%   Part of: Antenna Array Pattern Optimization Tool.

n_theta = size(element_patterns_stacked, 2);
n_phi   = size(element_patterns_stacked, 3);

% Broadcast weights over the spatial grid, then sum over the element axis.
% [PORT] Python w[:,None,None] * stack, sum(axis=0) -> reshape weights to (N,1,1).
w = reshape(weights_complex(:), [], 1, 1);                    % (N x 1 x 1)
array_factor = reshape(sum(w .* element_patterns_stacked, 1), n_theta, n_phi);
end
