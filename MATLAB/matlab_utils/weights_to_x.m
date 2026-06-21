function x = weights_to_x(weights_complex)
% WEIGHTS_TO_X  Encode N complex weights into the 2N real optimization vector.
%
%   x = WEIGHTS_TO_X(weights_complex)
%
%   Interleaves real and imaginary parts so that element n occupies positions
%   2n-1 (real) and 2n (imaginary) in the returned vector (1-based):
%       x = [Re(w_1), Im(w_1), Re(w_2), Im(w_2), ...]
%
%   Inputs:
%       weights_complex : complex weight vector, length N. Units: dimensionless (V/V).
%
%   Outputs:
%       x : real optimization variable column vector, 2N x 1. Units: dimensionless (V/V).
%
%   Ported from src/cost/cost_function.py:weights_to_x.
%   Part of: Antenna Array Pattern Optimization Tool.

% Stack real and imaginary side-by-side then flatten interleaved.
% [PORT] Python np.column_stack([re, im]).ravel() (row-major) -> fill odd/even slots.
weights_complex = weights_complex(:);          % force column
n_elements      = numel(weights_complex);

x = zeros(2 * n_elements, 1);
x(1:2:end) = real(weights_complex);
x(2:2:end) = imag(weights_complex);
end
