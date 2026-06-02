function weights_complex = x_to_weights(x)
% X_TO_WEIGHTS  Decode the 2N real optimization vector into N complex weights.
%
%   weights_complex = X_TO_WEIGHTS(x)
%
%   The encoding interleaves real and imaginary parts:
%       x = [Re(w_0), Im(w_0), Re(w_1), Im(w_1), ..., Re(w_{N-1}), Im(w_{N-1})]
%
%   Inputs:
%       x : real optimization variable vector, length 2N. Units: dimensionless (V/V).
%
%   Outputs:
%       weights_complex : complex weight column vector, N x 1. Units: dimensionless (V/V).
%
%   Ported from src/cost/cost_function.py:x_to_weights.
%   Part of: Antenna Array Pattern Optimization Tool.

% Real parts occupy odd positions, imaginary parts even positions.
% [PORT] Python x[0::2] (0-based) -> MATLAB x(1:2:end) (1-based).
weights_real = x(1:2:end);
weights_imag = x(2:2:end);

% Force column orientation so the result is independent of the input shape.
weights_complex = weights_real(:) + 1i * weights_imag(:);
end
