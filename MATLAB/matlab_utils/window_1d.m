function w = window_1d(name, n_side)
% WINDOW_1D  1-D aperture taper window matching numpy / scipy definitions.
%
%   w = WINDOW_1D(name, n_side)
%
%   Base-MATLAB reimplementation (no Signal Processing Toolbox) of the windows
%   used by compare_classical.py. Each matches its numpy/scipy counterpart:
%       uniform        -> ones
%       hamming        -> numpy.hamming
%       hanning        -> numpy.hanning
%       kaiser_b3      -> numpy.kaiser(n, 3.0)
%       kaiser_b6      -> numpy.kaiser(n, 6.0)
%       chebyshev_25dB -> scipy.signal.windows.chebwin(n, at=25.0)
%       chebyshev_40dB -> scipy.signal.windows.chebwin(n, at=40.0)
%       taylor_25dB    -> scipy.signal.windows.taylor(n, nbar=4, sll=25, norm=True)
%
%   Inputs:
%       name   : window name (char) as above.
%       n_side : number of samples (= array elements per dimension).
%
%   Outputs:
%       w : real window values, (n_side x 1). Units: dimensionless.
%
%   Ported from scripts/compare_classical.py:_window_1d (+ scipy/numpy sources).
%   Part of: Antenna Array Pattern Optimization Tool.

switch name
    case 'uniform'
        w = ones(n_side, 1);
    case 'hamming'
        w = hamming_win(n_side);
    case 'hanning'
        w = hanning_win(n_side);
    case 'kaiser_b3'
        w = kaiser_win(n_side, 3.0);
    case 'kaiser_b6'
        w = kaiser_win(n_side, 6.0);
    case 'chebyshev_25dB'
        w = chebwin_win(n_side, 25.0);
    case 'chebyshev_40dB'
        w = chebwin_win(n_side, 40.0);
    case 'taylor_25dB'
        w = taylor_win(n_side, 4, 25, true);
    otherwise
        error('window_1d:UnknownWindow', ...
            'Unknown classical technique ''%s''.', name);
end
w = w(:);
end


function w = hamming_win(M)
% numpy.hamming: 0.54 - 0.46*cos(2*pi*k/(M-1)).
if M == 1, w = 1; return; end
k = (0:M-1).';
w = 0.54 - 0.46 * cos(2 * pi * k / (M - 1));
end


function w = hanning_win(M)
% numpy.hanning: 0.5 - 0.5*cos(2*pi*k/(M-1)).
if M == 1, w = 1; return; end
k = (0:M-1).';
w = 0.5 - 0.5 * cos(2 * pi * k / (M - 1));
end


function w = kaiser_win(M, beta)
% numpy.kaiser: i0(beta*sqrt(1-((k-alpha)/alpha)^2)) / i0(beta).
if M == 1, w = 1; return; end
k = (0:M-1).';
alpha = (M - 1) / 2;
w = besseli(0, beta * sqrt(1 - ((k - alpha) / alpha) .^ 2)) / besseli(0, beta);
end


function w = chebwin_win(M, at)
% scipy.signal.windows.chebwin(M, at, sym=True). Dolph-Chebyshev window.
if M == 1, w = 1; return; end
order = M - 1.0;
beta  = cosh(1.0 / order * acosh(10 ^ (abs(at) / 20)));
k     = 0:M-1;
x     = beta * cos(pi * k / M);

p = zeros(1, M);
gt = x > 1;
lt = x < -1;
inr = abs(x) <= 1;
p(gt)  = cosh(order * acosh(x(gt)));
p(lt)  = (2 * mod(M, 2) - 1) * cosh(order * acosh(-x(lt)));
p(inr) = cos(order * acos(x(inr)));

if mod(M, 2) == 1
    w = real(fft(p));
    n = (M + 1) / 2;
    w = w(1:n);
    w = [w(n:-1:2), w];                 % [PORT] np w[n-1:0:-1] then full w[:n]
else
    p = p .* exp(1i * pi / M * (0:M-1));
    w = real(fft(p));
    n = M / 2 + 1;
    w = [w(n:-1:2), w(2:n)];            % [PORT] np w[n-1:0:-1] then w[1:n]
end
w = w / max(w);
w = w(:);
end


function w = taylor_win(M, nbar, sll, do_norm)
% scipy.signal.windows.taylor(M, nbar, sll, norm, sym=True).
if M == 1, w = 1; return; end
B  = 10 ^ (sll / 20);
A  = acosh(B) / pi;
s2 = nbar ^ 2 / (A ^ 2 + (nbar - 0.5) ^ 2);
ma = 1:(nbar - 1);          % row
m2 = ma .^ 2;

signs = ones(1, nbar - 1);
signs(2:2:end) = -1;

Fm = zeros(1, nbar - 1);
for mi = 1:(nbar - 1)
    numer  = signs(mi) * prod(1 - m2(mi) ./ s2 ./ (A ^ 2 + (ma - 0.5) .^ 2));
    others = (1:(nbar - 1)) ~= mi;
    denom  = 2 * prod(1 - m2(mi) ./ m2(others));
    Fm(mi) = numer / denom;
end

% W(n) = 1 + 2 * sum_m Fm_m cos(2*pi*m*(n - M/2 + 1/2)/M)
nvec = 0:(M - 1);
arg  = 2 * pi * (ma.') * ((nvec - M / 2 + 0.5) / M);   % (nbar-1 x M)
w    = 1 + 2 * (Fm * cos(arg));                        % (1 x M)

if do_norm
    nW   = (M - 1) / 2;
    argc = 2 * pi * (ma.') * ((nW - M / 2 + 0.5) / M);  % (nbar-1 x 1)
    Wc   = 1 + 2 * (Fm * cos(argc));
    w    = w / Wc;
end
w = w(:);
end
