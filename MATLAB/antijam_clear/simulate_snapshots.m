function snapshot_block = simulate_snapshots(scenario, step_index, array)
% ══════════════════════════════════════════════════════════════════
% SIMULATE_SNAPSHOTS
% Turn the truth into the only data an algorithm is allowed to see.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   snapshot_block = SIMULATE_SNAPSHOTS(scenario, step_index, array)
%
%   One step of a digital beamforming receiver: K simultaneous complex samples
%   from every element,
%
%       x = e_s * s  +  e_j * j  +  n
%
%   where e_s and e_j are the array's responses in the signal and jammer
%   directions, s and j are the two transmissions, and n is receiver noise,
%   independent on every element. With one polarization component each source is
%   RANK ONE: a single scalar amplitude times a single steering vector. That is
%   what makes the covariance structure, and every explanation of it, simple.
%
%   This function is the boundary. Truth goes in; only snapshots come out. The
%   jammer's direction is used here to BUILD the data and is not returned.
%
%   THE FACT THAT MATTERS MOST IN THIS WHOLE FOLDER: the wanted signal is inside
%   the snapshots. Real receivers cannot separate it out -- that is what you are
%   trying to receive. So the covariance an algorithm estimates contains the
%   signal as well as the jammer, which makes the beamformer MPDR rather than
%   MVDR, and it has one serious consequence:
%
%       if the assumed signal direction is even slightly wrong, the solver sees
%       wanted signal it cannot account for and spends degrees of freedom
%       CANCELLING YOUR OWN TRANSMISSION -- and the stronger your signal, the
%       more of it there is to cancel.
%
%   This is why diagonal loading in max_sinr_weights is structural rather than a
%   tuning knob: it is the thing standing between the solver and self-
%   cancellation. Previous measurements put the cost at up to 8.1 dB of a 20 dB
%   signal increase on this very array; it is measured here rather than assumed.
%
%   WHY K IS DERIVED. K comes from the array as array.snapshots_per_step: a
%   covariance estimated from K samples of an n-element array lands within about
%   3 dB of the ideal once K is roughly 2 * n_elements (the standard
%   Reed-Mallett-Brennan result), so it is 32 snapshots for a 16-element array.
%   The question "why 32?" has an answer that is not a shrug. See make_array.
%
%   RANDOMNESS. Draws come from scenario.random_stream, a seeded handle object
%   that advances in place, so a run is reproducible from its seed and two
%   consecutive steps are statistically independent.
%
%   [FUTURE] Two-component operation: e_s and e_j become (n_elements x 2) and
%   each source's power is split equally across its two columns.
%
%   Inputs:
%       scenario   : struct from make_scenario. TRUTH -- see that file.
%       step_index : 1-based step to simulate, 1 .. scenario.n_steps.
%       array      : struct from make_array. Supplies the element patterns, the
%                    angle grids and K = snapshots_per_step.
%
%   Outputs:
%       snapshot_block : (n_elements x K) complex received samples.
%                        Units: V/m times dimensionless amplitude.

if step_index < 1 || step_index > scenario.n_steps
    error('simulate_snapshots:StepOutOfRange', ...
        'Step %d requested; scenario has %d steps.', step_index, scenario.n_steps);
end

element_patterns = array.element_patterns;
theta_deg        = array.theta_deg;
phi_deg          = array.phi_deg;
n_elements       = array.n_elements;
n_snapshots      = array.snapshots_per_step;
stream           = scenario.random_stream;

% ────────────────────────── THE WANTED SIGNAL ─────────────────────

% Fixed direction throughout the run.
signal_steering = steering_vector(element_patterns, theta_deg, phi_deg, ...
                                  scenario.signal_theta_deg, scenario.signal_phi_deg);

signal_amplitude = circular_gaussian(stream, 1, n_snapshots, scenario.signal_power);
snapshot_block   = signal_steering * signal_amplitude;

% ────────────────────────── THE JAMMER ────────────────────────────

% Its direction is read at THIS step, so a drifting jammer moves between steps
% and an on/off jammer contributes nothing at all while it is switched off.
if scenario.jammer_on(step_index)
    jammer_steering = steering_vector(element_patterns, theta_deg, phi_deg, ...
                                      scenario.jammer_theta_deg(step_index), ...
                                      scenario.jammer_phi_deg(step_index));

    jammer_amplitude = circular_gaussian(stream, 1, n_snapshots, scenario.jammer_power);
    snapshot_block   = snapshot_block + jammer_steering * jammer_amplitude;
end

% ────────────────────────── RECEIVER NOISE ────────────────────────

% Independent on every element: this is the term that makes the covariance
% full rank and therefore invertible at all.
snapshot_block = snapshot_block ...
                 + circular_gaussian(stream, n_elements, n_snapshots, scenario.noise_power);
end


% ────────────────────────── HELPERS ───────────────────────────────

function samples = circular_gaussian(stream, n_rows, n_columns, power)
% Circular complex Gaussian samples of the given total power.
%
% The power is split equally between the real and imaginary parts, so
% E[|sample|^2] = power. This is the standard model for a narrowband
% communications or noise waveform observed in complex baseband.
samples = sqrt(power / 2) * (randn(stream, n_rows, n_columns) ...
                             + 1i * randn(stream, n_rows, n_columns));
end
