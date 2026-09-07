function d_dbi = kpi_quiescent_directivity(run_log, n_el, stack1, stack2, ...
                                           theta_deg, phi_deg)
% KPI_QUIESCENT_DIRECTIVITY  Directivity [dBi] toward the target of the quiescent beam.
%
%   d_dbi = KPI_QUIESCENT_DIRECTIVITY(run_log, n_el, stack1, stack2, theta_deg, phi_deg)
%
%   The LCMV/MVDR solution for R = I — the beam the algorithm starts from before
%   any covariance is accumulated (closed_loop_run builds exactly this as the
%   oracle's step-1 weight). It is the reference dir_loss_db_ss is measured
%   against.
%
%   [P12] Extracted from run_amplitude_sweep_script's local
%   quiescent_directivity_dbi. IMPORTANT: this is a property of the ARRAY and
%   the steering direction, not of the campaign. The sweep script cached one
%   value from its first oracle run, which was correct only because every run in
%   it used the same array. A campaign spanning several arrays must call this
%   once per array, or every dir_loss_db_ss on the later arrays is silently
%   measured against the wrong reference.
%
%   Inputs:
%       run_log        : any run log on this array (only .grid is read, for e_s).
%       n_el           : element count. Units: dimensionless.
%       stack1, stack2 : (N_el x N_theta x N_phi) far-field stacks; stack2 is []
%                        for single-component operation. Units: V/m.
%       theta_deg, phi_deg : grid axes. Units: degrees.
%
%   Outputs:
%       d_dbi : quiescent directivity toward the target. Units: dBi.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P11, P12].

w_q     = adapt_lcmv(eye(n_el), run_log.grid.e_s, 0);
ref_log = struct('W', w_q, 'grid', run_log.grid);
d_dbi   = compute_directivity_trace(ref_log, stack1, stack2, theta_deg, phi_deg);
end
