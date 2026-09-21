function a_db = tss_absorption_db(f_hz, cfg, h_refl, n_steps)
%   Layer:        L1
%   WP:           WP2
%   Claims:       C1
%   Status:       implemented
%   Requirements: INS-10
%
%   Two-way ionospheric absorption, dB, for a vertical sounding.
%
%   Inputs
%     f_hz     sounding frequency, Hz.  Scalar.
%     cfg      configuration struct.
%     h_refl   reflection altitude, m.  Optional; computed if omitted.
%     n_steps  integration steps.  Default 3000.
%
%   Outputs
%     a_db     two-way absorption, dB, positive.  NaN when the wave is
%              evanescent at the spacecraft.
%
%   Example
%     a = tss_absorption_db(5e6, tss_config());   % 0.783 dB
%
%   This integral declares its own tolerance in the transcription verifier.
%   It is one of the two cases named in docs/matlab-python-bridge.md that
%   cannot meet the default 1e-12 relative tolerance, because the two sides
%   use quadrature with their own convergence behaviour rather than a closed
%   form.  The declared tolerance is 1e-6 relative and the justification is
%   this sentence.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_ABSORPTION_DB  Two-way ionospheric absorption, dB, for a vertical sounding.
%
%   The absorption coefficient is kappa = (omega/c) * |Im(n)| in Np/m, with the
%   complex refractive index
%
%       n^2 = 1 - X / (1 - i Z),   X = (fN/f)^2,   Z = nu/omega
%
%   integrated over the physical path from the spacecraft to the reflection
%   height and doubled.  Because the complex index is used, the integrand stays
%   finite at the reflection point, where X -> 1 and |Im(n)| -> sqrt(Z/2).
%
%   This replaces the f^-2 scaling from a single 5 MHz anchor used in the
%   ECE 6930 paper.  That extrapolation diverges below about 2 MHz -- it returns
%   125 dB at 1 MHz where the true figure is a fraction of a decibel -- and it
%   was the source of the apparent low-frequency shortfall in earlier analyses.
%   Set cfg.absorption_legacy = true to reproduce it for comparison only.
%
%   The magnetic field is neglected, so the O and X modes are not separated
%   here; that belongs in the propagation model, not in a link budget.

if nargin < 4, n_steps = 3000; end
c = tss_constants();

if cfg.absorption_legacy
    a_db = 5.0 * (5.0 / (f_hz/1e6))^2;      % legacy anchor: 5 dB two-way at 5 MHz
    return
end

if nargin < 3 || isempty(h_refl) || ~isfinite(h_refl)
    [h_refl, status] = tss_reflection_height(f_hz, cfg);
    if strcmp(status, 'evanescent')
        a_db = NaN;  return
    end
    if strcmp(status, 'penetrates')
        h_refl = cfg.hmf2;
    end
end

% Non-uniform grid clustered toward the reflection point: |Im(n)| varies as the
% square root of the distance from reflection there, so a uniform grid converges
% slowly and non-monotonically.  A cubic stretch fixes it at the same cost.
s = linspace(0, 1, n_steps).^3;
h = h_refl + (cfg.h_sat - h_refl) .* s;

ne = tss_ne(h, cfg);
te = tss_te(h, cfg);
nu = tss_collision_frequency(ne, te);

omega = 2*pi*f_hz;
fn2   = (c.plasma_const^2) .* ne;            % fN^2
X     = fn2 ./ (f_hz^2);
Z     = nu ./ omega;
n_cpx = sqrt(1 - X ./ (1 - 1i*Z));
kappa = (omega / c.c_light) .* abs(imag(n_cpx));     % Np/m

one_way_np = trapz(h, kappa);
a_db = 2.0 * abs(one_way_np) * 8.685889638;
end
