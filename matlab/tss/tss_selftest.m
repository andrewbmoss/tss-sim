function tss_selftest()
%   Layer:        test
%   WP:           WP8
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   Verify the MATLAB link budget against the Python reference values.
%
%   Inputs
%     (none)
%
%   Outputs
%     (none)  prints PASS or FAIL per check and a summary line.
%
%   Example
%     tss_selftest
%
%   RETAINED AS THE AUTHOR WROTE IT.  Every assertion here is also ported to
%   matlab/tests/TestLinkBudget.m, which is what CI runs.  This file is kept
%   because it is runnable in one line from a bare MATLAB prompt with no test
%   framework, which is what it is for.  If a check is added in one place it
%   must be added in the other; TestLinkBudget names this file as its source.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_SELFTEST  Verify the MATLAB link budget against the Python reference.
%
%   The expected values below were produced by the Python implementation
%   (tss_link_budget.py) on the reference case defined in TSS_CONFIG.  The two
%   implementations are independent transcriptions of the same equations, so
%   agreement to the tolerances below is a real cross-check, not a tautology.
%
%   Run:  tss_selftest
%
%   A failure prints the expected and actual values and the difference.

fprintf('tss_selftest -- MATLAB vs Python reference\n');
fprintf('%s\n', repmat('-', 1, 66));

cfg = tss_config();
c   = tss_constants();
n_fail = 0;

% --- geometry and profile ---------------------------------------------------
n_fail = check('foF2 (MHz)', c.plasma_const*sqrt(cfg.nmf2)/1e6, 8.978663, 1e-5, n_fail);
n_fail = check('Ne at 600 km (m^-3)', tss_ne(600e3, cfg), 8.198333e10, 1e-4, n_fail);
n_fail = check('local fN at 600 km (MHz)', ...
    c.plasma_const*sqrt(tss_ne(cfg.h_sat, cfg))/1e6, 2.570836, 1e-5, n_fail);

% --- environment ------------------------------------------------------------
n_fail = check('sky temp, Franklin, 5 MHz (K)', ...
    tss_sky_temperature(5.0, 'franklin'), 1.254772e6, 1e-5, n_fail);
n_fail = check('sky temp, ITU-R, 5 MHz (K)', ...
    tss_sky_temperature(5.0, 'itu'), 1.134402e6, 1e-5, n_fail);
n_fail = check('collision freq, 5e11 m^-3, 2000 K (s^-1)', ...
    tss_collision_frequency(5e11, 2000.0), 304.664262, 1e-6, n_fail);

% --- reflection height and absorption ---------------------------------------
n_fail = check('h_refl at 3 MHz (km)', tss_reflection_height(3e6, cfg)/1e3, 569.0174, 1e-4, n_fail);
n_fail = check('h_refl at 5 MHz (km)', tss_reflection_height(5e6, cfg)/1e3, 465.2476, 1e-4, n_fail);
n_fail = check('h_refl at 8 MHz (km)', tss_reflection_height(8e6, cfg)/1e3, 357.1336, 1e-4, n_fail);

n_fail = check('absorption at 3 MHz (dB)', tss_absorption_db(3e6, cfg), 0.17303, 1e-3, n_fail);
n_fail = check('absorption at 5 MHz (dB)', tss_absorption_db(5e6, cfg), 0.78260, 1e-3, n_fail);
n_fail = check('absorption at 8 MHz (dB)', tss_absorption_db(8e6, cfg), 3.07232, 1e-3, n_fail);

% --- full budget ------------------------------------------------------------
n_fail = check('noise bandwidth (Hz)', ...
    2*cfg.sweep_rate*cfg.range_resolution/c.c_light, 37.359179, 1e-6, n_fail);

r3 = tss_link_budget(3e6, cfg);
r5 = tss_link_budget(5e6, cfg);
r8 = tss_link_budget(8e6, cfg);
n_fail = check('SNR at 3 MHz (dB)', r3.snr_db, 39.21783, 1e-3, n_fail);
n_fail = check('SNR at 5 MHz (dB)', r5.snr_db, 38.08692, 1e-3, n_fail);
n_fail = check('SNR at 8 MHz (dB)', r8.snr_db, 35.31559, 1e-3, n_fail);
n_fail = check('Pr at 5 MHz (dBW)',  r5.rx_power_dbw, -122.32589, 1e-3, n_fail);
n_fail = check('N at 5 MHz (dBW)',   r5.noise_dbw,    -160.41281, 1e-3, n_fail);
n_fail = check('Tsys at 5 MHz (K)',  r5.sys_temp_k,   1.762943e5, 1e-4, n_fail);

% --- legacy regression ------------------------------------------------------
[sl, pl, nl] = snr_ionosonde_legacy(5.6e6, 300e3, 500e3, 5e6, 1.0);
n_fail = check('legacy SNR at 5 MHz (dB)', sl, 17.65723, 1e-3, n_fail);
n_fail = check('legacy Pr at 5 MHz (dBW)', pl, -111.80021, 1e-3, n_fail);
n_fail = check('legacy N at 5 MHz (dBW)',  nl, -129.45744, 1e-3, n_fail);

% --- limiting-case behaviour ------------------------------------------------
[~, st_lo] = tss_reflection_height(1.0e6, cfg);
n_fail = check_str('1 MHz is below the local plasma frequency', st_lo, 'evanescent', n_fail);
[~, st_hi] = tss_reflection_height(12.0e6, cfg);
n_fail = check_str('12 MHz penetrates the F2 peak', st_hi, 'penetrates', n_fail);

r_nan = tss_link_budget(1.0e6, cfg);
n_fail = check_true('evanescent case returns NaN SNR', isnan(r_nan.snr_db), n_fail);

% Receive efficiency should very nearly cancel while the sky dominates: halving
% it must change the SNR by close to the transmit-side 3 dB only.
cfg_lo = cfg;  cfg_lo.eff_ref_db = cfg.eff_ref_db - 3.0;
d = tss_link_budget(5e6, cfg).snr_db - tss_link_budget(5e6, cfg_lo).snr_db;
n_fail = check('3 dB efficiency drop costs (dB)', d, 3.0, 0.05, n_fail);

fprintf('%s\n', repmat('-', 1, 66));
if n_fail == 0
    fprintf('ALL CHECKS PASS\n');
else
    fprintf('%d CHECK(S) FAILED\n', n_fail);
end
end

% ---------------------------------------------------------------------------
function n_fail = check(name, actual, expected, rtol, n_fail)
if expected == 0
    ok = abs(actual) <= rtol;
else
    ok = abs(actual - expected) <= rtol*max(1, abs(expected));
end
if ok
    fprintf('  PASS  %-42s %14.6g\n', name, actual);
else
    fprintf('  FAIL  %-42s %14.6g  expected %14.6g  (diff %.3g)\n', ...
            name, actual, expected, actual-expected);
    n_fail = n_fail + 1;
end
end

function n_fail = check_str(name, actual, expected, n_fail)
if strcmp(actual, expected)
    fprintf('  PASS  %-42s %14s\n', name, actual);
else
    fprintf('  FAIL  %-42s %14s  expected %s\n', name, actual, expected);
    n_fail = n_fail + 1;
end
end

function n_fail = check_true(name, cond, n_fail)
if cond
    fprintf('  PASS  %-42s %14s\n', name, 'true');
else
    fprintf('  FAIL  %-42s %14s\n', name, 'false');
    n_fail = n_fail + 1;
end
end
