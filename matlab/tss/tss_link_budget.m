function r = tss_link_budget(f_hz, cfg)
%   Layer:        L2
%   WP:           WP3
%   Claims:       C1
%   Status:       implemented
%   Requirements: INS-01, INS-10
%
%   Link budget for a spaceborne HF FMCW topside sounder.
%
%   Inputs
%     f_hz   sounding frequency, Hz.  Scalar.
%     cfg    configuration struct from TSS_CONFIG.
%
%   Outputs
%     r      struct of every budget term in dB plus the intermediate
%            temperatures, bandwidth and geometry.  r.snr_db is the figure
%            INS-10 is scored against.  r.status is 'ok', 'evanescent' or
%            'penetrates'; the last two return NaN power and SNR, which is a
%            result and not a failure.
%
%   Example
%     r = tss_link_budget(5e6, tss_config());
%     tss_print_budget(r);
%
%   TRANSCRIPTION.  Per docs/adr/0002, the Python src/tss/l2_instrument is the
%   reference for the link budget and this file is the transcription, even
%   though this file came first.  The two are compared to machine precision by
%   matlab/tests/TestTranscription.m and by tss verify --level V2.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_LINK_BUDGET  Link budget for a spaceborne HF FMCW topside sounder.
%
%   r = TSS_LINK_BUDGET(f_hz, cfg) returns a struct with every term in dB, so
%   the result can be printed as a budget table.  cfg comes from TSS_CONFIG.
%
%   This replaces snr_ionosonde.m.  Four things are done differently:
%
%   1. SKY NOISE.  At HF the noise floor is the galactic background, not the
%      receiver.  The default reference is Franklin and MacLean's in-situ
%      Alouette measurement.  A 300 K receiver temperature is about 25 dB
%      optimistic at 5 MHz; snr_ionosonde.m's effective 7900 K is about 20 dB
%      optimistic.
%
%   2. ANTENNA EFFICIENCY IS NOT A SIMPLE LOSS.  An electrically small antenna
%      couples to the sky with the same efficiency with which it radiates, so
%      T_a = eta*T_sky + (1-eta)*T_phys.  Where eta*T_sky dominates the receiver
%      contribution -- everywhere in the HF band for any efficiency above about
%      -30 dB -- the RECEIVE antenna efficiency cancels out of the SNR entirely.
%      Only the transmit efficiency costs anything.  Treating receive efficiency
%      as a straight loss double-counts it.
%
%   3. ABSORPTION is integrated properly (see TSS_ABSORPTION_DB) rather than
%      extrapolated as f^-2 from a single anchor.
%
%   4. NOISE BANDWIDTH is that of one range cell, B = 2*K_f*dR/c, not the swept
%      bandwidth.  This is the entire SWaP argument for the FMCW topology.
%      Using the swept bandwidth, as snr_ionosonde.m does, throws away about
%      27 dB of processing gain.
%
%   See also TSS_CONFIG, TSS_PRINT_BUDGET, SNR_IONOSONDE_V2.

c = tss_constants();
f_mhz = f_hz / 1e6;
lam   = c.c_light / f_hz;

% --- antenna: directivity and radiation efficiency kept separate -------------
eff_db = interp1(log10(cfg.eff_ref_mhz), cfg.eff_ref_db, log10(f_mhz), ...
                 'linear', 'extrap');
eta    = 10^(eff_db/10);
gain   = eff_db + cfg.directivity_dbi;

% --- geometry ---------------------------------------------------------------
[h_refl, status] = tss_reflection_height(f_hz, cfg);

% Specular reflection from a smooth overdense layer: the image source sits at
% twice the true range and the layer does not focus, so the geometric term is
% free-space spreading over a path of 2R.
path_db = 20*log10(lam / (4*pi*2*cfg.true_range));

if strcmp(status, 'ok')
    a_db = -tss_absorption_db(f_hz, cfg, h_refl);
else
    a_db = -tss_absorption_db(f_hz, cfg);
end

% --- signal -----------------------------------------------------------------
pt_dbw = 10*log10(cfg.tx_power);
eirp   = pt_dbw + gain;
pr_dbw = eirp + path_db + a_db + cfg.pol_loss_db + cfg.extra_loss_db ...
         + gain + cfg.feed_loss_db;

if ~strcmp(status, 'ok')
    % Evanescent: the wave never leaves the antenna.  Penetrating: there is no
    % topside reflection, only a ground echo, whose geometry is a different
    % calculation and is out of scope here.
    pr_dbw = NaN;
end

% --- noise ------------------------------------------------------------------
t_sky = tss_sky_temperature(f_mhz, cfg.sky_model);
t_ant = eta*t_sky + (1-eta)*cfg.ant_phys_temp;
l_feed = 10^(cfg.feed_loss_db/10);
if l_feed < 1
    t_feed = (1/l_feed - 1) * cfg.ant_phys_temp;
else
    t_feed = 0;
end
t_rx  = (10^(cfg.rx_noise_fig_db/10) - 1) * c.t0_ref;
t_sys = t_ant + t_feed + t_rx;

bw    = 2*cfg.sweep_rate*cfg.range_resolution / c.c_light;   % one range cell
n_dbw = 10*log10(c.k_boltz * t_sys * bw);

% --- result -----------------------------------------------------------------
gi = 10*log10(max(cfg.coherent_sweeps, 1));
r = struct( ...
    'f_mhz',            f_mhz, ...
    'status',           status, ...
    'tx_power_dbw',     pt_dbw, ...
    'tx_gain_dbi',      gain, ...
    'eirp_dbw',         eirp, ...
    'path_loss_db',     path_db, ...
    'absorption_db',    a_db, ...
    'pol_loss_db',      cfg.pol_loss_db, ...
    'extra_loss_db',    cfg.extra_loss_db, ...
    'rx_gain_dbi',      gain, ...
    'rx_power_dbw',     pr_dbw, ...
    'rx_efficiency_db', eff_db, ...
    'sky_temp_k',       t_sky, ...
    'ant_temp_k',       t_ant, ...
    'sys_temp_k',       t_sys, ...
    'noise_bw_hz',      bw, ...
    'dwell_s',          1/bw, ...
    'noise_dbw',        n_dbw, ...
    'snr_single_db',    pr_dbw - n_dbw, ...
    'integration_db',   gi, ...
    'snr_db',           pr_dbw - n_dbw + gi, ...
    'h_refl_km',        h_refl/1e3);
end
