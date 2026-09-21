function [snrR, pwrR, nR] = snr_ionosonde_v2(Kf, R, Rmax, f, Pt, cfg)
%   Layer:        L2
%   WP:           WP3
%   Claims:       C1
%   Status:       implemented
%   Requirements: INS-10
%
%   Drop-in replacement for the 2014 snr_ionosonde.m.
%
%   Inputs
%     Kf     sweep rate, Hz/s.
%     R      reflection range, m.
%     Rmax   maximum reflection range, m.  Ambiguity check only.
%     f      signal frequency, Hz.
%     Pt     transmit power, W.
%     cfg    optional configuration struct from TSS_CONFIG.
%
%   Outputs
%     snrR   signal to noise ratio, dB.
%     pwrR   received signal power, dBW.
%     nR     noise power, dBW.
%
%   Example
%     [s, p, n] = snr_ionosonde_v2(5.6e6, 300e3, 500e3, 5e6, 1.0);
%
%   Results are memoized on a persistent map.  Call CLEAR SNR_IONOSONDE_V2
%   after changing cfg.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%SNR_IONOSONDE_V2  Drop-in replacement for the 2014 snr_ionosonde.m.
%
%   [snrR, pwrR, nR] = SNR_IONOSONDE_V2(Kf, R, Rmax, f, Pt) keeps the original
%   call signature so existing scripts (TssModel3.m and friends) run unchanged:
%   replace the function name and nothing else.
%
%   [snrR, pwrR, nR] = SNR_IONOSONDE_V2(Kf, R, Rmax, f, Pt, cfg) uses the
%   configuration struct from TSS_CONFIG for everything the old signature could
%   not express (ionospheric profile, sky model, antenna description).
%
%   ARGUMENTS  (unchanged from the original)
%     Kf    sweep rate of the FMCW signal, Hz/s
%     R     reflection range, m
%     Rmax  maximum reflection range, m
%     f     signal frequency, Hz
%     Pt    transmit power, W
%
%   OUTPUT     (unchanged from the original)
%     snrR  signal-to-noise ratio of the received signal, dB
%     pwrR  received signal power, dBW
%     nR    noise power, dBW
%
%   *** WHAT CHANGED, AND WHY THE NUMBERS MOVE ***
%
%   The original computed the noise bandwidth as BW = 2*Kf*Rmax/c -- the whole
%   swept band -- and therefore discarded the FMCW processing gain that is the
%   entire reason for choosing this waveform.  Here the noise bandwidth is that
%   of one range cell, B = 2*Kf*dR/c, taken from cfg.range_resolution (1 km by
%   default).  Rmax is now used only to check for range ambiguity, and a warning
%   is issued if R exceeds it.  This alone moves the answer by about +27 dB.
%
%   Against that, the original assumed a 250 K antenna temperature plus a 15 dB
%   "environmental noise" term, an effective 7900 K.  At HF the sky is the
%   noise floor: Franklin and MacLean measured 2e7 K at 0.6 MHz falling to
%   4e5 K at 12 MHz from Alouette.  Correcting this costs about 20 dB.
%
%   The original also treated the receive antenna efficiency as a straight loss
%   while assuming a fixed antenna temperature.  Because an electrically small
%   antenna couples to the sky with the same efficiency with which it radiates,
%   that double-counts it; see TSS_LINK_BUDGET.
%
%   Net effect at 5 MHz, 1 W, 300 km range: the original reports 17.7 dB, this
%   reports about 38.1 dB.  The original was pessimistic, not optimistic.
%
%   PERFORMANCE.  The absorption integral is far too slow to run per sample
%   inside a simulation loop, so results are memoized on frequency (to 10 kHz)
%   and range (to 1 km).  Call CLEAR SNR_IONOSONDE_V2 to reset the cache after
%   changing cfg.
%
%   See also TSS_LINK_BUDGET, TSS_CONFIG.

persistent memo memo_cfg

if nargin < 6 || isempty(cfg)
    cfg = tss_config();
end
if isempty(memo) || ~isequal(memo_cfg, cfg)
    memo = containers.Map('KeyType', 'char', 'ValueType', 'any');
    memo_cfg = cfg;
end

if R > Rmax
    warning('snr_ionosonde_v2:ambiguous', ...
            'reflection range %.1f km exceeds Rmax %.1f km: the echo aliases', ...
            R/1e3, Rmax/1e3);
end

key = sprintf('%.0f_%.0f', round(f/1e4), round(R/1e3));
if isKey(memo, key)
    v = memo(key);
else
    c = cfg;
    c.sweep_rate = Kf;
    c.tx_power   = Pt;
    c.true_range = R;
    r = tss_link_budget(f, c);
    v = [r.snr_db, r.rx_power_dbw, r.noise_dbw];
    memo(key) = v;
end

snrR = v(1);  pwrR = v(2);  nR = v(3);
end
