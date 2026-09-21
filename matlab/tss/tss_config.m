function cfg = tss_config(varargin)
%   Layer:        L2
%   WP:           WP3
%   Claims:       C1, C4
%   Status:       implemented
%   Requirements: INS-01, INS-10
%
%   Default configuration struct for the link budget.
%
%   Inputs
%     varargin   optional 'name', value pairs overriding any field.  An
%                unknown field name is an error rather than a silent no-op.
%
%   Outputs
%     cfg        configuration struct.  Fields are documented below.
%
%   Example
%     cfg = tss_config();
%     cfg = tss_config('tx_power', 0.5, 'sky_model', 'itu');
%
%   This is not the YAML campaign configuration.  It is the small struct the
%   migrated link budget functions take.  A study that needs to be driven from
%   config/ reads the YAML with tss_yaml_lite and maps the fields it needs onto
%   this struct; the mapping is deliberately explicit rather than automatic, so
%   that an unmapped YAML key is visible instead of quietly ignored.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_CONFIG  Default configuration for the FMCW topside sounder link budget.
%
%   cfg = TSS_CONFIG() returns the reference case: a 600 km circular orbit, a
%   moderate daytime mid-latitude ionosphere, and the TISC point-design waveform
%   corrected to the science requirements (1 km range cell, 1000 km maximum
%   apparent range).
%
%   cfg = TSS_CONFIG('name', value, ...) overrides any field.
%
%   Fields
%     sweep_rate        Hz/s      linear FM sweep rate K_f
%     range_resolution  m         delta R; sets the noise bandwidth
%     coherent_sweeps   -         M sweeps integrated coherently
%     h_sat             m         spacecraft altitude
%     true_range        m         spacecraft to reflection point
%     tx_power          W
%     pol_loss_db       dB        O/X split plus elliptical-to-linear (negative)
%     extra_loss_db     dB        anything else (negative)
%     rx_noise_fig_db   dB
%     feed_loss_db      dB        antenna to LNA (negative)
%     directivity_dbi   dBi       small loop / short dipole
%     eff_ref_mhz       MHz       two-point radiation efficiency reference
%     eff_ref_db        dB        radiation efficiency at those frequencies
%     ant_phys_temp     K
%     sky_model         char      'franklin' | 'itu' | 'legacy300' | 'legacy_sni'
%     nmf2              m^-3      F2 peak density
%     hmf2              m         F2 peak height
%     plasma_scale_ht   m         topside PLASMA scale height Hp
%     te_peak, te_top   K         electron temperature at hmF2 and at 1000 km
%     absorption_legacy logical   use the discredited f^-2 scaling (comparison only)
%
%   See also TSS_LINK_BUDGET, SNR_IONOSONDE_V2.
%
%   Andrew B. Moss -- PhD dissertation, Utah State University.

cfg = struct( ...
    'sweep_rate',        5.6e6, ...
    'range_resolution',  1000.0, ...
    'coherent_sweeps',   1, ...
    'h_sat',             600e3, ...
    'true_range',        300e3, ...
    'tx_power',          1.0, ...
    'pol_loss_db',      -6.0, ...
    'extra_loss_db',     0.0, ...
    'rx_noise_fig_db',   3.0, ...
    'feed_loss_db',      0.0, ...
    'directivity_dbi',   1.76, ...
    'eff_ref_mhz',       [1.0 20.0], ...
    'eff_ref_db',        [-11.76 -5.76], ...
    'ant_phys_temp',     290.0, ...
    'sky_model',         'franklin', ...
    'nmf2',              1.0e12, ...
    'hmf2',              300e3, ...
    'plasma_scale_ht',   100e3, ...
    'te_peak',           1500.0, ...
    'te_top',            3000.0, ...
    'absorption_legacy', false);

for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~isfield(cfg, name)
        error('tss_config:unknownField', 'unknown configuration field "%s"', name);
    end
    cfg.(name) = varargin{k+1};
end
end
