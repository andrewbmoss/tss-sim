function t = tss_sky_temperature(f_mhz, model)
%   Layer:        L2
%   WP:           WP3, WP12
%   Claims:       C1
%   Status:       implemented
%   Requirements: INS-10
%
%   Brightness temperature of the HF sky, in kelvin.
%
%   Inputs
%     f_mhz   frequency or vector of frequencies, MHz.
%     model   'franklin' (default), 'itu', 'legacy300' or 'legacy_sni'.
%
%   Outputs
%     t       sky brightness temperature, K, the same size as f_mhz.
%
%   Example
%     t = tss_sky_temperature(5.0);            % 1.25e6 K
%     t = tss_sky_temperature(0.5:0.1:20, 'itu');
%
%   The Franklin and MacLean reference is measured from 0.6 to 12 MHz.  Below
%   0.6 MHz this function extrapolates a power law, and that extrapolation is
%   exactly the open item WP12 exists to close.  Any band edge result that
%   leans on it carries the extrapolated-noise flag.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_SKY_TEMPERATURE  Brightness temperature of the HF sky, in kelvin.
%
%   'franklin'   (default) Franklin and MacLean (1969), Proc. IEEE 57(6).
%                In-situ measurement from Alouette I and II: 2e7 K at 0.6 MHz
%                falling to 4e5 K at 12 MHz, +/- 2 dB.  Log-log interpolated.
%                This is the recommended reference for a spaceborne HF sounder:
%                it was measured in the exact environment in question.
%   'itu'        ITU-R P.372 galactic term, Fa = 52 - 23*log10(f_MHz) dB above
%                kT0.  A ground-based model, included for comparison.
%   'legacy300'  Flat 300 K -- the assumption in the ECE 6930 paper.
%   'legacy_sni' Flat 250 K raised 15 dB (~7900 K) -- snr_ionosonde.m.
%
%   Solar bursts (>30 dB) and ionospheric noise bands (up to 80 dB) are not
%   included; both are availability questions rather than link budget terms.

if nargin < 2, model = 'franklin'; end
c = tss_constants();
switch lower(model)
    case 'franklin'
        f1 = 0.6;  t1 = 2.0e7;
        f2 = 12.0; t2 = 4.0e5;
        slope = log10(t2/t1) / log10(f2/f1);      % ~ -1.306
        t = t1 .* (f_mhz ./ f1).^slope;
    case 'itu'
        fa_db = 52.0 - 23.0*log10(f_mhz);
        t = c.t0_ref .* 10.^(fa_db/10);
    case 'legacy300'
        t = 300.0 * ones(size(f_mhz));
    case 'legacy_sni'
        t = 250.0 * 10^(15/10) * ones(size(f_mhz));
    otherwise
        error('tss_sky_temperature:model', 'unknown sky model "%s"', model);
end
end
