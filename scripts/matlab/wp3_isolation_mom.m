function T = wp3_isolation_mom(outdir)
%WP3_ISOLATION_MOM  Moment-method anchor for the WP3 isolation profile (G1(b)).
%
%   T = wp3_isolation_mom(outdir) models two 5 m dipoles (radius 5 mm, the
%   baseline_6u antenna) in free space with the MATLAB Antenna Toolbox, takes
%   the two-port impedance matrix at 36 frequencies from 0.33 to 20 MHz, adds
%   each antenna's matching-coil loss |X|/Q (Q = 150, as tss ShortDipole), and
%   reports two isolations from the transmit port to the receive port:
%     iso_<geom>_db       inverse maximum available gain, i.e. both ports
%                         simultaneously conjugate-matched with the other
%                         antenna present: the least isolation the pair allows
%     iso_<geom>_self_db  inverse transducer gain with each port matched to
%                         its own self-impedance, as if tuned alone
%   and the same two without the coil loss (_lossless).
%
%   Two geometries, feeds 0.366 m apart (the long axis of a 6U bus):
%     parallel        both dipoles along x, feeds on the z axis: the worst case
%     crossed_offset  transmit along x, receive along y, receive feed displaced
%                     (0.03, 0.08) m from the axis, as a non-centred mount
%   The bus itself is not in the model (WP3 Phase 1 note).  Writes
%   isolation_mom.csv to outdir.  WP3 Phase 1, tss-sim, 2026-10-27.
if nargin < 1, outdir = pwd; end
f  = logspace(log10(0.33e6), log10(20e6), 36);
L  = 5; r = 5e-3; Q = 150; Z0 = 50; sep = 0.366;
w  = cylinder2strip(r);
dx = dipole('Length', L, 'Width', w, 'Tilt', 90, 'TiltAxis', [0 1 0]);
dy = dipole('Length', L, 'Width', w, 'Tilt', 90, 'TiltAxis', [1 0 0]);
geoms = {'parallel', {dx, dx}, [0 0 0; 0 0 sep]; ...
         'crossed_offset', {dx, dy}, [0 0 0; 0.03 0.08 sep]};
T = table(f(:), 'VariableNames', {'f_hz'});
for g = 1:size(geoms, 1)
    arr = conformalArray('Element', geoms{g, 2}, 'ElementPosition', geoms{g, 3});
    S = sparameters(arr, f, Z0);
    Z = s2z(S.Parameters, Z0);
    [isoM, isoS, isoM0, isoS0, eta] = deal(zeros(numel(f), 1));
    for k = 1:numel(f)
        Zk = Z(:, :, k);
        for lossy = [true false]
            Zp = Zk;
            if lossy
                Zp(1,1) = Zp(1,1) + abs(imag(Zk(1,1))) / Q;
                Zp(2,2) = Zp(2,2) + abs(imag(Zk(2,2))) / Q;
            end
            Zs = conj(Zp(1,1)); ZL = conj(Zp(2,2));
            Gs = 4 * real(Zs) * real(ZL) * abs(Zp(2,1))^2 / ...
                 abs((Zp(1,1) + Zs) * (Zp(2,2) + ZL) - Zp(1,2) * Zp(2,1))^2;
            m12 = Zp(1,2) * Zp(2,1);
            kk = (2 * real(Zp(1,1)) * real(Zp(2,2)) - real(m12)) / abs(m12);
            Gm = kk - sqrt(kk^2 - 1);
            if lossy, isoM(k) = -10*log10(Gm); isoS(k) = -10*log10(Gs);
            else,     isoM0(k) = -10*log10(Gm); isoS0(k) = -10*log10(Gs); end
        end
        eta(k) = 10 * log10(real(Zk(1,1)) / (real(Zk(1,1)) + abs(imag(Zk(1,1))) / Q));
    end
    T.(['iso_' geoms{g,1} '_db']) = isoM;
    T.(['iso_' geoms{g,1} '_self_db']) = isoS;
    T.(['iso_' geoms{g,1} '_lossless_db']) = isoM0;
    T.(['iso_' geoms{g,1} '_self_lossless_db']) = isoS0;
    if g == 1, T.eta_tx_db = eta; end
end
writetable(T, fullfile(outdir, 'isolation_mom.csv'));
end
