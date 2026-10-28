function [T, sens] = wp3_isolation_mom_bus(outdir)
%WP3_ISOLATION_MOM_BUS  Moment-method isolation of two 5 m dipoles on a 6U bus.
%
%   [T, sens] = wp3_isolation_mom_bus(outdir) extends wp3_isolation_mom with the
%   spacecraft body: a closed conducting box 0.100 x 0.2263 x 0.366 m (6U),
%   meshed at 1.5 cm (wp3_box_mesh), solved together with the two dipoles by
%   full MoM (installedAntenna, SolverType 'MoM').  Each dipole (5 m, radius
%   5 mm) is fed 5 cm off a bus face, as on an insulating standoff.  Isolation is
%   computed as in wp3_isolation_mom: coil loss |X|/Q at Q = 150 in series with
%   each port, then 1/MAG (both ports conjugate-matched together) and 1/G_T
%   with each port self-matched.
%
%   Mounts (bus from z = 0 to 0.366 m; x is the 10 cm side, y the 22.6 cm side):
%     crossed_ends_offset   tx along x below the bus, rx along y above it, rx
%                           feed displaced (3, 8) cm: the free-space anchor's
%                           geometry with the bus added
%     crossed_ends_centred  as above with the rx feed on the axis
%     parallel_ends         both along x, one below and one above the bus
%     crossed_colocated     both below the bus, tx along x at z = -5 cm, rx along
%                           y at z = -8 cm (a turnstile pair)
%   Sensitivity (sens), at 0.33, 1.26, 5 and 20 MHz, each imperfection alone
%   on the two symmetric crossed mounts (crossed_ends_centred and
%   crossed_colocated): the rx dipole rotated from y toward x by a deployment
%   error of 0 to 3 deg, and the rx feed displaced along y (along its own axis,
%   which breaks the mirror symmetry) by 0 to 10 cm.  The symmetric mounts with
%   no imperfection sit at the numerical floor (k ~ 5e-6).
%
%   Mesh: 5 cm on the end mounts (coupling converged to 2 % against 10 cm);
%   3 cm on the co-located pair, whose crossing 3 cm apart converges only to
%   about 3 dB in coupling (k) between 5 and 3 cm; its symmetric case is at the
%   numerical floor (k ~ 1e-5) and is not a physical number.
%
%   Writes isolation_mom_bus.csv and isolation_mom_bus_sensitivity.csv.
%   WP3, tss-sim, 2026-10-27.
if nargin < 1, outdir = pwd; end
here = fileparts(mfilename('fullpath')); addpath(here);
f  = logspace(log10(0.33e6), log10(20e6), 36);
L  = 5; r = 5e-3; Q = 150; Z0 = 50; bx = [0.1 0.2263 0.366]; g = 0.05;
stl = fullfile(tempdir, 'wp3_bus6u.stl');
stlwrite(wp3_box_mesh(bx, 0.015), stl);
plat = platform('FileName', stl, 'Units', 'm');
below = [0 0 -g]; above = [0 0 bx(3) + g];
mounts = {
  'crossed_ends_offset',  0, 90, below, above + [0.03 0.08 0], 0.05
  'crossed_ends_centred', 0, 90, below, above,                  0.05
  'parallel_ends',        0,  0, below, above,                  0.05
  'crossed_colocated',    0, 90, below, below + [0 0 -0.03],    0.03
};
T = table(f(:), 'VariableNames', {'f_hz'});
for m = 1:size(mounts, 1)
    [isoM, isoS, kc] = solve(mounts(m, :), f);
    T.(['iso_' mounts{m,1} '_db']) = isoM;
    T.(['iso_' mounts{m,1} '_self_db']) = isoS;
    T.(['k_' mounts{m,1}]) = kc;
end
writetable(T, fullfile(outdir, 'isolation_mom_bus.csv'));

fs = [0.33e6 1.26e6 5e6 20e6];
rows = {};
for base = [2 4]
    for err = [0 0.1 0.3 1 3]
        mt = mounts(base, :); mt{3} = 90 - err;
        [isoM, ~, kc] = solve(mt, fs);
        rows(end+1, :) = {mounts{base,1}, 'angle_deg', err, isoM(1), isoM(2), isoM(3), isoM(4), kc(1)}; %#ok<AGROW>
    end
    for off = [0 0.01 0.03 0.1]
        mt = mounts(base, :); mt{5} = mt{5} + [0 off 0];
        [isoM, ~, kc] = solve(mt, fs);
        rows(end+1, :) = {mounts{base,1}, 'feed_offset_m', off, isoM(1), isoM(2), isoM(3), isoM(4), kc(1)}; %#ok<AGROW>
    end
end
sens = cell2table(rows, 'VariableNames', {'mount', 'parameter', 'value', ...
    'iso_0p33MHz_db', 'iso_1p26MHz_db', 'iso_5MHz_db', 'iso_20MHz_db', 'k_0p33MHz'});
writetable(sens, fullfile(outdir, 'isolation_mom_bus_sensitivity.csv'));

    function [isoM, isoS, kc] = solve(mt, freqs)
        w = cylinder2strip(r);
        tx = dipole('Length', L, 'Width', w, 'Tilt', 90, 'TiltAxis', [0 1 0]);   % along x
        % rx: along x rotated about z by the given angle (90 = along y)
        rx = dipole('Length', L, 'Width', w, 'Tilt', [90 mt{3}], 'TiltAxis', [0 1 0; 0 0 1]);
        if mt{2} ~= 0, tx.Tilt = [90 mt{2}]; tx.TiltAxis = [0 1 0; 0 0 1]; end
        ia = installedAntenna('Platform', plat, 'Element', {tx, rx}, ...
            'ElementPosition', [mt{4}; mt{5}], 'SolverType', 'MoM');
        mesh(ia, 'MaxEdgeLength', mt{6});
        S = sparameters(ia, freqs, Z0);
        Z = s2z(S.Parameters, Z0);
        [isoM, isoS, kc] = deal(zeros(numel(freqs), 1));
        for k = 1:numel(freqs)
            Zp = Z(:, :, k);
            kc(k) = abs(Zp(1,2)) / sqrt(abs(imag(Zp(1,1))) * abs(imag(Zp(2,2))));
            Zp(1,1) = Zp(1,1) + abs(imag(Zp(1,1))) / Q;
            Zp(2,2) = Zp(2,2) + abs(imag(Zp(2,2))) / Q;
            Zs = conj(Zp(1,1)); ZL = conj(Zp(2,2));
            Gs = 4 * real(Zs) * real(ZL) * abs(Zp(2,1))^2 / ...
                 abs((Zp(1,1) + Zs) * (Zp(2,2) + ZL) - Zp(1,2) * Zp(2,1))^2;
            m12 = Zp(1,2) * Zp(2,1);
            kk = (2 * real(Zp(1,1)) * real(Zp(2,2)) - real(m12)) / abs(m12);
            isoM(k) = 10 * log10(kk + sqrt(kk^2 - 1));   % 1/MAG, stable for weak coupling
            isoS(k) = -10 * log10(Gs);
        end
    end
end
