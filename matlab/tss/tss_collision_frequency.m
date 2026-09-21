function nu = tss_collision_frequency(ne, te)
%   Layer:        L0
%   WP:           WP1
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   Electron collision frequency, s^-1.
%
%   Inputs
%     ne    electron density, m^-3.
%     te    electron temperature, K.  Same size as ne, or scalar.
%
%   Outputs
%     nu    collision frequency, s^-1.
%
%   Example
%     nu = tss_collision_frequency(5e11, 2000.0);   % 304.66 s^-1
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_COLLISION_FREQUENCY  Electron collision frequency, s^-1.
%
%   Electron-ion Coulomb collisions dominate in the topside F region.  The
%   Schunk and Nagy form nu_ei = 54.5 * n_i[cm^-3] / Te^1.5 is used with
%   quasi-neutrality n_i = n_e.  Electron-neutral collisions are negligible
%   above about 300 km and are omitted.
nu = 54.5 .* (ne * 1e-6) ./ te.^1.5;
end
