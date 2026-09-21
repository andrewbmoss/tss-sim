function c = tss_constants()
%   Layer:        support
%   WP:           WP3
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   Physical constants shared by the link budget functions.
%
%   Inputs
%     (none)
%
%   Outputs
%     c   struct with fields c_light (m/s), k_boltz (J/K), q_e (C), m_e (kg),
%         eps0 (F/m), t0_ref (K) and plasma_const, where
%         fN [Hz] = plasma_const * sqrt(Ne [m^-3]).
%
%   Example
%     c = tss_constants();
%     fof2_mhz = c.plasma_const * sqrt(1e12) / 1e6;
%
%   CODATA 2018 values.  These mirror src/tss/constants.py, which is the
%   reference implementation.  If the two ever disagree, Python is right and
%   this file is the defect.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_CONSTANTS  Physical constants shared by the link budget functions.
c.c_light  = 299792458.0;          % m/s
c.k_boltz  = 1.380649e-23;         % J/K
c.q_e      = 1.602176634e-19;      % C
c.m_e      = 9.1093837015e-31;     % kg
c.eps0     = 8.8541878128e-12;     % F/m
c.t0_ref   = 290.0;                % K
% fN = plasma_const * sqrt(Ne),  Ne in m^-3, fN in Hz
c.plasma_const = sqrt(c.q_e^2 / (c.eps0 * c.m_e)) / (2*pi);
end
