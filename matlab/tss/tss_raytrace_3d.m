function [rays, info] = tss_raytrace_3d(state, platform, freq_hz, varargin)
%   Layer:        L1
%   WP:           WP2
%   Claims:       C1, C2, C3
%   Status:       stub
%   Requirements: INS-01, INS-05
%
%   Ray-traced L1 echo generation over PHaRLAP's raytrace_3d.  NOT IMPLEMENTED.
%
%   Inputs
%     state      L0 truth, either a struct with fields alt_m, ne, nu, te,
%                b_field_t, lat_deg, lon_deg, epoch, source, or a path to the
%                CSV that TSS_IRI_PROFILES writes.  The struct fields are the
%                MATLAB image of tss.interfaces.IonosphereState.
%     platform   struct with fields alt_m, lat_deg, lon_deg,
%                velocity_ecef_mps, epoch.  The image of
%                tss.interfaces.PlatformState.
%     freq_hz    vector of sounding frequencies, Hz.
%
%     'Elevations'   elevation angles to launch, degrees.  Default 90, which
%                    is the vertical sounding case and the one that must agree
%                    with the analytic tier.
%     'Bearings'     ray bearings, degrees from north.  Default 0.
%     'Mode'         'O', 'X' or 'no-field'.  Default 'O'.
%     'NHops'        default 1.  A topside sounding is one hop up and back.
%     'Tolerance'    integrator tolerance passed to PHaRLAP.
%     'GridSpec'     struct describing the ionospheric grid PHaRLAP wants
%                    (start latitude and longitude, steps, counts, height
%                    start and step) plus the matching geomagnetic grid.
%     'OutFile'      .mat v7.3 path for the ray set.  When given, the result
%                    is saved through TSS_EXPORT_STRUCT so the Python side can
%                    read it with a sidecar rather than by guessing.
%
%   Outputs
%     rays       struct array, one entry per (frequency, mode, ray), with
%                fields freq_hz, mode, group_delay_s, amplitude_db,
%                absorption_db, doppler_hz, aoa_deg, reflection_alt_m and
%                flags.  These are the columns of tss.interfaces.EchoSet, so
%                the Python ray-traced L1 tier can build an EchoSet from the
%                exported .mat without a translation layer.
%     info       struct: pharlap_root, pharlap_version, grid spec actually
%                used, elapsed_s, and the path written.
%
%   Example
%     % Once WP2 lands:
%     st = tss_read_csv_table(profiles_csv);
%     rays = tss_raytrace_3d(st, plat, (1:0.1:12)*1e6, 'Elevations', 90);
%
%   WHY THIS IS A STUB AND NOT A WRAPPER
%
%   The signature above is the intended one for this repository, and it is
%   settled.  What is not settled is the exact argument list of PHaRLAP's own
%   raytrace_3d in 4.7.4: the order of the ionospheric grid arrays, the
%   geomagnetic field arrays and the two grid parameter vectors, and which of
%   the three return values carries which quantity.  That has to be read off
%   the installed toolkit's own help text and checked against a run, not
%   guessed at from memory.
%
%   A wrapper written from a guessed signature is worse than no wrapper.  It
%   would either fail at the call, which wastes the time of whoever runs it, or
%   silently pass arguments into the wrong slots and return a ray set that
%   looks plausible.  A ray set that looks plausible and is wrong would then be
%   checksummed, recorded in a manifest and compared against the analytic tier,
%   and the disagreement would be blamed on the analytic tier.
%
%   So this errors, clearly, and names the work package that closes it.
%
%   WHAT WP2 HAS TO DO
%
%   1. Read the raytrace_3d help in the installed PHaRLAP 4.7.4 and write the
%      argument list into this file as a comment before writing any code.
%   2. Build the ionospheric grid from an L0 state.  PHaRLAP wants electron
%      density on a latitude, longitude and height grid, plus a second grid
%      offset in time for the Doppler calculation, plus collision frequency,
%      plus the three geomagnetic field components on their own coarser grid.
%      TSS_IRI_PROFILES already produces density on a grid; the collision
%      frequency and the field are the missing pieces.
%   3. Run vertical incidence first and check group delay against
%      TSS_REFLECTION_HEIGHT and absorption against TSS_ABSORPTION_DB.  That
%      comparison is the V2 case for L1 and it is named as such in
%      docs/work-packages.md.
%   4. Only then add off-vertical rays and the O and X split.
%
%   See also TSS_IRI_PROFILES, TSS_EXPORT_STRUCT, TSS_REFLECTION_HEIGHT.

narginchk(3, inf);

error('tss:notImplemented', ...
    ['tss_raytrace_3d is a stub.  The ray-traced L1 tier is WP2 and has not ' ...
     'been written.\n' ...
     'The signature is settled; the PHaRLAP raytrace_3d argument list is not, ' ...
     'and it must be read off the installed toolkit rather than guessed.\n' ...
     'Until WP2 lands, use the analytic L1 tier on the Python side, or ' ...
     'tss_reflection_height and tss_absorption_db here.\n' ...
     'See the WHAT WP2 HAS TO DO block in this file and ' ...
     'docs/work-packages.md.']);
end
