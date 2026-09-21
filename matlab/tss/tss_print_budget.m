function tss_print_budget(r)
%   Layer:        L2
%   WP:           WP3
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   Print a link budget result struct as a budget table.
%
%   Inputs
%     r   result struct from TSS_LINK_BUDGET.
%
%   Outputs
%     (none)  prints to stdout.
%
%   Example
%     tss_print_budget(tss_link_budget(5e6, tss_config()));
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_PRINT_BUDGET  Print a link budget result struct as a budget table.
rows = {
    'Frequency',                        sprintf('%.3f', r.f_mhz),        'MHz'
    'Transmit power',                   sprintf('%+.2f', r.tx_power_dbw), 'dBW'
    'Transmit antenna gain',            sprintf('%+.2f', r.tx_gain_dbi),  'dBi'
    'EIRP',                             sprintf('%+.2f', r.eirp_dbw),     'dBW'
    'Path loss (2R free space)',        sprintf('%+.2f', r.path_loss_db), 'dB'
    'Ionospheric absorption (2-way)',   sprintf('%+.2f', r.absorption_db),'dB'
    'Polarization loss',                sprintf('%+.2f', r.pol_loss_db),  'dB'
    'Other losses',                     sprintf('%+.2f', r.extra_loss_db),'dB'
    'Receive antenna gain',             sprintf('%+.2f', r.rx_gain_dbi),  'dBi'
    'Received signal power',            sprintf('%+.2f', r.rx_power_dbw), 'dBW'
    '',                                 '',                                ''
    'Sky temperature',                  sprintf('%.3e', r.sky_temp_k),    'K'
    'Receive efficiency',               sprintf('%+.2f', r.rx_efficiency_db), 'dB'
    'Antenna temperature',              sprintf('%.3e', r.ant_temp_k),    'K'
    'System temperature',               sprintf('%.3e', r.sys_temp_k),    'K'
    'Noise bandwidth (one range cell)', sprintf('%.2f',  r.noise_bw_hz),  'Hz'
    'Noise power',                      sprintf('%+.2f', r.noise_dbw),    'dBW'
    '',                                 '',                                ''
    'SNR, single sweep',                sprintf('%+.2f', r.snr_single_db),'dB'
    'Coherent integration gain',        sprintf('%+.2f', r.integration_db),'dB'
    'SNR',                              sprintf('%+.2f', r.snr_db),       'dB'
    'Status',                           r.status,                          ''
};
w = max(cellfun(@numel, rows(:,1)));
for k = 1:size(rows,1)
    if isempty(rows{k,1})
        fprintf('\n');
    else
        fprintf('%-*s  %12s  %s\n', w, rows{k,1}, rows{k,2}, rows{k,3});
    end
end
end
