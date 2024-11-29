% -_-_-_-_-_-_-_-_- imaging_scan_nyquist_ghost_corr -_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
%
% run Nyquist ghost (N/2) correction of EPI readout.
%
% Inputs:   twix: output of mapVBVD file.
% ------    v: reads the number of volume scanned.
%           params: scan parameters.
%
% Outputs:  data: corrected k-space data.
% -------
%
% Article:
% -------
%
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function data = imaging_scan_nyquist_ghost_corr(twix,v,params)
if(~params.silent)
    fprintf("    nyquist ghost correction....")
end
data=imaging_scan_phase_correction(twix,v);
data=padarray(data,[(params.grappa_win-1)/2,0,0,0],'both');
if(~params.silent)
    fprintf("Done!\n")
end
end