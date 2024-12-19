% -_-_-_-_-_-_-_-_-_- ref_scan_nyquist_ghost_corr -_-_-__-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Nyquist ghost (N/2) correction of EPI readout for reference scans.
%
% Inputs:   twix: output of mapVBVD file.
% ------
% 
% Outputs:  ref_data: corrected k-space data.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
function ref_data = ref_scan_nyquist_ghost_corr(twix)
fprintf("    nyquist ghost correction of reference scans....")
data=ref_scan_phase_correction(twix);
data=reshape(data,[size(data,1),size(data,2)*size(data,3),size(data,4)]);
for k=1:size(data,2)
    ref_data(:,twix.refscan.Lin(k)-min(twix.refscan.Lin)+1,twix.refscan.Par(k),:)=data(:,k,:);
end
fprintf("Done!\n")
end