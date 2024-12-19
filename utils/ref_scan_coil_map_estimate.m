% -_-_-_-_-_-_-_-_-_- ref_scan_coil_map_estimate -_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% estimate coil sensitivity of each PE plane.
%
% Inputs:   ref_kdata_cor: k-space data of reference scans after phase
% ------                   correction.
%           traj: radial trajectory at each PE plane.
%           params: scan paramters.
% 
% Outputs:  sens_map: coil sensitivity map of each PE plane.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function sens_map = ref_scan_coil_map_estimate(ref_kdata_cor,traj,params)
sens_map=zeros(params.Nrecon,params.Npe,params.Nrecon,params.Ncoil);
fprintf("    estimate coil sensitivities....PE direction 000 of 000")
for p=1:params.Npe
    fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d',p,params.Npe)
    [~,sens]=bart('nlinv -a 32 -b 16  -S -d0 -i11 -g -x 32:32:1 -t',traj',reshape(ref_kdata_cor(:,p,:,:),1,params.Ncol*params.Nproj,1,params.Ncoil));
    sens=bart('fft -u 7', sens);
    sens=bart(strjoin({'resize -c 0 ',num2str(params.Nrecon*2),' 1 ',num2str(params.Nrecon*2)}), sens);
    sens=bart('fft -u -i 7', sens);
    sens=bart(strjoin({'resize -c 0 ',num2str(params.Nrecon),' 1 ',num2str(params.Nrecon)}), sens);
    sens_map(:,p,:,:)=bart('normalize 8',sens);
end
fprintf("...Done!\n")
disp("estimate coil sensitivity....Done!")
end