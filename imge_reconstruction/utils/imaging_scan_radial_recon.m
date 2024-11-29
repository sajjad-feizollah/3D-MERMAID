% -_-_-_-_-_-_-_-_-_-_-_- imaging_scan_radial_recon -_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% reconstructs a volume by reconstructing each PE plane.
%
% Inputs:   kdata_cor: final corrected k-space data.
% ------    traj: radial trajectory of each PE plane.
%           sens_map: coil sensitivity map of each PE plane.
%           params: scan parameters.
%
% Outputs:  img: reconstructed volume.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function img=imaging_scan_radial_recon(kdata_cor,traj,sens_map,params)
fprintf("    Reconstructs each PE plane....PE plane 000 of 000")
img=zeros(params.Npe,params.Nrecon,params.Nrecon);
for p=1:params.Npe
    fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d',p,params.Npe)
    switch params.recon
        case 'L2'   % L2 norm regularization
            img(params.Npe-p+1,:,:)=bart('pics -g -S -d0 -R Q:0.0001 -i25 -t',traj',reshape(kdata_cor(:,p,:,:),1,[],1,params.Ncoil),reshape(sens_map(:,p,:,:),params.Nrecon,params.Nrecon,1,params.Ncoil));
        case 'CS'   % CS with FISTA and 100 iterations
            img(params.Npe-p+1,:,:)=bart('pics -g -S -d0 -R W:3:0:0.001 -e --fista -i100 -t',traj',reshape(kdata_cor(:,p,:,:),1,[],1,params.Ncoil),reshape(sens_map(:,p,:,:),params.Nrecon,params.Nrecon,1,params.Ncoil));
        otherwise   % regular SENSE
            img(params.Npe-p+1,:,:)=bart('pics -g -S -d0 -i10 -t',traj',reshape(kdata_cor(:,p,:,:),1,[],1,params.Ncoil),reshape(sens_map(:,p,:,:),params.Nrecon,params.Nrecon,1,params.Ncoil));
    end
end

fprintf(" Done!\n")
end