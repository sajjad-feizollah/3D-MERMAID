% -_-_-_-_-_-_-_-_-_-_- imaging_scan_generate_trajectory-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% generates radial trajectory in each PE plane from acquired spokes.
% 
% Inputs:   twix: output of mapVBVD that includes header.
% ------    ind_mc: indices of spokes to be removed from magnitude motion correction.
%           params: scan parameters.
%
% Outputs:  traj: radial trajectory
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function traj = imaging_scan_generate_trajectory(twix,ind_mc,params)
x=(-params.Ncol/2:1:params.Ncol/2-1)'/params.fovOS;
x=[x,zeros(params.Ncol,1)];
n=1;
proj_indx=twix.image.Par(1:twix.hdr.Dicom.EchoTrainLength:twix.hdr.Dicom.EchoTrainLength*params.Nproj)-1;
for p=1:params.Nproj
    theta=(pi/params.Nproj_total)*(proj_indx(p));
    traj(:,:,n)=x*[cos(theta),-sin(theta);sin(theta),cos(theta)];
    n=n+1;
end
if(params.bvalue>50)
    traj=traj(:,:,ind_mc);
end
traj=permute(traj,[1,3,2]);
traj=reshape(traj,[],2);
traj=single([traj,zeros(size(traj,1),1)]);
end