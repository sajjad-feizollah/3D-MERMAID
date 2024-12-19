% -_-_-_-_-_-_-_-_-_-_- ref_scan_generate_trajectory _-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% generate radial trajectory of each PE plane.
%
% Inputs:   params: scan parameters.
% ------
% 
% Outputs:  traj: radial trajectory.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function traj = ref_scan_generate_trajectory(params)
x=(-params.Ncol/2:1:params.Ncol/2-1)'/params.fovOS;
x=[x,zeros(params.Ncol,1)];
traj=zeros(params.Ncol,2,params.Nproj);
for p=1:params.Nproj
    theta=(pi/params.Nproj_total)*(p-1);
    traj(:,:,p)=x*[cos(theta),-sin(theta);sin(theta),cos(theta)];
end
traj=permute(traj,[1,3,2]);
traj=reshape(traj,[],2);
traj=single([traj,zeros(size(traj,1),1)]);
end