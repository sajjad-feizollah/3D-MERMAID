% -_-_-_-_-_-_-_-_-_-_- ref_scan_GRAPPA_kernel_estimate -_-_-_-_-_-_-_-_-_-_
%
% Description: 
% -----------
% 
% estimates GRPPA kernel from reference scans and source and target
% indices.
%
% Inputs:   ref_data: reference k-space data after Nyquiest ghost
% ------              correction.
%           source_indx: indices for source samples
%           target_indx: indices for target samples
%           params: scan parameters.
%
% Outputs:  Wplane: GRAPPA weights of each projection.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function Wplane = ref_scan_GRAPPA_kernel_estimate(ref_data,source_indx,target_indx,params)
fprintf("    in-plane weight estimation....projection 000 of 000")
ref_data=reshape(ref_data,[params.Ncol*params.NpeRef,params.Nproj,params.Ncoil]);
Wplane=zeros(params.Ncoil*params.grappa_win*params.grappa_line,params.Ncoil,params.Nproj,params.Rpe-1);
for p=1:params.Nproj
    fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d',p,params.Nproj)
    source=reshape(squeeze(ref_data(source_indx,p,:)),[params.grappa_win*params.grappa_line,size(source_indx,1)/params.grappa_win/params.grappa_line,params.Ncoil]);
    source=permute(source,[1,3,2]);
    source=reshape(source,[size(source,1)*size(source,2),size(source,3)]).';
    for n=1:params.Rpe-1
        target=squeeze(ref_data(target_indx(:,n),p,:));
        Wplane(:,:,p,n)=pinv(source,params.grappa_tol) * target;

        % Wplane(:,:,p,n)=pinv([source;lambda*eye(size(source,2))]) *
        % [target;zeros(size(source,2),size(target,2))]; % uncomment for regularized kernel estimation

    end
end
fprintf("...Done!")
end