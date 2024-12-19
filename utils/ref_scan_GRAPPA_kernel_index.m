% -_-_-_-_-_-_-_-_-_-_- ref_scan_GRAPPA_kernel_index _-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% generates source and target indices for GRAPPA kernel estimation from
% reference scans.
%
% Inputs:   params: scan parameters.
% ------
% 
% Outputs:  source_indx: source indicies.
% -------   target_indx: target indicies.
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [source_indx,target_indx] = ref_scan_GRAPPA_kernel_index(params)
% generate source indices for kernel estimation
source_indx_x=(repmat(1:params.grappa_win,[params.Ncol-params.grappa_win+1,1])+repmat((0:params.Ncol-params.grappa_win)',[1,params.grappa_win]))';
source_indx_x=repmat(source_indx_x,[params.grappa_line,1]);
source_indx_y=repmat((1:params.Rpe:(params.grappa_line-1)*params.Rpe+1)',[1,(params.NpeRef-(params.Rpe*params.grappa_line))/params.Rpe+1])+repmat(0:params.Rpe:params.NpeRef-(params.Rpe*params.grappa_line),[params.grappa_line,1]);
source_indx_y=(repmat(source_indx_y(:),[1,params.grappa_win]))';
source_indx_y=source_indx_y(:);
source_indx_y=reshape(source_indx_y,[params.grappa_line*params.grappa_win,(params.NpeRef-(params.Rpe*params.grappa_line))/params.Rpe+1]);
source_indx_y=repmat(source_indx_y,[size(source_indx_x,2),1]);
source_indx_y=source_indx_y(:);
source_indx_x=source_indx_x(:);
source_indx_x=repmat(source_indx_x,[(params.NpeRef-(params.Rpe*params.grappa_line))/params.Rpe+1,1]);
source_indx_y=source_indx_y(:);
source_indx=sub2ind([params.Ncol,params.NpeRef],source_indx_x,source_indx_y);

% generate target indices for kernel estimation
target_indx_x=(params.grappa_win-1)/2+1:params.Ncol-(params.grappa_win-1)/2;
target_indx_x=repmat(target_indx_x',[(params.NpeRef-(params.Rpe*params.grappa_line))/params.Rpe+1,1]);
target_indx_x=target_indx_x(:);
target_indx_y=(params.grappa_line/2-1)*params.Rpe+2:params.Rpe:params.NpeRef-params.Rpe*(params.grappa_line/2)+1;
target_indx_y=repmat(target_indx_y,[params.Ncol-params.grappa_win+1,1]);
target_indx_y=target_indx_y(:);
for k=0:params.Rpe-2
    target_indx(:,k+1)=sub2ind([params.Ncol,params.NpeRef],target_indx_x,target_indx_y+k);
end
end