% -_-_-_-_-_-_-_-_-_- imaging_scan_GRAPPA_imaging_index -_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% 
%
% Inputs:   twix: output of mapVBVD file that includes indices.
% ------    params: scan parameters.
% 
% Outputs:  source_indx: indices of acquired scans from header.
% -------   target_indx: indices of missing k-space samples.
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [source_indx,target_indx] = imaging_scan_GRAPPA_imaging_index(twix,params)

% generate source indices of imaging scans
source_indx_x=(repmat(1:params.grappa_win,[params.Ncol,1])+repmat((0:params.Ncol-1)',[1,params.grappa_win]))';
source_indx_x=repmat(source_indx_x,[params.grappa_line,1]);
source_indx_y=twix.image.Lin(1)+params.NpadPf-(params.grappa_line/2-1)*params.Rpe:params.Rpe:twix.image.Lin(1)+params.NpadPf-(params.grappa_line/2-1)*params.Rpe+params.Rpe*(params.grappa_line-1);
source_indx_y=repmat(source_indx_y',[1,params.Ntrain])+repmat(params.Rpe*(0:params.Ntrain-1),[params.grappa_line,1]);
source_indx_y=(repmat(source_indx_y(:),[1,params.grappa_win]))';
source_indx_y=source_indx_y(:);
source_indx_y=reshape(source_indx_y,[params.grappa_line*params.grappa_win,params.Ntrain]);
source_indx_y=repmat(source_indx_y,[size(source_indx_x,2),1]);
source_indx_y=source_indx_y(:);
source_indx_x=source_indx_x(:);
source_indx_x=repmat(source_indx_x,[params.Ntrain,1]);
source_indx=sub2ind([params.Ncol+params.grappa_win-1,twix.image.Lin(params.Ntrain)+params.NpadPf+(params.grappa_line/2)*params.Rpe],source_indx_x,source_indx_y);

% generate target indices of imaging scans
target_indx_x=(1:params.Ncol)+(params.grappa_win-1)/2;
target_indx_x=repmat(target_indx_x',[params.Ntrain,1]);
target_indx_x=target_indx_x(:);
target_indx_y=twix.image.Lin(1:params.Ntrain)+params.NpadPf+1;
target_indx_y=repmat(target_indx_y,[params.Ncol,1]);
target_indx_y=target_indx_y(:);
target_indx=[];
for k=0:params.Rpe-2
    target_indx(:,k+1)=sub2ind([params.Ncol+(params.grappa_win-1),twix.image.Lin(params.Ntrain)+params.NpadPf+(params.grappa_line/2)*params.Rpe],target_indx_x,target_indx_y+k);
end

end