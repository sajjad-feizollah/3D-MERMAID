% -_-_-_-_-_-_-_-_- imaging_scan_magnitude_motion_corr -_-_-_-_-_-_-_-_-_-_
%
% Description: 
% -----------
% 
% magnitude motion correction of projections.
%
% Inputs:   kdata_cor: k-data of projections after phase correction.
% ------    im_cor: projections after phase correction.
%           params: scan parameters.
% Outputs:  kdata_cor: k-data of projections after magnitude correction.
% -------   ind_mc: indices of spokes to be reconstructed.
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [kdata_cor,ind_mc]=imaging_scan_magnitude_motion_corr(kdata_cor,im_cor,params)
if(params.doMC)
    if(params.bvalue>50)
        fprintf("    motion correction of projections....")
        img_cor=squeeze(sum(abs(im_cor).^2,4));
        diff_proj=diff(squeeze(sum(sum(img_cor,1),2)));
        ind_mc=diff_proj<=mean(abs(diff_proj)); % the difference is less than average means there was a void
        kdata_cor=kdata_cor(:,:,ind_mc,:);
        fprintf(" Done!\n")
    else
        ind_mc=1:size(kdata_cor,3);
    end
else
    ind_mc=1:size(kdata_cor,3);
end
end
