% _-_-_-_-_-_-_-_-_-_-_- imaging_scan_GRAPPA_recon _-_-_-_-_-_-_-_-_-_-_-_-
%
% Description:
% -----------
%
% GRAPPA reconstruction of each projection using calculated weights and
% indices.
%
% Inputs:   twix: output of mapVBVD file including header.
% ------    data: acquired k-space data after Nyquist ghost correction.
%           Wplane: GRAPPA weights of each projection.
%           source_indx: indices of acquired data.
%           target_indx: indices of missing k-space data.
%           filt: triangle filter.
%           params: scan parameters.
%
% Outputs:  im: reconstructred projections.
% -------   im_filt: reconstructed projections with triangle filter.
%
% Article:
% -------
%
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [im,im_filt] = imaging_scan_GRAPPA_recon(twix,data,Wplane,source_indx,target_indx,filt,params)
if(~params.silent)
    fprintf("    in-plane GRAPPA....projection 000 of 000")
end
im=single(zeros(params.Ncol,params.Npe,params.Nproj,params.Ncoil));
im_filt=single(zeros(params.Ncol,params.Npe,params.Nproj,params.Ncoil));
kdata=zeros(size(data,1),twix.image.Lin(params.Ntrain)+params.NpadPf+(params.grappa_line/2)*params.Rpe,params.Nproj,params.Ncoil);

kdata(:,twix.image.Lin(1:params.Ntrain)+params.NpadPf,:,:)=data;
kdata=reshape(kdata,[],params.Nproj,params.Ncoil);
for p=1:params.Nproj
    if(~params.silent)
        fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d',p,params.Nproj)
    end
    % fill missing k-space data
    source=reshape(squeeze(kdata(source_indx,p,:)),params.grappa_win*params.grappa_line,[],params.Ncoil);
    source=permute(source,[1,3,2]);
    source=reshape(source,[size(source,1)*size(source,2),size(source,3)]).';
    for k=0:params.Rpe-2
        kdata(target_indx(:,k+1),p,:)=source*squeeze(Wplane(:,:,p,k+1));
    end
    kdata_tmp=reshape(kdata(:,p,:),params.Ncol+params.grappa_win-1,[],params.Ncoil);
    kdata_tmp=kdata_tmp((params.grappa_win-1)/2+1:end-(params.grappa_win-1)/2,1:params.Npe,:);
    
    % reconstruct original and filtered projections
    switch params.doPOCS
        case true
            im(:,:,p,:)=imaging_scan_pocs_recon(kdata_tmp,params);
            im_filt(:,:,p,:)=imaging_scan_pocs_recon(kdata_tmp.*filt,params);
        case false
            im(:,:,p,:)=fftshift(ifft2(fftshift(kdata_tmp)));
            im_filt(:,:,p,:)=fftshift(ifft2(fftshift(kdata_tmp.*filt)));
    end
end
if(~params.silent)
    fprintf(" ...Done!\n")
end
end