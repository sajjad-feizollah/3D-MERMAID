% _-_-_-_-_-_-_-_-_-_-_-_-_- ref_scan_remove_phase _-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% phase correction of reference scans due to motion.
%
% Inputs:   ref_data: reference data after Nyquist ghost correction.
% ------    filt: triangle filter.
%           params: scan parameters.
%
% Outputs:  ref_kdata_cor: k-space of phase corrected projections of
%                          reference data.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function ref_kdata_cor = ref_scan_remove_phase(ref_data,filt,params)
ref_data=reshape(ref_data,[params.Ncol,params.NpeRef,params.Nproj,params.Ncoil]);
ref_im=zeros(size(ref_data));
ref_im_filt=zeros(size(ref_data));
for c=1:params.Ncoil
    ref_im(:,:,:,c)=fftshift(ifft2(fftshift(ref_data(:,:,:,c))));
    ref_im_filt(:,:,:,c)=fftshift(ifft2(fftshift(ref_data(:,:,:,c).*filt)));
end
ref_im_cor=abs(ref_im).*exp(1i*(angle(ref_im)-angle(ref_im_filt)));

fprintf("    resample into k-space....projection 000 of 000")
ref_kdata_cor=single(zeros([params.Ncol,params.NpeRef,params.Nproj,params.Ncoil]));
for p=1:params.Nproj
    fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d',p,params.Nproj)
    ref_kdata_cor(:,:,p,:)=fftshift(fft2(fftshift(squeeze(ref_im_cor(:,:,p,:)))));
end
fprintf("...Done!\n")

end