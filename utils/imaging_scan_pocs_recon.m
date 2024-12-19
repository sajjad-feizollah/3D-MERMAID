% -_-_-_-_-_-_-_-_-_-_-_- imaging_scan_pocs_recon -_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% POCS image reconstruction of each projection.
%
% Inputs:   kdata: k-space data after GRAPPA recon.
% ------    params: scan parametrs.
% 
% Outputs:  im_pocs: reconstructed image.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function im_pocs = imaging_scan_pocs_recon(kdata,params)

im_pocs = zeros(params.Ncol,params.Npe,params.Ncoil);
phs = zeros(params.Ncol,params.Npe,params.Ncoil);
phs(:,params.NpadPf+1:params.Npe-params.NpadPf,:) = ...
    kdata(:,params.NpadPf+1:params.Npe-params.NpadPf,:).*reshape(hann(params.Npe-2*params.NpadPf),1,[]);

phs = exp(1j*angle(fftshift(ifft2(fftshift(phs)))));
kdata=fftshift(ifft(fftshift(kdata,1),[],1),1);
for i=1:params.itrPOCS
    tmp = im_pocs.*phs;
    tmp = fftshift(fft(fftshift(tmp,2),[],2),2);
    tmp(:,params.NpadPf+1:end,:) = kdata(:,params.NpadPf+1:end,:);
    tmp = fftshift(ifft(fftshift(tmp,2),[],2),2).*conj(phs);
    tmp = max(real(tmp),0);
    im_pocs = tmp;
end
im_pocs = im_pocs.*phs;
end