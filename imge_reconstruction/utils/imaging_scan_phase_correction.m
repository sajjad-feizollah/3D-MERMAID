% _-_-_-_-_-_-_-_-_-_- imaging_scan_phase_correction -_-_-_-_-_-_-_-_-_-_-_
%
% Description: 
% -----------
% 
% corrects Nyquist ghost (N/2) correction of EPI readout using the method
% explained in https://patents.google.com/patent/US6043651A/en
%
% Inputs:   twix: output of mapVBVD including phase correction scans.
% ------    v: number of volume to be reconstructed.
% 
% Outputs:  data: corrected k-spce data.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function data=imaging_scan_phase_correction(twix,v)

TrainLength = twix.hdr.Dicom.EchoTrainLength;
NCol= twix.image.NCol;
Nproj=twix.image.NAcq/twix.hdr.Dicom.EchoTrainLength/twix.image.NRep;
Nread=Nproj*TrainLength;
phaseCor=twix.phasecor.unsorted((v-1)*3*Nproj+1:3*Nproj*v);

phaseCor=permute(phaseCor,[1,3,2]);
phaseCor=reshape(phaseCor,[size(phaseCor,1),3,size(phaseCor,2)/3,size(phaseCor,3)]);

if(rem(twix.hdr.Dicom.EchoTrainLength,2))
    FTphaseCorEven=ifftshift(ifft(ifftshift(squeeze((phaseCor(:,1,:,:)+phaseCor(:,3,:,:))/2)),[],1));
    FTphaseCorOdd=ifftshift(ifft(ifftshift(squeeze(phaseCor(:,2,:,:))),[],1));
else
    FTphaseCorOdd=ifftshift(ifft(ifftshift(squeeze((phaseCor(:,1,:,:)+phaseCor(:,3,:,:))/2)),[],1));
    FTphaseCorEven=ifftshift(ifft(ifftshift(squeeze(phaseCor(:,2,:,:))),[],1));
end
dPhiOdd=squeeze(angle(sum(FTphaseCorOdd(2:end,:,:).*conj(FTphaseCorOdd(1:end-1,:,:)),1)));
dPhiEven=squeeze(angle(sum(FTphaseCorEven(2:end,:,:).*conj(FTphaseCorEven(1:end-1,:,:)),1)));
dPhiOdd=reshape(dPhiOdd,[1,size(dPhiOdd,1),size(dPhiOdd,2)]);
dPhiEven=reshape(dPhiEven,[1,size(dPhiEven,1),size(dPhiEven,2)]);

FTphaseCorOdd=FTphaseCorOdd.*exp(-1i*(-NCol/2:NCol/2-1)'.*dPhiOdd);
FTphaseCorEven=FTphaseCorEven.*exp(-1i*(-NCol/2:NCol/2-1)'.*dPhiEven);
phi0=squeeze(angle(sum(FTphaseCorEven.*conj(FTphaseCorOdd),3)));

clear phaseCor FTphaseCorEven FTphaseCorOdd;

data=twix.image.unsorted((v-1)*Nread+1:Nread*v);
data=permute(data,[1,3,2]);
data=reshape(data,[NCol,TrainLength,size(data,2)/TrainLength,size(data,3)]);

data=ifftshift(ifft(ifftshift(data),[],1));

dPhiEven=reshape(dPhiEven,[1,1,size(dPhiEven,2),size(dPhiEven,3)]);
dPhiOdd=reshape(dPhiOdd,[1,1,size(dPhiOdd,2),size(dPhiOdd,3)]);
phi0=reshape(phi0,[size(phi0,1),1,size(phi0,2)]);

data(:,1:2:end,:,:)=data(:,1:2:end,:,:).*exp(-1i*(-NCol/2:NCol/2-1)'.*dPhiOdd);
data(:,2:2:end,:,:)=data(:,2:2:end,:,:).*exp(-1i*(-NCol/2:NCol/2-1)'.*dPhiEven);
data(:,1:2:end,:,:)=data(:,1:2:end,:,:).*exp(1i*phi0);

data=fftshift(fft(fftshift(data),[],1));
end